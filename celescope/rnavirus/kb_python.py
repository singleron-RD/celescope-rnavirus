import os
import shutil

import pandas as pd
import scipy
from celescope.tools import utils
from celescope.tools.step import Step, s_common
import subprocess
import sys
from celescope.tools import parse_chemistry, matrix


def add_underscore(barcodes_file, outfile):
    with utils.generic_open(barcodes_file, "r") as f:
        lines = f.readlines()
        barcode_length = len(lines[0].strip()) // 3
        with utils.generic_open(outfile, "wt") as f:
            for line in lines:
                line = line.strip()
                line = "_".join(
                    [
                        line[i : i + barcode_length]
                        for i in range(0, len(line), barcode_length)
                    ]
                )
                f.write(line + "\n")


def remove_underscore(barcodes_file, outfile):
    with utils.generic_open(barcodes_file, "rt") as f:
        lines = f.readlines()
        with utils.generic_open(outfile, "wt") as f:
            for line in lines:
                line = line.strip()
                line = line.replace("_", "")
                f.write(line + "\n")


class Kb_python(Step):
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        self.chemistry = self.get_slot_key(
            slot="metrics", step_name="sample", key="Chemistry"
        )
        self.pattern_dict, _bc = parse_chemistry.get_pattern_dict_and_bc(self.chemistry)

        self.kb_whitelist_file = f"{self.outdir}/kb_whitelist.txt"
        remove_underscore(args.barcodes_file, self.kb_whitelist_file)

        # mapping
        self.mapping_df = pd.read_csv(
            f"{args.kbDir}/ID_to_taxonomy_mapping.csv", na_values=["."]
        )
        self.rep_df = self.mapping_df[
            self.mapping_df["ID"] == self.mapping_df["rep_ID"]
        ]
        self.rep_df = self.rep_df.drop("ID", axis=1)
        self.rep_df = self.rep_df.dropna(
            subset=["phylum", "class", "order", "family", "genus", "species"], how="all"
        )

        # outfile
        self.filtered = f"{self.outdir}/virus.filtered"
        self.positive_cell_csv = f"{self.outdir}/positive_cell.csv"
        self.outs = [self.filtered, self.positive_cell_csv]

    def generate_kb_whitelist(self, bc_list, kb_whitelist_file):
        bc_segments = [open(bc_file, "r").read().splitlines() for bc_file in bc_list]
        # remove underscore
        bcs = [
            "".join([bc1, bc2, bc3])
            for bc1 in bc_segments[0]
            for bc2 in bc_segments[1]
            for bc3 in bc_segments[2]
        ]
        with open(self.kb_whitelist_file, "wt") as f:
            for bc in bcs:
                f.write(bc + "\n")

    def get_kb_pattern_args(self):
        cb_pos = ",".join([f"0,{x.start},{x.stop}" for x in self.pattern_dict["C"]])
        umi_pos = (
            f"0,{self.pattern_dict['U'][0].start},{self.pattern_dict['U'][0].stop}"
        )
        return ":".join([cb_pos, umi_pos, "1,0,0"])

    def run_kb(self):
        pattern_args = self.get_kb_pattern_args()
        cmd = (
            f"kb count --aa -i {self.args.kbDir}/index.idx -g {self.args.kbDir}/palmdb_clustered_t2g.txt "
            f"-x {pattern_args} "
            f"-w {self.kb_whitelist_file} "
            f"-t {self.args.thread} "
            f"-o {self.outdir}  --overwrite "
            f"{self.args.fq1} {self.args.fq2}"
        )
        sys.stderr.write(cmd + "\n")
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def get_filtered_matrix(self):
        os.makedirs(self.filtered, exist_ok=True)
        add_underscore(
            f"{self.outdir}/counts_unfiltered/cells_x_genes.barcodes.txt",
            f"{self.filtered}/barcodes.tsv",
        )
        features_file = f"{self.filtered}/features.tsv"
        shutil.copy(
            f"{self.outdir}/counts_unfiltered/cells_x_genes.genes.txt",
            features_file,
        )
        features_df = pd.read_csv(features_file, sep="\t", header=None)
        features_df[1] = features_df[0]
        features_df[2] = "RNA Virus"
        features_df.to_csv(features_file, sep="\t", header=False, index=False)
        matrix = scipy.io.mmread(f"{self.outdir}/counts_unfiltered/cells_x_genes.mtx")
        matrix = matrix.transpose()
        with open(f"{self.filtered}/matrix.mtx", "wb") as f:
            scipy.io.mmwrite(f, matrix)
        cmd = f"gzip {self.filtered}/barcodes.tsv {self.filtered}/features.tsv {self.filtered}/matrix.mtx "
        sys.stderr.write(cmd + "\n")
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def get_positive_cell_count(self):
        filtered = matrix.CountMatrix.from_matrix_dir(self.filtered)
        df = filtered.to_df(gene_id=True)
        sum_df = df.sum(axis=1)
        sum_df = sum_df[sum_df > 0]
        sum_df.name = "positive_cell"
        sum_df = sum_df.to_frame()
        sum_df["positive_cell"] = sum_df["positive_cell"].astype(int)

        merge_df = sum_df.merge(
            self.rep_df, left_index=True, right_on="rep_ID", how="inner"
        )
        merge_df.to_csv(self.positive_cell_csv)

        table_dict = self.get_table_dict(
            title="Positive Cell Count", table_id="rnavirus", df_table=merge_df
        )
        self.add_data(table_dict=table_dict)

    @utils.add_log
    def run(self):
        self.run_kb()
        self.get_filtered_matrix()
        self.get_positive_cell_count()


@utils.add_log
def kb_python(args):
    with Kb_python(args, display_title="RNA Virus") as runner:
        runner.run()


def get_opts_kb_python(parser, sub_program):
    parser.add_argument(
        "--kbDir",
        help="kb reference directory path. Must contain index.idx, palmdb_clustered_t2g.txt and ID_to_taxonomy_mapping.csv.",
        required=True,
    )
    parser.add_argument("--barcodes_file", required=True)
    if sub_program:
        parser.add_argument(
            "--fq1",
            help="R1 fastq file.",
            required=True,
        )
        parser.add_argument(
            "--fq2",
            help="R2 fastq file.",
            required=True,
        )
        parser = s_common(parser)

    return parser
