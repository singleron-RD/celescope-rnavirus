import os
import shutil

import numpy as np
import pandas as pd
import scipy
from celescope.tools import utils
from celescope.tools.step import Step, s_common
import subprocess
import sys
from celescope.tools.matrix import CountMatrix


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

        self.kb_whitelist_file = f"{self.outdir}/kb_whitelist.txt"
        remove_underscore(args.barcodes_file, self.kb_whitelist_file)

        # outfile
        self.filtered = f"{self.outdir}/transcript.filtered"
        self.outs = [self.filtered]

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

    def run_kb(self):
        fq_line = " ".join(
            [
                f"{fq1} {fq2}"
                for fq1, fq2 in zip(self.args.fq1.split(","), self.args.fq2.split(","))
            ]
        )
        cmd = (
            f"kb count -i {self.args.kbDir}/index.idx -g {self.args.kbDir}/t2g.txt "
            f"-x 0,0,9,0,9,18,0,18,27:0,27,35:1,0,0 "
            f"-w {self.kb_whitelist_file} "
            f"-t {self.args.thread} "
            f"-o {self.outdir} --tcc --quant-umis --overwrite "
            f"{fq_line} "
        )
        sys.stderr.write(cmd + "\n")
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def get_transcript_gene_dict(self):
        df = pd.read_csv(f"{self.args.kbDir}/t2g.txt", sep="\t", header=None)
        transcript_gene_dict = dict(zip(df[0], df[2]))
        return transcript_gene_dict

    @utils.add_log
    def get_filtered_matrix(self):
        os.makedirs(self.filtered, exist_ok=True)
        add_underscore(
            f"{self.outdir}/counts_unfiltered/cells_x_tcc.barcodes.txt",
            f"{self.filtered}/barcodes.tsv",
        )
        features_file = f"{self.filtered}/features.tsv"
        shutil.copy(
            f"{self.outdir}/quant_unfiltered/transcripts.txt",
            features_file,
        )
        transcript_gene_dict = self.get_transcript_gene_dict()
        features_df = pd.read_csv(features_file, sep="\t", header=None)
        features_df[1] = [
            transcript_gene_dict[transcript] for transcript in features_df[0]
        ]
        features_df[2] = "Gene Expression"
        features_df.to_csv(features_file, sep="\t", header=False, index=False)
        matrix = scipy.io.mmread(f"{self.outdir}/quant_unfiltered/matrix.abundance.mtx")
        matrix = matrix.transpose()
        with open(f"{self.filtered}/matrix.mtx", "wb") as f:
            scipy.io.mmwrite(f, matrix)
        cmd = f"gzip {self.filtered}/barcodes.tsv {self.filtered}/features.tsv {self.filtered}/matrix.mtx "
        sys.stderr.write(cmd + "\n")
        subprocess.check_call(cmd, shell=True)

        count_matrix = CountMatrix.from_matrix_dir(self.filtered)
        bc_transcriptNum, total_transcript = count_matrix.get_bc_geneNum()
        self.add_metric(
            name="Total Transcripts",
            value=total_transcript,
            help_info="The number of transcripts with at least one UMI count in any cell.",
        )
        self.add_metric(
            name="Median Transcripts per Cell",
            value=round(np.median(list(bc_transcriptNum.values())), 2),
            help_info="median number of transcripts per cell.",
        )

    @utils.add_log
    def run(self):
        if not self.args.kbDir:
            sys.stderr.write("--kbDir not provided. Skip kb_python")
            return
        self.run_kb()
        self.get_filtered_matrix()


@utils.add_log
def kb_python(args):
    with Kb_python(args, display_title="Transcripts") as runner:
        runner.run()


def get_opts_kb_python(parser, sub_program):
    parser.add_argument(
        "--kbDir",
        help="kb reference directory path. Must contain index.idx and t2g.txt.",
    )
    parser.add_argument(
        "--chemistry",
        default="mobiu-1",
        choices=["mobiu-1", "mobiu-2", "mobiu-3"],
        help="chemistry version",
    )
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
        parser.add_argument("--barcodes_file", required=True)
        parser = s_common(parser)

    return parser
