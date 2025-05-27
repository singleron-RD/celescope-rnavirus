"""
convert BAM to R1 and R2
"""

from celescope.tools import utils, parse_chemistry
from celescope.tools.step import Step, s_common
import subprocess
import gzip
import pysam


class Convert(Step):
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        self.r1_file = f"{self.outdir}/{self.sample}_R1.fq.gz"
        self.r2_file = f"{self.outdir}/{self.sample}_R2.fq.gz"
        self.fq1_list = args.fq1.split(",")
        self.fq2_list = args.fq2.split(",")
        self.chemistry = parse_chemistry.get_chemistry("rna", "auto", self.fq1_list)

    @utils.add_log
    def cat_r1_files(self):
        fq1_str = " ".join(self.fq1_list)
        cmd1 = f"cat {fq1_str} > {self.r1_file}"
        subprocess.check_call(cmd1, shell=True)

    @utils.add_log
    def cat_r2_files(self):
        fq2_str = " ".join(self.fq2_list)
        cmd2 = f"cat {fq2_str} > {self.r2_file}"
        subprocess.check_call(cmd2, shell=True)

    @utils.add_log
    def convert_v3(self):
        auto_runner = parse_chemistry.AutoRNA(self.fq1_list)
        with gzip.open(self.r1_file, mode="wt", compresslevel=1) as r1:
            for fq1_file in self.fq1_list:
                with pysam.FastqFile(fq1_file) as fq1:
                    for read in fq1:
                        seq = read.sequence
                        qual = read.quality
                        offset = auto_runner.v3_offset(seq)
                        if offset != -1:
                            seq = seq[offset:]
                            qual = qual[offset:]
                        r1.write(utils.fastq_line(read.name, seq, qual))

    def run(self):
        self.cat_r2_files()
        if self.chemistry != "GEXSCOPE-V3":
            self.cat_r1_files()
        else:
            self.convert_v3()


@utils.add_log
def convert(args):
    with Convert(args) as runner:
        runner.run()


def get_opts_convert(parser, sub_program=True):
    if sub_program:
        parser.add_argument(
            "--fq1",
            required=True,
        )
        parser.add_argument(
            "--fq2",
            required=True,
        )
        parser = s_common(parser)

    return parser
