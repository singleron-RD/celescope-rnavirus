"""
convert BAM to R1 and R2
"""

from celescope.tools import utils, parse_chemistry
from celescope.tools.step import Step, s_common
import subprocess
import pysam


class Convert(Step):
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        self.bc_set = set(utils.one_col_to_list(self.args.bc))

    def run(self):
        with pysam.AlignmentFile(self.args.bam) as bam:
            for read in bam:
                cb = read.get_tag("CB")
                ub = read.get_tag("UB")
                if cb not in self.bc_set:
                    continue
                seq = read.seq


@utils.add_log
def convert(args):
    with Convert(args) as runner:
        runner.run()


def get_opts_convert(parser, sub_program=True):
    parser.add_argument(
        "--bam",
        required=True,

    )
    parser.add_argument(
        "--bc",
        help='scRNA-Seq barcode.tsv.gz file path.',
        required=True,

    )
    if sub_program:
        parser = s_common(parser)

    return parser
