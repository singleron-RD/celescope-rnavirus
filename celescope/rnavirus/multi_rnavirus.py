from celescope.rnavirus.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_rnavirus(Multi):
    def convert(self, sample):
        step = "convert"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        fq1 = ",".join(arr["fq1"])
        fq2 = ",".join(arr["fq2"])
        cmd = f"{cmd_line} " f"--fq1 {fq1} --fq2 {fq2} "
        self.process_cmd(cmd, step, sample)

    def kb_python(self, sample):
        step = "kb_python"
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f"{cmd_line} "
            f"--fq1 {self.outdir_dic[sample]['convert']}/{sample}_R1.fq.gz "
            f"--fq2 {self.outdir_dic[sample]['convert']}/{sample}_R2.fq.gz "
        )
        self.process_cmd(cmd, step, sample, m=10, x=self.args.thread)


def main():
    multi = Multi_rnavirus(__ASSAY__)
    multi.run()


if __name__ == "__main__":
    main()
