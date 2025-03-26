import os
import shutil

import numpy as np
import scipy
from celescope.tools import utils
from celescope.tools.step import Step, s_common
import subprocess
import sys
from celescope.tools.parse_chemistry import get_pattern_dict_and_bc
from celescope.tools.matrix import CountMatrix


def add_underscore(barcodes_file, outfile):
    """
    add underscore to barcodes
    """

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


class Kb_python(Step):
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        _pattern_dict, bc_list = get_pattern_dict_and_bc(self.args.chemistry)
        self.kb_whitelist_file = f"{self.outdir}/kb_whitelist.txt"
        self.generate_kb_whitelist(bc_list, self.kb_whitelist_file)

        # outfile
        self.raw = f"{self.outdir}/transcript.raw"
        self.filtered = f"{self.outdir}/transcript.filtered"
        self.outs = [self.raw, self.filtered]

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
    def get_raw_matrix(self):
        os.makedirs(self.raw, exist_ok=True)
        add_underscore(
            f"{self.outdir}/counts_unfiltered/cells_x_tcc.barcodes.txt",
            f"{self.raw}/barcodes.tsv",
        )
        shutil.copy(
            f"{self.outdir}/quant_unfiltered/transcripts.txt",
            f"{self.raw}/features.tsv",
        )
        matrix = scipy.io.mmread(f"{self.outdir}/quant_unfiltered/matrix.abundance.mtx")
        matrix = matrix.transpose()
        with open(f"{self.raw}/matrix.mtx", "wb") as f:
            scipy.io.mmwrite(f, matrix)
        cmd = f"gzip {self.raw}/barcodes.tsv {self.raw}/features.tsv {self.raw}/matrix.mtx "
        sys.stderr.write(cmd + "\n")
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def get_filtered_matrix(self):
        raw_matrix = CountMatrix.from_matrix_dir(self.raw)
        barcodes, _ = utils.get_barcode_from_matrix_dir(self.args.gene_filtered)
        filtered_matrix = raw_matrix.slice_matrix_bc(barcodes)
        filtered_matrix.to_matrix_dir(self.filtered)

        bc_transcriptNum, total_transcript = filtered_matrix.get_bc_geneNum()
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
        # self.run_kb()
        self.get_raw_matrix()
        self.get_filtered_matrix()


@utils.add_log
def kb_python(args):
    with Kb_python(args, display_title="Transcripts") as runner:
        runner.run()


def get_opts_kb_python(parser, sub_program):
    parser.add_argument(
        "--kbDir",
        help="kb reference directory path. Must contain index.idx and t2g.txt.",
        required=True,
    )
    parser.add_argument(
        "--chemistry",
        required=True,
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
        parser.add_argument("--gene_filtered", required=True)
        parser = s_common(parser)

    return parser
