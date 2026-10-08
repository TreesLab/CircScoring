#! /usr/bin/env python

"""
A python script to get the sequences consist of the flanking +-Nbp
 from the circRNA junction.

Note:
  The `bedtools` is required.
  We use `bedtools getfasta` to get the needed sequences from the genome.

Usage:
  get_flanking_seq.py [Genome] [CircRNAs] [Output] ([N])

  Genome:
    The reference genome file.

  CircRNAs:
    CircRNAs data in the 4-columns format with (chr, pos1, pos2, strand).

  Output:
    The output file in FASTA format.

  N:
    (Optional)
    Get the sequences consist of the flanking +-Nbp from the circRNA junction. (default: 100)
"""

import sys
import io
import re
import tempfile as tp
import subprocess as sp
import logging
from pathlib import Path


logging.basicConfig(
    format="{asctime} - {message}",
    level=logging.INFO,
    style='{'
)


class Bed:
    BED_TITLE = (
        'chr',
        'start',
        'end',
        'name',
        'score',
        'strand',
        'thickStart',
        'thickEnd',
        'itemRGB',
        'blockCount',
        'blockSizes',
        'blockStarts'
    )

    def __init__(self, regions, name='.'):
        self.regions = regions
        self.name = name

        self._parse_regions(regions)

    def _parse_regions(self, regions):
        self._parse_region(regions[0])
        for region in regions[1:]:
            self._add_block(region)

        if self.strand == '-':
            self._reverse_region()

    def _parse_region(self, region):
        self.chrom = region[0]
        self.start = int(region[1]) - 1
        self.end = int(region[2])
        self.strand = region[3]

        self.block_count = 1
        self.block_sizes = [self.end - self.start]
        self.block_starts = [0]

    def _add_block(self, region):
        other = Bed([region])

        assert self.strand == other.strand

        self.block_count += other.block_count
        self.block_sizes += other.block_sizes

        all_starts = [
            start + self.start for start in self.block_starts
        ] + [other.start]

        if other.start < self.start:
            self.start = other.start

        self.block_starts = [start - self.start for start in all_starts]

        if other.end > self.end:
            self.end = other.end

    def _reverse_region(self):
        self.block_sizes = list(reversed(self.block_sizes))
        self.block_starts = list(reversed(self.block_starts))

    def get_data(self, all_fields=False):
        fields = [
            self.chrom,
            self.start,
            self.end,
            self.name,
            '.',
            self.strand
        ]

        if all_fields:
            fields += [
                self.start,
                self.end,
                0,
                self.block_count,
                self._list_to_str(self.block_sizes, sep=','),
                self._list_to_str(self.block_starts, sep=',')
            ]

        return fields

    def to_string(self, all_fields=False):
        data = self.get_data(all_fields=all_fields)
        bed_txt = self._list_to_str(data)
        return bed_txt

    @staticmethod
    def _list_to_str(list_, sep='\t', end=''):
        with io.StringIO() as tmp:
            print(*list_, sep=sep, end=end, file=tmp)
            return tmp.getvalue()


def get_fasta(bed_file, ref_file, use_blocks=False, bedtools_bin='bedtools'):
    cmd = [bedtools_bin, 'getfasta']
    cmd += ['-fi', ref_file]
    cmd += ['-bed', bed_file]
    cmd += ['-fo', '/dev/stdout']
    cmd += ['-nameOnly', '-s']

    if use_blocks:
        cmd += ['-split']

    res = sp.run(cmd, stdout=sp.PIPE, encoding='utf-8')

    return res.stdout


class CircRNAs:
    def __init__(self, chr_, pos1, pos2, strand):
        self.chr_ = chr_
        self.pos1 = int(pos1)
        self.pos2 = int(pos2)
        self.strand = strand

        self.get_donor_acceptor()

    def get_donor_acceptor(self):
        pos1, pos2 = sorted([self.pos1, self.pos2])

        if self.strand == '+':
            self.donor_site = pos2
            self.acceptor_site = pos1
        elif self.strand == '-':
            self.donor_site = pos1
            self.acceptor_site = pos2

    def get_flanking_region(self, N):
        donor_flanking = self._get_flanking_region_of_donor(N)
        acceptor_flanking = self._get_flanking_region_of_acceptor(N)
        return [donor_flanking, acceptor_flanking]

    def _get_flanking_region_of_donor(self, N):
        if self.strand == '+':
            donor_start = self.donor_site - (N - 1)
            donor_end = self.donor_site

        elif self.strand == "-":
            donor_start = self.donor_site
            donor_end = self.donor_site + (N - 1)

        donor_flanking = [
            self.chr_,
            donor_start,
            donor_end,
            self.strand
        ]

        return donor_flanking

    def _get_flanking_region_of_acceptor(self, N):
        if self.strand == '+':
            acceptor_start = self.acceptor_site
            acceptor_end = self.acceptor_site + (N - 1)

        elif self.strand == "-":
            acceptor_start = self.acceptor_site - (N - 1)
            acceptor_end = self.acceptor_site

        acceptor_flanking = [
            self.chr_,
            acceptor_start,
            acceptor_end,
            self.strand
        ]

        return acceptor_flanking

    @property
    def id(self):
        return f"{self.chr_}:{self.donor_site}|{self.acceptor_site}({self.strand})"


def get_flanking_seq(genome_file, circRNAs_file, out_file, N=100):
    logging.info("Start to get the flanking sequences.")

    output_path = Path(out_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_bed = tp.NamedTemporaryFile(
        prefix='flanking_region.',
        suffix='.bed',
        dir=output_path.parent
    )

    with open(circRNAs_file) as circRNA_data, \
            open(tmp_bed.name, 'w') as bed_out:
        for line in circRNA_data:
            circRNA = CircRNAs(*line.rstrip('\n').split('\t'))
            flanking_region = circRNA.get_flanking_region(int(N))
            region_bed = Bed(flanking_region, name=circRNA.id)
            print(region_bed.to_string(all_fields=True), file=bed_out)

    fasta_data = get_fasta(tmp_bed.name, genome_file, use_blocks=True)

    with output_path.open('w') as out:
        out.write(
            re.sub(r'(?<=\([+-]\))\([+-]\)$', '', fasta_data, flags=re.M)
        )

    tmp_bed.close()

    logging.info("Process completed!")


def print_usage():
    print(__doc__, file=sys.stderr)


if __name__ == "__main__":

    if len(sys.argv) < 4:
        print_usage()
        exit(1)

    get_flanking_seq(*sys.argv[1:])
