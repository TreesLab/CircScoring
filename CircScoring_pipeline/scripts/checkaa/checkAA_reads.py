#! /usr/bin/env python

"""
A python script to check for reads if there exist colinear alignments
 or multiple alignments.

External tools:
  - blat and mp_blat.py, or pblat

Usage:
  checkAA_reads.py -rg [Genome] -ro [OtherRefs] -p [NumProc] [Reads]
"""


import argparse
import glob
import tempfile as tp
import subprocess as sp
import csv
import logging
import os
import shlex
import textwrap
from itertools import groupby, chain
from functools import partial
from collections import namedtuple, deque, OrderedDict

logging.basicConfig(
    format="{asctime} - {message}",
    level=logging.INFO,
    style='{'
)


WEB_BLAT_OPTIONS = "-tileSize=9 -stepSize=9 -repMatch=32768"
COLINEAR_THRESHOLD = 0.8
END_DISTANCE = 10
CHIMERA_MIN_MATCHES = 30
CHIMERA_MIN_DIFF = 3
RESULT_HEADER = [
    'circRNA_id',
    'with an alternative co-linear explanation',
    'with multiple_hits',
    'alignment ambiguity '
    '(with an alternative co-linear explanation or multiple hits)'
]

FastaRecord = namedtuple('FastaRecord', ('id', 'sequence'))


class Blat:
    def __init__(self, work_dir='.', blat_bin='blat', blat_opts="",
                 num_proc=1, mp_blat_bin='mp_blat.py'):
        self.work_dir = work_dir
        self.blat_bin = blat_bin
        self.blat_opts = blat_opts
        self.num_proc = num_proc
        self.mp_blat_bin = mp_blat_bin

    def _run(self, ref_file, fa_file):
        fd, out_file = tp.mkstemp(dir=self.work_dir, suffix='.psl')
        os.close(fd)

        try:
            cmd = self._generate_cmd(ref_file, fa_file, out_file)
            sp.run(cmd, check=True)
            return out_file
        except Exception:
            if os.path.exists(out_file):
                os.remove(out_file)
            raise

    def _generate_cmd(self, ref_file, fa_file, out_file):
        cmd = [self.mp_blat_bin, ref_file, fa_file, out_file]
        cmd += ['-p', str(self.num_proc)]
        cmd += ['--tmp_path', self.work_dir]
        cmd += ['--blat_bin', self.blat_bin]
        cmd += ['--blat_options', self.blat_opts]
        return cmd

    def __call__(self, ref_file, fa_file):
        return self._run(ref_file, fa_file)


class Pblat:
    def __init__(self, work_dir='.', pblat_bin='pblat', blat_opts="",
                 num_proc=1):
        self.work_dir = work_dir
        self.pblat_bin = pblat_bin
        self.blat_opts = blat_opts
        self.num_proc = num_proc

    def _run(self, ref_file, fa_file):
        fd, out_file = tp.mkstemp(dir=self.work_dir, suffix='.psl')
        os.close(fd)

        try:
            cmd = self._generate_cmd(ref_file, fa_file, out_file)
            logging.info(
                f'Starting pblat with {self.num_proc} thread'
                f'{"" if self.num_proc == 1 else "s"} on {fa_file}'
            )
            sp.run(cmd, check=True)
            logging.info(f'pblat completed. ({fa_file})')
            return out_file
        except Exception:
            for path in [out_file] + glob.glob(out_file + '.tmp.*'):
                if os.path.exists(path):
                    os.remove(path)
            raise

    def _generate_cmd(self, ref_file, fa_file, out_file):
        cmd = [self.pblat_bin, f'-threads={self.num_proc}']
        cmd += shlex.split(self.blat_opts)
        cmd += [ref_file, fa_file, out_file]
        return cmd

    def __call__(self, ref_file, fa_file):
        return self._run(ref_file, fa_file)


class BlatPsl:
    PSL_TITLE = (
        'matches',
        'misMatches',
        'repMatches',
        'nCount',
        'qNumInsert',
        'qBaseInsert',
        'tNumInsert',
        'tBaseInsert',
        'strand',
        'qName',
        'qSize',
        'qStart',
        'qEnd',
        'tName',
        'tSize',
        'tStart',
        'tEnd',
        'blockCount',
        'blockSizes',
        'qStarts',
        'tStarts'
    )
    _psl = namedtuple('Psl', PSL_TITLE)

    def __init__(self, psl_text):
        self.psl_data = self._parse(psl_text)

    @classmethod
    def _parse(cls, psl_text):
        psl_data = []

        for line in cls.iter_text(psl_text):
            psl_data.append(line)

        return psl_data

    @classmethod
    def iter_text(cls, psl_text):
        no_header_psl = cls._remove_header(psl_text)
        for line in no_header_psl.splitlines():
            if not line:
                continue
            yield cls._parse_line(line)

    @classmethod
    def iter_file(cls, psl_file):
        with open(psl_file) as f_in:
            first_line = f_in.readline()
            if not first_line:
                return

            if first_line.startswith('psLayout'):
                for _ in range(4):
                    f_in.readline()
            else:
                first_line = first_line.rstrip('\n')
                if first_line:
                    yield cls._parse_line(first_line)

            for line in f_in:
                line = line.rstrip('\n')
                if not line:
                    continue
                yield cls._parse_line(line)

    @classmethod
    def _parse_line(cls, line):
        fields = line.split('\t')
        return cls._psl(
            int(fields[0]),
            int(fields[1]),
            int(fields[2]),
            int(fields[3]),
            int(fields[4]),
            int(fields[5]),
            int(fields[6]),
            int(fields[7]),
            fields[8],
            fields[9],
            int(fields[10]),
            int(fields[11]),
            int(fields[12]),
            fields[13],
            int(fields[14]),
            int(fields[15]),
            int(fields[16]),
            int(fields[17]),
            cls._parse_int_list(fields[18]),
            cls._parse_int_list(fields[19]),
            cls._parse_int_list(fields[20]),
        )

    @staticmethod
    def _parse_int_list(value):
        return tuple(map(int, value.rstrip(',').split(',')))

    @staticmethod
    def _remove_header(psl_text):
        if psl_text.startswith('psLayout'):
            no_header_psl = psl_text.split('\n', 5)[5]
        else:
            no_header_psl = psl_text

        return no_header_psl


class PslFilters:
    @staticmethod
    def _is_colinear(match, threshold=COLINEAR_THRESHOLD):
        colinear_score = (match.matches + match.repMatches) / match.qSize
        return colinear_score > threshold

    @classmethod
    def colinear_filter(cls, psl_data, threshold=COLINEAR_THRESHOLD):
        for match in psl_data:
            if cls._is_colinear(match, threshold=threshold):
                yield match

    @staticmethod
    def _is_chimera(psl_data, threshold):
        chimera = PslChimera()
        chimera.parse(psl_data)

        if not chimera.linear:
            if (threshold.check(chimera.five_end) is True) and \
               (threshold.check(chimera.three_end) is True):
                return True

        return False

    @classmethod
    def chimera_filter(cls, psl_data,
                       min_matches=CHIMERA_MIN_MATCHES,
                       min_diff=CHIMERA_MIN_DIFF):
        threshold = _Threshold(min_matches=min_matches, min_diff=min_diff)
        _is_chimera = partial(cls._is_chimera, threshold=threshold)

        for _, gp in groupby(psl_data, key=lambda match: match.qName):
            gp = list(gp)
            if _is_chimera(gp):
                yield gp


class _Threshold:
    def __init__(self, min_matches, min_diff):
        self.min_matches = min_matches
        self.min_diff = min_diff

    def check(self, matches):
        num_matches = len(matches)
        if num_matches == 0:
            best_match_score = 0
            second_best_score = 0
        elif num_matches == 1:
            best_match_score = self._get_score(matches[0])
            second_best_score = 0
        elif num_matches == 2:
            best_match_score = self._get_score(matches[0])
            second_best_score = self._get_score(matches[1])

        if best_match_score < self.min_matches:
            return False

        if best_match_score - second_best_score < self.min_diff:
            return False

        return True

    @staticmethod
    def _get_score(match):
        return match.matches if hasattr(match, 'matches') else match


class PslChimera:
    def __init__(self):
        self._linear = False
        self._five_end = deque(maxlen=2)
        self._three_end = deque(maxlen=2)

    def parse(self, psl_data):
        for match in psl_data:
            if self._is_colinear(match):
                self._linear = True

            if self._is_five_end(match):
                self._add_best_match(self._five_end, match)

            if self._is_three_end(match):
                self._add_best_match(self._three_end, match)

    @property
    def linear(self):
        return self._linear

    @property
    def five_end(self):
        return self._five_end

    @property
    def three_end(self):
        return self._three_end

    @staticmethod
    def _is_colinear(match):
        return (match.qStart < END_DISTANCE) and \
            (match.qSize - match.qEnd < END_DISTANCE)

    @staticmethod
    def _is_five_end(match):
        return match.qStart < END_DISTANCE

    @staticmethod
    def _is_three_end(match):
        return match.qSize - match.qEnd < END_DISTANCE

    @staticmethod
    def _is_better_match(match, other):
        return match.matches >= other.matches

    @classmethod
    def _add_best_match(cls, q, match):
        len_q = len(q)

        if len_q == 0:
            q.append(match)
        elif len_q == 1:
            if cls._is_better_match(match, q[0]):
                q.appendleft(match)
            else:
                q.append(match)
        else:
            if cls._is_better_match(match, q[0]):
                q.appendleft(match)
            elif cls._is_better_match(match, q[1]):
                q.appendleft(match)
                q.reverse()
            else:
                pass


class PslChimeraScores:
    def __init__(self):
        self._linear = False
        self._five_end = deque(maxlen=2)
        self._three_end = deque(maxlen=2)

    def parse(self, match):
        if PslChimera._is_colinear(match):
            self._linear = True

        if PslChimera._is_five_end(match):
            self._add_best_score(self._five_end, match.matches)

        if PslChimera._is_three_end(match):
            self._add_best_score(self._three_end, match.matches)

    @property
    def linear(self):
        return self._linear

    @property
    def five_end(self):
        return self._five_end

    @property
    def three_end(self):
        return self._three_end

    @staticmethod
    def _add_best_score(q, score):
        len_q = len(q)

        if len_q == 0:
            q.append(score)
        elif len_q == 1:
            if score >= q[0]:
                q.appendleft(score)
            else:
                q.append(score)
        else:
            if score >= q[0]:
                q.appendleft(score)
            elif score >= q[1]:
                q.appendleft(score)
                q.reverse()


class PslFileSummary:
    """Single-pass summary of all query alignments in one PSL file."""

    def __init__(self):
        self.query_scores = {}
        self.colinear_ids = set()

    @classmethod
    def from_file(cls, psl_file):
        summary = cls()
        for match in BlatPsl.iter_file(psl_file):
            scores = summary.query_scores.setdefault(
                match.qName,
                PslChimeraScores()
            )
            scores.parse(match)
            if PslFilters._is_colinear(match):
                summary.colinear_ids.add(match.qName)
        return summary

    @property
    def hit_ids(self):
        return set(self.query_scores)

    def get_chimera_ids(self,
                        min_matches=CHIMERA_MIN_MATCHES,
                        min_diff=CHIMERA_MIN_DIFF):
        threshold = _Threshold(min_matches=min_matches, min_diff=min_diff)
        return {
            qname
            for qname, scores in self.query_scores.items()
            if not scores.linear and
            threshold.check(scores.five_end) and
            threshold.check(scores.three_end)
        }


class PslUtils:
    @staticmethod
    def get_uniq_qname(psl_data):
        uniq_qnames = sorted({match.qName for match in psl_data})
        return uniq_qnames


class AmbAlnChecker:
    def __init__(self,
                 ref_genome,
                 ref_others,
                 work_dir='.',
                 num_proc=1,
                 aligner='blat',
                 blat_bin='blat',
                 pblat_bin='pblat',
                 mp_blat_bin='mp_blat.py'):

        self.ref_genome = ref_genome
        self.ref_others = ref_others
        self.work_dir = work_dir
        self.num_proc = num_proc
        self.aligner = aligner
        self.blat_bin = blat_bin
        self.pblat_bin = pblat_bin
        self.mp_blat_bin = mp_blat_bin

        if self.aligner == 'blat':
            runner_options = {
                'work_dir': self.work_dir,
                'num_proc': self.num_proc,
                'blat_bin': self.blat_bin,
                'mp_blat_bin': self.mp_blat_bin
            }
            runner_class = Blat
        elif self.aligner == 'pblat':
            runner_options = {
                'work_dir': self.work_dir,
                'num_proc': self.num_proc,
                'pblat_bin': self.pblat_bin
            }
            runner_class = Pblat
        else:
            raise ValueError(f'Unsupported aligner: {self.aligner}')

        self.blat = runner_class(**runner_options)
        self.web_blat = runner_class(
            blat_opts=WEB_BLAT_OPTIONS,
            **runner_options
        )

    @property
    def infos(self):
        infos = {
            'ref_genome': self.ref_genome,
            'ref_others': self.ref_others,
            'work_dir': self.work_dir,
            'num_proc': self.num_proc,
            'aligner': self.aligner,
            'blat_bin': self.blat_bin,
            'pblat_bin': self.pblat_bin,
            'mp_blat_bin': self.mp_blat_bin
        }
        return infos

    @staticmethod
    def _get_colinear_ids(psl_file):
        return sorted({
            match.qName
            for match in BlatPsl.iter_file(psl_file)
            if PslFilters._is_colinear(match)
        })

    @staticmethod
    def _get_chimera_ids(psl_file,
                         min_matches=CHIMERA_MIN_MATCHES,
                         min_diff=CHIMERA_MIN_DIFF):
        threshold = _Threshold(min_matches=min_matches, min_diff=min_diff)
        chimera_by_qname = {}

        for match in BlatPsl.iter_file(psl_file):
            chimera = chimera_by_qname.setdefault(
                match.qName,
                PslChimeraScores()
            )
            chimera.parse(match)

        return sorted(
            qname
            for qname, chimera in chimera_by_qname.items()
            if not chimera.linear and
            threshold.check(chimera.five_end) and
            threshold.check(chimera.three_end)
        )

    @staticmethod
    def _get_multiple_hit_ids(psl_file, excluded_ids):
        excluded_ids = set(excluded_ids)
        return sorted({
            match.qName
            for match in BlatPsl.iter_file(psl_file)
            if match.qName not in excluded_ids
        })

    def check(self, reads_file):
        result = AmbAlnResult(reads_file=reads_file)
        result.set_infos(**self.infos)
        psl_files = []

        try:
            # run blat
            psl_rG_1 = self.blat(self.ref_genome, reads_file)
            psl_files.append(psl_rG_1)
            psl_rG_2 = self.web_blat(self.ref_genome, reads_file)
            psl_files.append(psl_rG_2)
            psl_rO_1 = self.blat(self.ref_others, reads_file)
            psl_files.append(psl_rO_1)
            psl_rO_2 = self.web_blat(self.ref_others, reads_file)
            psl_files.append(psl_rO_2)

            # Each PSL is parsed exactly once. The sets below retain the
            # original classification and merge semantics.
            summary_rG_1 = PslFileSummary.from_file(psl_rG_1)
            summary_rG_2 = PslFileSummary.from_file(psl_rG_2)
            summary_rO_1 = PslFileSummary.from_file(psl_rO_1)
            summary_rO_2 = PslFileSummary.from_file(psl_rO_2)

            result._psl_rG_1_colinear_ids = sorted(
                summary_rG_1.colinear_ids
            )
            result._psl_rG_2_colinear_ids = sorted(
                summary_rG_2.colinear_ids
            )
            result._psl_rO_1_colinear_ids = sorted(
                summary_rO_1.colinear_ids
            )
            result._psl_rO_2_colinear_ids = sorted(
                summary_rO_2.colinear_ids
            )

            result.colinear_ids = sorted(set(chain(
                result._psl_rG_1_colinear_ids,
                result._psl_rG_2_colinear_ids,
                result._psl_rO_1_colinear_ids,
                result._psl_rO_2_colinear_ids
            )))

            # multiple hits part
            result._psl_rG_1_chimera_ids = sorted(
                summary_rG_1.get_chimera_ids()
            )
            result._psl_rG_2_chimera_ids = sorted(
                summary_rG_2.get_chimera_ids()
            )

            psl_rG_1_excluded_ids = set(chain(
                result._psl_rG_1_colinear_ids,
                result._psl_rG_1_chimera_ids
            ))
            psl_rG_2_excluded_ids = set(chain(
                result._psl_rG_2_colinear_ids,
                result._psl_rG_2_chimera_ids
            ))

            result._psl_rG_1_multiple_hits_ids = sorted(
                summary_rG_1.hit_ids.difference(psl_rG_1_excluded_ids)
            )
            result._psl_rG_2_multiple_hits_ids = sorted(
                summary_rG_2.hit_ids.difference(psl_rG_2_excluded_ids)
            )

            result.multiple_hits_ids = sorted(
                set(chain(
                    result._psl_rG_1_multiple_hits_ids,
                    result._psl_rG_2_multiple_hits_ids
                )).difference(
                    set(result.colinear_ids)
                )
            )

            result.update_result()
        finally:
            for psl_file in psl_files:
                if os.path.exists(psl_file):
                    os.remove(psl_file)

        return result


class AmbAlnResult:
    def __init__(self, reads_file):
        self.reads_file = reads_file
        self.infos = {}

        self._init_result()

    def _init_result(self):
        self.read_ids = self._get_read_ids(self.reads_file)
        self.result = {read_id: [0, 0, 0] for read_id in self.read_ids}

    @classmethod
    def _get_read_ids(cls, reads_file):
        read_ids = []
        with open(reads_file) as f_in:
            for line in f_in:
                if line.startswith('>'):
                    read_ids.append(line[1:].split(None, 1)[0])
        return read_ids

    def set_infos(self, **kwargs):
        self.infos.update(kwargs)

    def update_result(self):
        for id_ in self.colinear_ids:
            self.result[id_][0] = 1

        for id_ in self.multiple_hits_ids:
            self.result[id_][1] = 1

        for id_, checking_result in self.result.items():
            if (checking_result[0] == 1) or (checking_result[1] == 1):
                self.result[id_][2] = 1

    def save_result(self, out_file):
        with open(out_file, 'w') as out:
            csv_writer = csv.writer(out, delimiter='\t')

            csv_writer.writerow(RESULT_HEADER)

            for read_id, check_items in self.result.items():
                csv_writer.writerow([read_id] + check_items)


def read_fasta_records(fasta_file):
    """Read FASTA records and validate IDs used by PSL and the cache."""
    records = OrderedDict()
    current_id = None
    sequence_parts = []

    def add_record(record_id, parts):
        if record_id is None:
            return
        sequence = ''.join(parts)
        if not sequence:
            raise ValueError(f'FASTA record has an empty sequence: {record_id}')
        if record_id in records:
            if records[record_id].sequence != sequence:
                raise ValueError(
                    f'Duplicate FASTA ID has different sequences: {record_id}'
                )
            return
        records[record_id] = FastaRecord(record_id, sequence)

    with open(fasta_file) as f_in:
        for line_number, line in enumerate(f_in, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                add_record(current_id, sequence_parts)
                header = line[1:].strip()
                if not header:
                    raise ValueError(
                        f'FASTA header has no ID at line {line_number}'
                    )
                current_id = header.split(None, 1)[0]
                sequence_parts = []
            else:
                if current_id is None:
                    raise ValueError(
                        f'FASTA sequence appears before a header at '
                        f'line {line_number}'
                    )
                sequence_parts.append(line)

    add_record(current_id, sequence_parts)
    if not records:
        raise ValueError('Input FASTA contains no records')
    return list(records.values())


def write_fasta_records(records, fasta_file):
    with open(fasta_file, 'w') as out:
        for record in records:
            print(f'>{record.id}', file=out)
            print(record.sequence, file=out)


def create_parser():

    parser = argparse.ArgumentParser(

        description=textwrap.dedent("""
            A python script to check for reads if there exist colinear alignments
             or multiple alignments.
            """),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-rg', dest='ref_genome', required=True)
    parser.add_argument('-ro', dest='ref_others', required=True)
    parser.add_argument('in_file', help='input reads file.')
    parser.add_argument('out_file')
    parser.add_argument(
        '-p', '--num-proc', dest='num_proc', type=int, default=1
    )
    parser.add_argument(
        '--aligner', choices=('blat', 'pblat'), default='blat',
        help='Alignment backend (default: blat).'
    )
    parser.add_argument('--blat-bin', default='blat')
    parser.add_argument('--pblat-bin', default='pblat')
    parser.add_argument('--mp-blat-bin', default='mp_blat.py')
    parser.add_argument(
        '--tmp-path',
        default='.',
        help='Directory for temporary FASTA and PSL files.'
    )

    return parser


def cli():
    parser = create_parser()
    args = parser.parse_args()

    if args.num_proc < 1:
        parser.error('--num-proc must be at least 1')
    if not os.path.isdir(args.tmp_path):
        parser.error(f'--tmp-path is not a directory: {args.tmp_path}')

    checker = AmbAlnChecker(
        ref_genome=args.ref_genome,
        ref_others=args.ref_others,
        work_dir=args.tmp_path,
        num_proc=args.num_proc,
        aligner=args.aligner,
        blat_bin=args.blat_bin,
        pblat_bin=args.pblat_bin,
        mp_blat_bin=args.mp_blat_bin
    )

    result = checker.check(args.in_file)
    result.save_result(args.out_file)


if __name__ == "__main__":
    cli()
