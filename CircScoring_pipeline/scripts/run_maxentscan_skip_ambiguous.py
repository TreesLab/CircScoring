import argparse
import csv
import os
import subprocess
import sys
import tempfile


VALID_BASES = set('ACGT')


def parse_args():
    parser = argparse.ArgumentParser(
        description='Run MaxEntScan while assigning NA to FASTA records with ambiguous bases.'
    )
    parser.add_argument('--mode', choices=('5', '3'), required=True)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    return parser.parse_args()


def fasta_records(path):
    header = None
    sequence_parts = []

    with open(path, encoding='utf-8') as f_in:
        for line in f_in:
            line = line.rstrip('\n')
            if not line:
                continue

            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(sequence_parts)
                header = line
                sequence_parts = []
            else:
                sequence_parts.append(line.strip())

    if header is not None:
        yield header, ''.join(sequence_parts)


def is_valid_sequence(sequence):
    return bool(sequence) and set(sequence.upper()) <= VALID_BASES


def write_wrapped_fasta_record(handle, header, sequence):
    handle.write(f'{header}\n')
    handle.write(f'{sequence.upper()}\n')


def write_valid_inputs(input_fasta, valid_fasta, status_tsv):
    total_count = 0
    valid_count = 0
    invalid_count = 0

    with open(valid_fasta, 'w', encoding='utf-8') as valid_out, \
            open(status_tsv, 'w', newline='', encoding='utf-8') as status_out:
        writer = csv.writer(status_out, delimiter='\t', lineterminator='\n')

        for header, sequence in fasta_records(input_fasta):
            total_count += 1
            sequence = sequence.upper()

            if is_valid_sequence(sequence):
                valid_count += 1
                writer.writerow(['valid'])
                write_wrapped_fasta_record(valid_out, header, sequence)
            else:
                invalid_count += 1
                writer.writerow(['invalid', sequence])

    return total_count, valid_count, invalid_count


def run_maxentscan(mode, input_fasta, output_path):
    command = f'maxentscan_score{mode}.pl'
    with open(output_path, 'w', encoding='utf-8') as out:
        subprocess.run([command, input_fasta], stdout=out, check=True)


def read_next_result(results_handle):
    for line in results_handle:
        line = line.rstrip('\n')
        if line:
            return line
    return None


def rebuild_full_output(status_tsv, valid_results, output_path):
    valid_results_used = 0

    with open(status_tsv, newline='', encoding='utf-8') as status_in, \
            open(valid_results, encoding='utf-8') as results_in, \
            open(output_path, 'w', encoding='utf-8') as out:
        reader = csv.reader(status_in, delimiter='\t')

        for row_num, row in enumerate(reader, start=1):
            if not row:
                raise ValueError(f'Invalid empty status row {row_num}')

            status = row[0]
            if status == 'valid':
                if len(row) != 1:
                    raise ValueError(f'Invalid valid status row {row_num}: {row!r}')
                result_line = read_next_result(results_in)
                if result_line is None:
                    raise ValueError(
                        f'MaxEntScan returned fewer rows than expected; missing valid result at status row {row_num}.'
                    )
                valid_results_used += 1
                out.write(f'{result_line}\n')
            elif status == 'invalid':
                if len(row) != 2:
                    raise ValueError(f'Invalid invalid status row {row_num}: {row!r}')
                sequence = row[1]
                out.write(f'{sequence}\tNA\n')
            else:
                raise ValueError(f'Invalid status value at row {row_num}: {status!r}')

        extra_result = read_next_result(results_in)
        if extra_result is not None:
            raise ValueError(
                'MaxEntScan returned more rows than expected; first extra result line: '
                f'{extra_result!r}'
            )

    return valid_results_used


def main():
    args = parse_args()
    output_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(output_dir, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix='.maxentscan_tmp.', dir=output_dir) as temp_dir:
        valid_fasta = os.path.join(temp_dir, 'valid.fa')
        status_tsv = os.path.join(temp_dir, 'status.tsv')
        valid_results = os.path.join(temp_dir, 'valid.MaxEntScan_results')

        total_count, valid_count, invalid_count = write_valid_inputs(
            args.input,
            valid_fasta,
            status_tsv
        )

        if valid_count > 0:
            run_maxentscan(args.mode, valid_fasta, valid_results)
        else:
            open(valid_results, 'w', encoding='utf-8').close()

        valid_results_used = rebuild_full_output(status_tsv, valid_results, args.output)

        if valid_results_used != valid_count:
            raise ValueError(
                f'Expected {valid_count} valid MaxEntScan results, used {valid_results_used}.'
            )

    print(
        f'Processed {total_count} FASTA records: '
        f'{valid_count} valid, {invalid_count} ambiguous records assigned NA.',
        file=sys.stderr
    )


if __name__ == '__main__':
    main()
