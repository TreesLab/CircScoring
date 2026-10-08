import csv
import re
import logging


logging.basicConfig(filename=snakemake.log[0], encoding='utf-8', level=logging.INFO)


def normalize_fasta_header(header):
    event_id = header.removeprefix('>').split()[0]
    event_id = event_id.split('::', 1)[0]
    event_id = re.sub(r'\([+-]\)$', '', event_id)
    return event_id


def read_fasta_event_ids(path):
    with open(path, encoding='utf-8') as f_in:
        for line in f_in:
            if line.startswith('>'):
                yield normalize_fasta_header(line.rstrip('\n'))


def read_maxentscan_results(path):
    with open(path, newline='', encoding='utf-8') as f_in:
        reader = csv.reader(f_in, delimiter='\t')
        for line_num, row in enumerate(reader, start=1):
            if not row:
                continue

            score = row[-1]
            yield line_num, score


donor_event_ids = list(read_fasta_event_ids(snakemake.input.five_prime_fa))
acceptor_event_ids = list(read_fasta_event_ids(snakemake.input.three_prime_fa))
donor_scores = list(read_maxentscan_results(snakemake.input.five_prime_results))
acceptor_scores = list(read_maxentscan_results(snakemake.input.three_prime_results))

expected_counts = {
    'donor FASTA records': len(donor_event_ids),
    'acceptor FASTA records': len(acceptor_event_ids),
    'donor MaxEntScan results': len(donor_scores),
    'acceptor MaxEntScan results': len(acceptor_scores),
}

if len(set(expected_counts.values())) != 1:
    raise ValueError("MaxEntScan input/result count mismatch: " + repr(expected_counts))

with open(snakemake.output[0], 'w', newline='', encoding='utf-8') as out:
    writer = csv.writer(out, delimiter='\t', lineterminator='\n')
    writer.writerow(['event_id', 'MAXENT(donor)', 'MAXENT(acceptor)'])

    for record_num, (donor_event_id, acceptor_event_id, donor_result, acceptor_result) in enumerate(
        zip(donor_event_ids, acceptor_event_ids, donor_scores, acceptor_scores),
        start=1
    ):
        if donor_event_id != acceptor_event_id:
            raise ValueError(
                f"MaxEntScan FASTA event ID mismatch at record {record_num}: "
                f"donor has {donor_event_id!r}, acceptor has {acceptor_event_id!r}"
            )

        _, donor_score = donor_result
        _, acceptor_score = acceptor_result
        writer.writerow([donor_event_id, donor_score, acceptor_score])
