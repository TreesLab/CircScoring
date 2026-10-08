import subprocess
from pathlib import Path


pairs = (
    (snakemake.input.five_prime_bed, snakemake.output.five_prime_fa),
    (snakemake.input.three_prime_bed, snakemake.output.three_prime_fa),
)
states = [Path(bed).stat().st_size > 0 for bed, _ in pairs]
if len(set(states)) != 1:
    raise ValueError("Donor and acceptor BED files have inconsistent content")
for _, output in pairs:
    Path(output).parent.mkdir(parents=True, exist_ok=True)
if not states[0]:
    for _, output in pairs:
        Path(output).write_text("", encoding="utf-8")
else:
    with open(snakemake.log[0], "wb") as log_handle:
        for bed, output in pairs:
            subprocess.run(
                ["bedtools", "getfasta", "-fi", snakemake.input.genome,
                 "-bed", bed, "-s", "-name", "-fo", output],
                stderr=log_handle,
                check=True,
            )
