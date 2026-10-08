MOUSE_T1_DATABASES = sorted(MOUSE_T1)
MOUSE_T1_REGEX = "(?:" + "|".join(MOUSE_T1_DATABASES) + ")"
MOUSE_T2_ROWS = read_tsv(MOUSE_T2_MANIFEST)
MOUSE_T2_FINALS = [result_path(row["final_output"]) for row in MOUSE_T2_ROWS]
MOUSE_T3_ROWS = read_tsv(MOUSE_T3_MANIFEST)
MOUSE_T3_FINALS = [result_path(row["final_output"]) for row in MOUSE_T3_ROWS]
MOUSE_T4_ROWS = read_tsv(MOUSE_T4_MANIFEST)
MOUSE_T4_FINALS = [result_path(row["final_output"]) for row in MOUSE_T4_ROWS]


rule mouse_tier1_final:
    input:
        clean=lambda wc: repo_path(MOUSE_T1[wc.database]["clean_input"])
    output:
        converted=result_path("{database}/mouse/processed/{database}.circRNAs.mm10_1based.tsv"),
        unmapped=result_path("{database}/mouse/processed/{database}.circRNAs.mm10_1based.unmapped.tsv"),
        summary=result_path("{database}/mouse/processed/{database}.circRNAs.mm10_1based.summary.tsv"),
        final=result_path("{database}/mouse/processed/{database}.circRNAs.final.tsv")
    params:
        coordinate_system=lambda wc: MOUSE_T1[wc.database]["coordinate_system"],
        chain=repo_path(config["chains"]["mm9_to_mm10"])
    wildcard_constraints:
        database=MOUSE_T1_REGEX
    shell:
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/liftover_to_1based.py "
        "--input {input.clean:q} --source-coordinate-system {params.coordinate_system:q} "
        "--target-assembly mm10 --output {output.converted:q} "
        "--unmapped-output {output.unmapped:q} --summary-output {output.summary:q} "
        "--chain {params.chain:q} --liftover-bin {LIFTOVER:q} && "
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/write_final.py "
        "--input {output.converted:q} --output {output.final:q}"


rule mouse_tier2_finals:
    input:
        clean=[repo_path(row["clean_input"]) for row in MOUSE_T2_ROWS],
        tier1=base_outputs("mouse", 1)
    output:
        MOUSE_T2_FINALS
    threads: RUNNER_JOBS
    params:
        chain=repo_path(config["chains"]["mm9_to_mm10"])
    shell:
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/build_mouse_tier2_finals.py "
        "--manifest {MOUSE_T2_MANIFEST:q} --input-root {REPO:q} --output-root {RESULTS:q} "
        "--chain {params.chain:q} --liftover-bin {LIFTOVER:q} --jobs {threads}"


rule mouse_tier3_finals:
    input:
        clean=[repo_path(row["clean_input"]) for row in MOUSE_T3_ROWS],
        locked=[
            repo_path(row["coordinate_check"])
            for row in MOUSE_T3_ROWS
            if row["coordinate_strategy"] == "per_row"
        ],
        presence=str(base_dir("mouse") / "circRNA_presence_master.v2.tsv"),
        ids=str(base_dir("mouse") / "circRNA_presence_master_circ_ids.v2.tsv")
    output:
        MOUSE_T3_FINALS
    threads: RUNNER_JOBS
    params:
        chain=repo_path(config["chains"]["mm9_to_mm10"])
    shell:
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/build_mouse_tier3_finals.py "
        "--manifest {MOUSE_T3_MANIFEST:q} --input-root {REPO:q} --output-root {RESULTS:q} "
        "--presence {input.presence:q} --ids {input.ids:q} --chain {params.chain:q} "
        "--liftover-bin {LIFTOVER:q} --jobs {threads}"


rule mouse_tier4_finals:
    input:
        clean=[repo_path(row["clean_input"]) for row in MOUSE_T4_ROWS],
        locked=[
            repo_path(row["coordinate_check"])
            for row in MOUSE_T4_ROWS
            if row["coordinate_strategy"] == "precomputed_row_level"
        ],
        presence=str(base_dir("mouse") / "circRNA_presence_master.v3.tsv"),
        ids=str(base_dir("mouse") / "circRNA_presence_master_circ_ids.v3.tsv")
    output:
        MOUSE_T4_FINALS
    threads: RUNNER_JOBS
    params:
        chain=repo_path(config["chains"]["mm9_to_mm10"])
    shell:
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/build_mouse_tier4_finals.py "
        "--manifest {MOUSE_T4_MANIFEST:q} --input-root {REPO:q} --output-root {RESULTS:q} "
        "--presence {input.presence:q} --ids {input.ids:q} --chain {params.chain:q} "
        "--liftover-bin {LIFTOVER:q} --jobs {threads}"


rule mouse_v1:
    input: base_outputs("mouse", 1)

rule mouse_v2:
    input: base_outputs("mouse", 2)

rule mouse_v3:
    input: base_outputs("mouse", 3)

rule mouse_v4:
    input: base_outputs("mouse", 4)
