HUMAN_T1_DATABASES = sorted(HUMAN_T1)
HUMAN_T1_REGEX = "(?:" + "|".join(HUMAN_T1_DATABASES) + ")"
HUMAN_T2_ROWS = read_tsv(HUMAN_T2_MANIFEST)
HUMAN_T2_NORMAL = [row for row in HUMAN_T2_ROWS if row["database"] not in {"21_CircTarget", "43_circRNADisease"}]
HUMAN_T2_NORMAL_FINALS = [result_path(row["final_output"]) for row in HUMAN_T2_NORMAL]
HUMAN_T2_NORMAL_ARGS = " ".join(f"--database {row['database']}" for row in HUMAN_T2_NORMAL)
HUMAN_T3_ROWS = read_tsv(HUMAN_T3_MANIFEST)
HUMAN_T3_DATABASES = sorted({row["database"] for row in HUMAN_T3_ROWS})
HUMAN_T3_FINALS = [result_path(f"{database}/processed/{database}.circRNAs.final.tsv") for database in HUMAN_T3_DATABASES]
HUMAN_T4_ROWS = read_tsv(HUMAN_T4_MANIFEST)
HUMAN_T4_DATABASES = sorted({row["database"] for row in HUMAN_T4_ROWS})
HUMAN_T4_FINALS = [result_path(f"{database}/processed/{database}.circRNAs.final.tsv") for database in HUMAN_T4_DATABASES]


rule human_tier1_final:
    input:
        clean=lambda wc: repo_path(HUMAN_T1[wc.database]["clean_input"])
    output:
        converted=result_path("{database}/processed/{database}.circRNAs.hg38_1based.tsv"),
        unmapped=result_path("{database}/processed/{database}.circRNAs.hg38_1based.unmapped.tsv"),
        summary=result_path("{database}/processed/{database}.circRNAs.hg38_1based.summary.tsv"),
        final=result_path("{database}/processed/{database}.circRNAs.final.tsv")
    params:
        coordinate_system=lambda wc: HUMAN_T1[wc.database]["coordinate_system"],
        chain=repo_path(config["chains"]["hg19_to_hg38"])
    wildcard_constraints:
        database=HUMAN_T1_REGEX
    shell:
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/liftover_to_hg38_1based.py "
        "--input {input.clean:q} --coordinate-system {params.coordinate_system:q} "
        "--output {output.converted:q} --unmapped-output {output.unmapped:q} "
        "--summary-output {output.summary:q} --chain {params.chain:q} --liftover-bin {LIFTOVER:q} && "
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/write_final.py "
        "--input {output.converted:q} --output {output.final:q}"


rule resolve_circtarget_ids:
    input:
        ids=repo_path("21_CircTarget/clean/21_CircTarget.circRNA_ids.clean.tsv"),
        circbase=result_path("01_circBase/processed/01_circBase.circRNAs.final.tsv"),
        circatlas=result_path("02_circAtlas/processed/02_circAtlas.circRNAs.final.tsv")
    output:
        clean=result_path("21_CircTarget/clean/21_CircTarget.circRNAs.clean.tsv"),
        audit=result_path("21_CircTarget/clean/21_CircTarget.circRNAs.clean.audit.tsv"),
        summary=result_path("21_CircTarget/clean/21_CircTarget.circRNAs.clean.summary.tsv")
    shell:
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/resolve_circrna_ids.py "
        "--input {input.ids:q} --reference circBase={input.circbase:q} "
        "--reference circAtlas={input.circatlas:q} --output {output.clean:q} "
        "--audit {output.audit:q} --summary-output {output.summary:q} --force"


rule build_circtarget_final:
    input:
        clean=rules.resolve_circtarget_ids.output.clean
    output:
        converted=result_path("21_CircTarget/processed/21_CircTarget.circRNAs.hg38_1based.tsv"),
        unmapped=result_path("21_CircTarget/processed/21_CircTarget.circRNAs.hg38_1based.unmapped.tsv"),
        summary=result_path("21_CircTarget/processed/21_CircTarget.circRNAs.hg38_1based.summary.tsv"),
        final=result_path("21_CircTarget/processed/21_CircTarget.circRNAs.final.tsv")
    shell:
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/liftover_to_hg38_1based.py "
        "--input {input.clean:q} --coordinate-system hg38_1based "
        "--output {output.converted:q} --unmapped-output {output.unmapped:q} "
        "--summary-output {output.summary:q} && "
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/write_final.py "
        "--input {output.converted:q} --output {output.final:q}"


rule resolve_circrnadisease_ids:
    input:
        ids=repo_path("43_circRNADisease/clean/43_circRNADisease.circRNA_ids.clean.tsv"),
        circbase=result_path("01_circBase/processed/01_circBase.circRNAs.final.tsv")
    output:
        clean=result_path("43_circRNADisease/clean/43_circRNADisease.circRNAs.clean.tsv"),
        audit=result_path("43_circRNADisease/clean/43_circRNADisease.circRNAs.clean.audit.tsv"),
        summary=result_path("43_circRNADisease/clean/43_circRNADisease.circRNAs.clean.summary.tsv")
    shell:
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/resolve_circrna_ids.py "
        "--input {input.ids:q} --reference circBase={input.circbase:q} "
        "--output {output.clean:q} --audit {output.audit:q} "
        "--summary-output {output.summary:q} --force"


rule build_circrnadisease_final:
    input:
        clean=rules.resolve_circrnadisease_ids.output.clean
    output:
        final=result_path("43_circRNADisease/processed/43_circRNADisease.circRNAs.final.tsv")
    shell:
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/write_final.py "
        "--input {input.clean:q} --output {output.final:q}"


rule human_tier2_normal_finals:
    input:
        clean=[repo_path(row["clean_input"]) for row in HUMAN_T2_NORMAL],
        tier1=base_outputs("human", 1)
    output:
        HUMAN_T2_NORMAL_FINALS
    threads: RUNNER_JOBS
    params:
        chain=repo_path(config["chains"]["hg19_to_hg38"]),
        databases=HUMAN_T2_NORMAL_ARGS
    shell:
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/build_human_tier2_finals.py "
        "--manifest {HUMAN_T2_MANIFEST:q} --input-root {REPO:q} --output-root {RESULTS:q} "
        "--chain {params.chain:q} --liftover-bin {LIFTOVER:q} --jobs {threads} "
        "{params.databases}"


rule human_tier3_finals:
    input:
        clean=[repo_path(row["clean_input"]) for row in HUMAN_T3_ROWS],
        presence=str(base_dir("human") / "circRNA_presence_master.v2.tsv"),
        ids=str(base_dir("human") / "circRNA_presence_master_circ_ids.v2.tsv")
    output:
        HUMAN_T3_FINALS
    threads: RUNNER_JOBS
    params:
        chain=repo_path(config["chains"]["hg19_to_hg38"])
    shell:
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/build_human_tier3_finals.py "
        "--manifest {HUMAN_T3_MANIFEST:q} --input-root {REPO:q} --output-root {RESULTS:q} "
        "--presence {input.presence:q} --ids {input.ids:q} --chain {params.chain:q} "
        "--liftover-bin {LIFTOVER:q} --jobs {threads}"


rule human_tier4_finals:
    input:
        clean=[repo_path(row["clean_input"]) for row in HUMAN_T4_ROWS],
        presence=str(base_dir("human") / "circRNA_presence_master.v3.tsv"),
        ids=str(base_dir("human") / "circRNA_presence_master_circ_ids.v3.tsv")
    output:
        HUMAN_T4_FINALS
    threads: RUNNER_JOBS
    params:
        chain=repo_path(config["chains"]["hg19_to_hg38"])
    shell:
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/build_human_tier4_finals.py "
        "--manifest {HUMAN_T4_MANIFEST:q} --input-root {REPO:q} --output-root {RESULTS:q} "
        "--presence {input.presence:q} --ids {input.ids:q} --chain {params.chain:q} "
        "--liftover-bin {LIFTOVER:q} --jobs {threads}"


rule human_v1:
    input: base_outputs("human", 1)

rule human_v2:
    input: base_outputs("human", 2)

rule human_v3:
    input: base_outputs("human", 3)

rule human_v4:
    input: base_outputs("human", 4)
