QC_ROOT = RESULTS / "qc" / "coordinate_system"


rule coordinate_qc:
    input:
        human_v1=str(base_dir("human") / "circRNA_presence_master.v1.tsv"),
        human_v2=str(base_dir("human") / "circRNA_presence_master.v2.tsv"),
        human_v3=str(base_dir("human") / "circRNA_presence_master.v3.tsv"),
        mouse_v1=str(base_dir("mouse") / "circRNA_presence_master.v1.tsv"),
        mouse_v2=str(base_dir("mouse") / "circRNA_presence_master.v2.tsv"),
        mouse_v3=str(base_dir("mouse") / "circRNA_presence_master.v3.tsv")
    output:
        str(QC_ROOT / ".complete")
    params:
        human_chain=repo_path(config["chains"]["hg19_to_hg38"]),
        mouse_chain=repo_path(config["chains"]["mm9_to_mm10"])
    shell:
        "{PYTHON:q} {REPO}/workflow/scripts/run_coordinate_qc.py --species human --tier 2 "
        "--manifest {HUMAN_T2_MANIFEST:q} --input-root {REPO:q} --base-dataset {input.human_v1:q} "
        "--output-dir {QC_ROOT}/human/tier2 --chain {params.human_chain:q} --liftover-bin {LIFTOVER:q} && "
        "{PYTHON:q} {REPO}/workflow/scripts/run_coordinate_qc.py --species human --tier 3 "
        "--manifest {HUMAN_T3_MANIFEST:q} --input-root {REPO:q} --base-dataset {input.human_v2:q} "
        "--output-dir {QC_ROOT}/human/tier3 --chain {params.human_chain:q} --liftover-bin {LIFTOVER:q} && "
        "{PYTHON:q} {REPO}/workflow/scripts/run_coordinate_qc.py --species human --tier 4 "
        "--manifest {HUMAN_T4_MANIFEST:q} --input-root {REPO:q} --base-dataset {input.human_v3:q} "
        "--output-dir {QC_ROOT}/human/tier4 --chain {params.human_chain:q} --liftover-bin {LIFTOVER:q} && "
        "{PYTHON:q} {REPO}/workflow/scripts/run_coordinate_qc.py --species mouse --tier 2 "
        "--manifest {MOUSE_T2_MANIFEST:q} --input-root {REPO:q} --base-dataset {input.mouse_v1:q} "
        "--output-dir {QC_ROOT}/mouse/tier2 --chain {params.mouse_chain:q} --liftover-bin {LIFTOVER:q} && "
        "{PYTHON:q} {REPO}/workflow/scripts/run_coordinate_qc.py --species mouse --tier 3 "
        "--manifest {MOUSE_T3_MANIFEST:q} --input-root {REPO:q} --base-dataset {input.mouse_v2:q} "
        "--output-dir {QC_ROOT}/mouse/tier3 --chain {params.mouse_chain:q} --liftover-bin {LIFTOVER:q} && "
        "{PYTHON:q} {REPO}/workflow/scripts/run_coordinate_qc.py --species mouse --tier 4 "
        "--manifest {MOUSE_T4_MANIFEST:q} --input-root {REPO:q} --base-dataset {input.mouse_v3:q} "
        "--output-dir {QC_ROOT}/mouse/tier4 --chain {params.mouse_chain:q} --liftover-bin {LIFTOVER:q} && "
        "mkdir -p {QC_ROOT:q} && printf 'completed\\n' > {output:q}"


rule workflow_check:
    shell:
        "{PYTHON:q} {REPO}/workflow/scripts/validate_workflow.py "
        "--repo {REPO:q} --config {REPO}/workflow/config/config.yaml"
