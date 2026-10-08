MOUSE_LIFTED = result_path("00_master/base_construction/mouse/liftOver_to_hg38/mouse_circRNA_presence_master.v4.hg38_1based.tsv")
MOUSE_LIFTED_UNMAPPED = result_path("00_master/base_construction/mouse/liftOver_to_hg38/mouse_circRNA_presence_master.v4.hg38_1based.unmapped.tsv")
MOUSE_LIFTED_SUMMARY = result_path("00_master/base_construction/mouse/liftOver_to_hg38/mouse_circRNA_presence_master.v4.hg38_1based.summary.tsv")
CONSERVATION = result_path("00_master/base_construction/mouse_conservation/circRNA_presence_master.v4.mouse_conservation.tsv")
CONSERVATION_SUMMARY = result_path("00_master/base_construction/mouse_conservation/circRNA_presence_master.v4.mouse_conservation.summary.tsv")


rule liftover_mouse_v4_to_hg38:
    input:
        mouse=str(base_dir("mouse") / "circRNA_presence_master.v4.tsv"),
        overrides=repo_path(config["mouse_partial_mapping_overrides"])
    output:
        mapped=MOUSE_LIFTED,
        unmapped=MOUSE_LIFTED_UNMAPPED,
        summary=MOUSE_LIFTED_SUMMARY
    params:
        chain=repo_path(config["chains"]["mm10_to_hg38"])
    shell:
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/liftover_mouse_to_hg38_1based.py "
        "--input {input.mouse:q} --output {output.mapped:q} "
        "--unmapped-output {output.unmapped:q} --summary-output {output.summary:q} "
        "--chain {params.chain:q} --liftover-bin {LIFTOVER:q} "
        "--partial-mapping-overrides {input.overrides:q} --liftover-rows-without-history"


rule build_mouse_conservation:
    input:
        human=str(base_dir("human") / "circRNA_presence_master.v4.tsv"),
        mouse=MOUSE_LIFTED
    output:
        table=CONSERVATION,
        summary=CONSERVATION_SUMMARY
    shell:
        "{PYTHON:q} {REPO}/00_master/pipeline_scripts/build_mouse_conservation.py "
        "--human {input.human:q} --mouse-lifted {input.mouse:q} "
        "--output {output.table:q} --summary-output {output.summary:q}"


rule mouse_conservation:
    input:
        CONSERVATION,
        CONSERVATION_SUMMARY
