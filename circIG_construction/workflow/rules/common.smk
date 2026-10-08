BASE_NAMES = {
    "master": "master_circRNAs.v{tier}.txt",
    "presence": "circRNA_presence_master.v{tier}.tsv",
    "ids": "circRNA_presence_master_circ_ids.v{tier}.tsv",
    "compact": "circRNA_presence_master_compact.v{tier}.tsv",
    "compact_ids": "circRNA_presence_master_compact_circ_ids.v{tier}.tsv",
    "sources": "database_sources.v{tier}.tsv",
}


def base_dir(species):
    suffix = "00_master/base_construction" if species == "human" else "00_master/base_construction/mouse"
    return RESULTS / suffix


def base_outputs(species, tier):
    root = base_dir(species)
    return [str(root / BASE_NAMES[key].format(tier=tier)) for key in ("master", "presence", "ids", "compact", "compact_ids")]


def expected_source_outputs(species, tier):
    manifests = HUMAN_SOURCE_MANIFESTS if species == "human" else MOUSE_SOURCE_MANIFESTS
    return source_outputs(manifests[tier])


rule build_human_base:
    input:
        manifest=lambda wc: str(HUMAN_SOURCE_MANIFESTS[int(wc.tier)]),
        finals=lambda wc: expected_source_outputs("human", int(wc.tier))
    output:
        master=str(base_dir("human") / "master_circRNAs.v{tier}.txt"),
        presence=str(base_dir("human") / "circRNA_presence_master.v{tier}.tsv"),
        ids=str(base_dir("human") / "circRNA_presence_master_circ_ids.v{tier}.tsv"),
        compact=str(base_dir("human") / "circRNA_presence_master_compact.v{tier}.tsv"),
        compact_ids=str(base_dir("human") / "circRNA_presence_master_compact_circ_ids.v{tier}.tsv"),
        sources=str(base_dir("human") / "database_sources.v{tier}.tsv")
    wildcard_constraints:
        tier="[1-4]"
    params:
        outdir=str(base_dir("human"))
    shell:
        "mkdir -p {params.outdir:q} && cp {input.manifest:q} {output.sources:q} && "
        "{PYTHON:q} {REPO}/00_master/base_construction/build_circRNA_presence.py "
        "--config {output.sources:q} --master-out {output.master:q} "
        "--presence-out {output.presence:q} --ids-out {output.ids:q} "
        "--compact-out {output.compact:q} --compact-ids-out {output.compact_ids:q}"


rule build_mouse_base:
    input:
        manifest=lambda wc: str(MOUSE_SOURCE_MANIFESTS[int(wc.tier)]),
        finals=lambda wc: expected_source_outputs("mouse", int(wc.tier))
    output:
        master=str(base_dir("mouse") / "master_circRNAs.v{tier}.txt"),
        presence=str(base_dir("mouse") / "circRNA_presence_master.v{tier}.tsv"),
        ids=str(base_dir("mouse") / "circRNA_presence_master_circ_ids.v{tier}.tsv"),
        compact=str(base_dir("mouse") / "circRNA_presence_master_compact.v{tier}.tsv"),
        compact_ids=str(base_dir("mouse") / "circRNA_presence_master_compact_circ_ids.v{tier}.tsv"),
        sources=str(base_dir("mouse") / "database_sources.v{tier}.tsv")
    wildcard_constraints:
        tier="[1-4]"
    params:
        outdir=str(base_dir("mouse"))
    shell:
        "mkdir -p {params.outdir:q} && cp {input.manifest:q} {output.sources:q} && "
        "{PYTHON:q} {REPO}/00_master/base_construction/build_circRNA_presence.py "
        "--config {output.sources:q} --master-out {output.master:q} "
        "--presence-out {output.presence:q} --ids-out {output.ids:q} "
        "--compact-out {output.compact:q} --compact-ids-out {output.compact_ids:q}"
