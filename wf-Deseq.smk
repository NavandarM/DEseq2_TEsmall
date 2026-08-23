import yaml
import os
import csv
from itertools import combinations

configfile: "deseq2_config.yaml"

Indir = config["Indir"]
Outdir = config["Outdir"]
MetaData = config["MetaData"]
smallestGroupSize = config['smallestGroupSize']
LType = config['Type']
RefCondition = config.get('RefCondition', "")

DE_dir = os.path.join(Outdir, 'DEoutput')

vsc = os.path.join(DE_dir, 'VarianceStabilizedCounts.txt')
nrc = os.path.join(DE_dir, 'Normalized_Read_Counts.txt')
pca = os.path.join(DE_dir, 'Pcaplot.pdf')
flag = os.path.join(DE_dir, 'deseq_results.done')


def get_conditions(metadata_path, ref_condition=""):
    """Unique conditions in the same level order DESeq2's factor()/relevel()
    would use: alphabetical, with ref_condition (if any) moved to the front."""
    with open(metadata_path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        conditions = sorted({row['condition'] for row in reader})
    if ref_condition:
        if ref_condition not in conditions:
            raise WorkflowError(
                f"RefCondition '{ref_condition}' not found in metadata conditions: {conditions}"
            )
        conditions = [ref_condition] + [c for c in conditions if c != ref_condition]
    return conditions


Conditions = get_conditions(MetaData, RefCondition)

if len(Conditions) < 2:
    raise WorkflowError(f"Need at least 2 conditions in metadata, found: {Conditions}")

IsMulti = len(Conditions) > 2

if IsMulti:
    Pairs = list(combinations(Conditions, 2))
    pair_outputs = []
    for ref_level, comp_level in Pairs:
        tag = f"condition_{comp_level}_vs_{ref_level}"
        pair_outputs += [
            os.path.join(DE_dir, f"volcanoplot_{tag}.pdf"),
            os.path.join(DE_dir, f"DEG_{tag}.txt"),
        ]

    rule all:
        input:
            vsc, nrc, pca, flag, *pair_outputs

    rule deseq2_analysis_multi:
        input:
            counts=Indir,
            metadata=MetaData
        output:
            vsc=vsc,
            nrc=nrc,
            pca=pca,
            pairs=pair_outputs,
            flag=touch(flag)
        params:
            DE=DE_dir,
            groupSz=smallestGroupSize,
            Type=LType,
            ref=RefCondition
        conda:
            config["DeseqEnv"]
        shell:
            """
            mkdir -p Robj
            Rscript scripts/runDeseq_multi.R "{input.counts}" "{input.metadata}" "{output.vsc}" "{output.nrc}" "{output.pca}" "{params.DE}" "{params.groupSz}" "{params.Type}" "{params.ref}"
            """
else:
    vol = os.path.join(DE_dir, 'volcanoplot.pdf')

    rule all:
        input:
            vsc, nrc, pca, vol, flag

    rule deseq2_analysis:
        input:
            counts=Indir,
            metadata=MetaData
        output:
            vsc=vsc,
            nrc=nrc,
            pca=pca,
            vol=vol,
            flag=touch(flag)
        params:
            DE=DE_dir,
            groupSz=smallestGroupSize,
            Type=LType,
            ref=RefCondition
        conda:
            config["DeseqEnv"]
        shell:
            """
            mkdir -p Robj
            Rscript scripts/runDeseq.R "{input.counts}" "{input.metadata}" "{output.vsc}" "{output.nrc}" "{output.pca}" "{output.vol}" "{params.DE}" "{params.groupSz}" "{params.Type}" "{params.ref}"
            """
