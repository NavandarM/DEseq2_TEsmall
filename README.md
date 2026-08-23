# DEseq2_workflow
Differential analysis: 
- Differential gene expression analysis: between two conditions
- Metagenomics: differential microbial community
- small RNA (output specificially from <a href="https://github.com/mhammell-laboratory/TEsmall">TEsmall</a> or any tubular raw expression count data)
<br><br>
Usages: snakemake -s wf-Deseq.smk --use-conda

snakemake pipeline for the DE small RNAs

Requirements:
Edit the config file and update below variables:
-   Indir
-   Outdir
-   MetaData
-   smallestGroupSize (Smallest number of samples per group)
-   Type: "" # specify if it is smallRNAs otherwise keept --> ""
-   RefCondition: "" # optional reference/control condition name; defaults to the
    alphabetically-first condition if left empty

The workflow automatically detects how many conditions are present in
`MetaData`'s `condition` column:
-   **2 conditions** -> a single comparison (RefCondition vs the other condition),
    same outputs as before (`volcanoplot.pdf`, `DEG_<contrast>.txt`, etc.)
-   **3+ conditions** -> every pairwise contrast between conditions is run
    automatically (using one DESeq2 fit on the full dataset), producing
    `volcanoplot_condition_<B>_vs_<A>.pdf` and `DEG_condition_<B>_vs_<A>.txt`
    per pair, alongside the shared PCA plot, normalized counts, and VST counts.
