# Evaluation of assembly-based SV calling methods

This snakemake pipeline was used to produce the results in our manuscript [SVIM-asm: Structural variant detection from haploid and diploid genome assemblies](https://doi.org/10.1101/2020.10.27.356907).

To run the pipeline:

```
#Clone repository
git clone https://github.com/eldariont/assembly_sv_calling_evaluation.git

#Change to dir
cd assembly_sv_calling_evaluation

#Add paths to input data
mv config/config_example.yaml config/config.yaml
vim config/config.yaml

#Run pipeline with x jobs running at the same time
snakemake -j x
```
