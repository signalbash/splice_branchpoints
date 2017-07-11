## Genome-wide predictions of branchpoints

**1. Install branchpointer
Open R
# version 0.99.13
> devtoools::install_github("betsig/branchpointer")

**2. Download Genome Annotations
```
cd data/

wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/GRCh38.primary_assembly.genome.fa.gz

gunzip gencode.v26.annotation.gtf.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz

gunzip gencode.v19.annotation.gtf.gz
gunzip GRCh37.p13.genome.fa.gz

cd ../

```
**3. Run Branchpointer
split into 10000 size chunks
279034 (gencode.v26) sites
274013 (gencode.v19) sites

```
for i in `seq 1 27`;
do Rscript scripts/genome_predictions/predictBranchpointsGenomeAnnotation.R --fa data/GRCh38.primary_assembly.genome.fa --gtf data/gencode.v26.annotation.gtf --bedtools /Apps/bedtools/bin --cores 1 --splitsize 10000 --split $i;
done;

for i in `seq 1 27`;
do Rscript scripts/genome_predictions/predictBranchpointsGenomeAnnotation.R --fa data/GRCh37.p13.genome.fa --gtf data/gencode.v19.annotation.gtf --bedtools /Apps/bedtools/bin --cores 1 --splitsize 10000 --split $i;
done;

**4. Combine split predictions into single RData files
```
Rscript scripts/genome_predictions/processPredictions.R
```