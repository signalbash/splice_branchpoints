# Download annotations

## Gencode annotations
*Gencodev24*
```
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/GRCh38.p5.genome.fa.gz
gunzip gencode.v24.annotation.gtf.gz
gunzip GRCh38.p5.genome.fa.gz
```

*Gencodev12*
```
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_12/gencode.v12.annotation.gtf.gz
#Gencode has removed the non-patched v12 version
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz
gunzip gencode.v12.annotation.gtf.gz
gunzip GRCh37.p13.genome.fa.gz
```

*Gencodev19*
```
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip gencode.v19.annotation.gtf.gz
```

*GencodevM9*
```
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gtf.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M10/GRCm38.p4.genome.fa.gz
gunzip gencode.vM10.annotation.gtf.gz
gunzip GRCm38.p4.genome.fa.gz
```

Format GTF as exon attribute table
```
../../scripts/genome_predictions/format_exon_annotations.sh gencode.v24.annotation.gtf gencode.v24.exons.txt
../../scripts/genome_predictions/format_exon_annotations.sh gencode.v19.annotation.gtf gencode.v19.exons.txt
Rscript ../../scripts/genome_predictions/gtf_to_exon_annotation.R -a gencode.v12.annotation.gtf -t type
../../scripts/genome_predictions/format_exon_annotations.sh gencode.vM10.annotation.gtf gencode.vM10.exons.txt
```

## Ensembl annotations
*Danio rerio*
```
wget ftp://ftp.ensembl.org/pub/release-85/gtf/danio_rerio/Danio_rerio.GRCz10.85.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-85/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz
gunzip Danio_rerio.GRCz10.85.gtf.gz
gunzip Danio_rerio.GRCz10.dna.toplevel.fa.gz
```

*Xenopus_tropicalis*
```
wget ftp://ftp.ensembl.org/pub/release-85/fasta/xenopus_tropicalis/dna/Xenopus_tropicalis.JGI_4.2.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-85/gtf/xenopus_tropicalis/Xenopus_tropicalis.JGI_4.2.85.gtf.gz
gunzip Xenopus_tropicalis.JGI_4.2.dna.toplevel.fa.gz
gunzip Xenopus_tropicalis.JGI_4.2.85.gtf.gz
```

*Drosophila_melanogaster*
```
wget ftp://ftp.ensembl.org/pub/release-85/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-85/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.85.gtf.gz
gunzip Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz
gunzip Drosophila_melanogaster.BDGP6.85.gtf.gz
```

*Gallus_gallus*
```
wget ftp://ftp.ensembl.org/pub/release-85/fasta/gallus_gallus/dna/Gallus_gallus.Galgal4.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-85/gtf/gallus_gallus/Gallus_gallus.Galgal4.85.gtf.gz
gunzip Gallus_gallus.Galgal4.dna.toplevel.fa.gz
gunzip Gallus_gallus.Galgal4.85.gtf.gz
```

Format GTF as exon attribute table
```
Rscript ../../scripts/genome_predictions/gtf_to_exon_annotation.R -a Danio_rerio.GRCz10.85.gtf
Rscript ../../scripts/genome_predictions/gtf_to_exon_annotation.R -a Xenopus_tropicalis.JGI_4.2.85.gtf
Rscript ../../scripts/genome_predictions/gtf_to_exon_annotation.R -a Drosophila_melanogaster.BDGP6.85.gtf
Rscript ../../scripts/genome_predictions/gtf_to_exon_annotation.R -a Gallus_gallus.Galgal4.85.gtf
```

#Prediction of intronic branchpoints from exon annotations

```
Rscript ../../scripts/genome_predictions/predict_branchpoints_in_genome.R --split 1 --splitsize 99999999 -c 1 -b /apps/bedtools/2.19.1/bedtools -fa GRCh38.p5.genome.fa -e gencode.v24.exons.txt -p gencode.v24
Rscript ../../scripts/genome_predictions/predict_branchpoints_in_genome.R --split 1 --splitsize 99999999 -c 1 -b /apps/bedtools/2.19.1/bedtools -fa GRCh37.p13.genome.fa -e gencode.v12.exons.txt -p gencode.v12
Rscript ../../scripts/genome_predictions/predict_branchpoints_in_genome.R --split 1 --splitsize 99999999 -c 1 -b /apps/bedtools/2.19.1/bedtools -fa GRCm38.p4.genome.fa -e gencode.vM9.exons.txt -p gencode.vM9

Rscript ../../scripts/genome_predictions/predict_branchpoints_in_genome.R --split 1 --splitsize 99999999 -c 1 -b /apps/bedtools/2.19.1/bedtools -fa Danio_rerio.GRCz10.dna.toplevel.fa -e Danio_rerio.GRCz10.85.exons.txt -p Danio_rerio
Rscript ../../scripts/genome_predictions/predict_branchpoints_in_genome.R --split 1 --splitsize 99999999 -c 1 -b /apps/bedtools/2.19.1/bedtools -fa Xenopus_tropicalis.JGI_4.2.dna.toplevel.fa -e Xenopus_tropicalis.JGI_4.2.85.exons.txt -p Xenopus_tropicalis
Rscript ../../scripts/genome_predictions/predict_branchpoints_in_genome.R --split 1 --splitsize 99999999 -c 1 -b /apps/bedtools/2.19.1/bedtools -fa Drosophila_melanogaster.BDGP6.dna.toplevel.fa -e Drosophila_melanogaster.BDGP6.85.exons.txt -p Drosophila_melanogaster
Rscript ../../scripts/genome_predictions/predict_branchpoints_in_genome.R --split 1 --splitsize 99999999 -c 1 -b /apps/bedtools/2.19.1/bedtools -fa Gallus_gallus.Galgal4.dna.toplevel.fa -e Gallus_gallus.Galgal4.85.exons.txt -p Gallus_gallus
```

*Process non-human species*
`Rscript ../../scripts/genome_predictions/process_predictions_nonhuman.R`

*Process human branchpoint predcitions*
`Rscript ../../scripts/genome_predicitons/process_predictions.R`




