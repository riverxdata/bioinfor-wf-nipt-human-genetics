<p align="center"> <h1 align="center">Bioinfor-wf-nipt-human-genetics</h1> </p> <p align="center"> <a href=""> <img height="300" src="./docs/assets/images/NIPT-Human-Genetics.png" alt="NIPT Human Genetics" align="center"> </a> </p>

## Overview
Bioinfor-wf-nipt-human-genetics is an end-to-end NIPT (Non-Invasive Prenatal Testing) GWAS (Genome-Wide Association Study) pipeline designed for processing fastq files. This pipeline facilitates the analysis of human genetics for NIPT applications and can be easily adapted for various research scenarios.

## Key Features
End-to-end workflow for processing fastq files into NIPT GWAS results.

Designed to handle more than 50 samples for real-world applications.

Collects necessary parameters for easy configuration.

Simplifies the creation of test profiles for easy setup and execution.

## Usages
### Install river-utils:
Install micromamba and river-utils
```bash
version="v1.2.0"
bash <(curl -Ls https://raw.githubusercontent.com/riverxdata/river-utils/${version}/install/setup.sh) $HOME $version
source ~/.river.sh
```

Check required softwares
```bash
which nextflow
which singularity
# examples
# /home/river/.river/images/micromamba/envs/river/bin/nextflow
# /home/river/.river/images/micromamba/envs/river/bin/singularity
```

### Setup
#### Get latest workflow
```bash
git clone https://github.com/riverxdata/bioinfor-wf-nipt-human-genetics
cd bioinfor-wf-nipt-human-genetics
```

#### Prepare the reference genome and database
```bash
# gcs 1KPG with known variants db and reference genomes
micromamba install conda-forge::google-cloud-sdk
mkdir -p data/gcs
gsutil -m cp -r \
  "gs://genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf" \
  "gs://genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx" \
  "gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz" \
  "gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi" \
  "gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz" \
  "gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" \
  "gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz" \
  "gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi" \
  "gs://genomics-public-data/resources/broad/hg38/v0/scattered_calling_intervals" \
  "gs://genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list" \
  .
# gmap files
cd ..
mkdir -p gmap && cd gmap
git clone https://github.com/odelaneau/GLIMPSE.git

# phasing db
cd ..
mkdir -p phasing && cd phasing
# example to run on small chromosomes
for i in {21..22}
do
    wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.vcf.gz
    wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.vcf.gz.tbi
done
```

#### Data
Open this on your browser, https://zenodo.org/records/13382182
Download and unzip in the data folder
```bash
data
├── gcs
├── gmap
├── gwas
├── phasing
└── samples
    ├── CL10001L0250.fq.gz
    ├── CL10002L0141.fq.gz
    ├── CL10003L0260.fq.gz
    ├── CL10004L0265.fq.gz
    ├── CL10005L0280.fq.gz
    ├── CL10006L0260.fq.gz
    ├── CL10007L0165.fq.gz
    ├── CL10008L0241.fq.gz
    ├── CL10009L0160.fq.gz
    ├── CL10010L0282.fq.gz
    ├── files-archive
    ├── fq.list
    └── sample.csv
```

### Execute the pipeline
```bash
nextflow run main.nf -profile docker -resume --outdir result
# Launching `main.nf` [insane_franklin] DSL2 - revision: cea5343a10

# ================== Pipeline Parameters =======================================================
# Input file:                 /home/river/bioinfor-wf-nipt/data/samples/fq.list
# Phenotype file:            /home/river/bioinfor-wf-nipt/data/gwas/phenotype.txt
# Reference genome:          /home/river/bioinfor-wf-nipt/data/gcs/Homo_sapiens_assembly38
# 1k Gold Indel DB:          /home/river/bioinfor-wf-nipt/data/gcs/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
# Known Assembly Indel DB:   /home/river/bioinfor-wf-nipt/data/gcs/Homo_sapiens_assembly38.known_indels.vcf.gz
# dbSNP DB:                  /home/river/bioinfor-wf-nipt/data/gcs/Homo_sapiens_assembly38.dbsnp138.vcf
# Genetic Maps:              /home/river/bioinfor-wf-nipt/data/gmap/GLIMPSE/maps/genetic_maps.b38
# Reference Panel DB:        /home/river/bioinfor-wf-nipt/data/phasing
# ==============================================================================================
# Workflow maintained by Thanh-Giang (River) Tan Nguyen
# For inquiries, contact: giangnguyen@riverxdata.com
```