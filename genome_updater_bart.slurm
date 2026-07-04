#!/bin/bash
#
#SBATCH --job-name=GenomeUpdater
#SBATCH -p lo_mem
#SBATCH -w bart
#SBATCH -o genome_updater.%j.out # STDOUT
#SBATCH -e genome_updater.%j.err # STDERR
#SBATCH --cpus-per-task=16

# Log start time
date

# Log the hostname
hostname

# Download GenBank viral genomes
../genome_updater.sh \
    -d "genbank" \
    -g "viral" \
    -c "all" \
    -l "complete genome,chromosome" \
    -f "genomic.fna.gz" \
    -o "genbank-viral" \
    -t 16 \
    -A "species:1" \
    -m -a -p

# Log end time
date
