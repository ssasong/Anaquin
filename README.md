# Anaquin

**Anaquin** is a C++/R bioinformatics framework for quantitative controls in next-generation sequencing experiments.

Repository for the R-package is hosted by Bioconductor and available at: https://github.com/sequinstandards/RAnaquin.

This is a **beta software** as we are trying to work with the bioinformatics community. Please send us your suggestions (e.g. what do you want Anaquin to do?).

You may need to rename the `data` directory to `resources` for resource bundle.

## Credits

[Kallisto](https://github.com/pachterlab/kallisto) for indexing FASTA files.

## Docker

Docker image for Anaquin is available at https://hub.docker.com/r/sequins/anaquin.

## Compilation

A modern (supporting C++11) C++ compiler is required.

    make

## Overview

Sequins are synthetic spike-ins that are added to an RNA or DNA sample, and act as internal controls during next generation sequencing (NGS) experiments. Sequins can be used in RNA sequencing, human and cancer genomics, and metagenomics.

`Anaquin` is a software toolkit that can identify sequin reads within any NGS library (.FASTQ or BAM), calibrate coverage and dilution between samples, and issue detailed reports on the performance of sequins within the NGS library. Anaquin accepts standardized file formats (FASTQ, BAM, VCF, BED) allowing easy integration into standard NGS workflows and common bioinformatic tools.

For more information, please see www.sequinstandards.com
