# Anaquin

**Anaquin** is a C++/R bioinformatics framework for quantitative controls in next-generation sequencing experiments.

You may need to rename the `data` directory to `resources` for resource bundle.

## Credits

[Kallisto](https://github.com/pachterlab/kallisto) for indexing FASTA files.

## Compilation

A modern (supporting C++11) C++ compiler is required.

    make

## Overview

Sequins are synthetic spike-ins that are added to an RNA or DNA sample, and act as internal controls during next generation sequencing (NGS) experiments. Sequins can be used in RNA sequencing, human and cancer genomics, and metagenomics.

`Anaquin` is a software toolkit that can identify sequin reads within any NGS library (.FASTQ or BAM), calibrate coverage and dilution between samples, and issue detailed reports on the performance of sequins within the NGS library. Anaquin accepts standardized file formats (FASTQ, BAM, VCF, BED) allowing easy integration into standard NGS workflows and common bioinformatic tools.

For more information, please see www.sequinstandards.com
