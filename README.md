# MGCplotter: Microbial Genome Circular plotter

![Python3](https://img.shields.io/badge/Language-Python3-steelblue)
![OS](https://img.shields.io/badge/OS-_Mac_|_Linux-steelblue)
![License](https://img.shields.io/badge/License-GPL3-steelblue)
[![Latest PyPI version](https://img.shields.io/pypi/v/mgcplotter.svg)](https://pypi.python.org/pypi/mgcplotter)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/mgcplotter.svg?color=green)](https://anaconda.org/bioconda/mgcplotter)
[![CI](https://github.com/moshi4/MGCplotter/actions/workflows/ci.yml/badge.svg)](https://github.com/moshi4/MGCplotter/actions/workflows/ci.yml)

## Currently Under Construction (Not Released)

## Overview

MGCplotter is easy-to-use circular layout plotting tool for microbial genome.
MGCplotter requires Genbank format genome file and implements following three main functions.

1. **`Plot microbial genome basic features`**  
  Basic Features = *Forward/Reverse CDS*, *rRNA*, *tRNA*, *GC content*, *GC skew*.
  MGCplotter can fully control plot result of feature's color/size/visibility by command options.  

2. **`Assign & Plot COG functional classification of reference genome CDS`**  
  Assign COG functional classification to reference species genome CDS using [COGclassifier](https://github.com/moshi4/COGclassifier).
  COG functional classification colors are used in plot result of forward/reverse CDS.

    <details>
    <summary>List of COG Functional Classification Color</summary>

    ![COG_definition_fig](https://github.com/moshi4/MGCplotter/blob/main/images/cog_definition_legend.png?raw=true)  

    </details>

3. **`Search & Plot conserved CDS between reference and query species`**  
  Conserved CDS of query genome relative to reference genome is searched by [MMseqs2](https://github.com/moshi4/COGclassifier) RBH method.
  Each query conserved CDS is plotted with gradient color based on CDS identity.

![MGCplotter_example_fig](https://github.com/moshi4/MGCplotter/blob/main/images/02_mycoplasma.png?raw=true)  
Fig.1: Plot result of Mycoplasma Gallisepticum genome  
Outer to inner tracks mean (1) Forward CDS (2) Reverse CDS (3) rRNA (4) tRNA (5) GC content (6) GC skew, respectively.
COG functional classification color is assigned to Forward/Reverse CDSs.

![MGCplotter_example_fig](https://github.com/moshi4/MGCplotter/blob/main/images/03_mycoplasma.png?raw=true)  
Fig.2: Add 3 query species conserved CDS track from Fig.1.
Conserved CDS of query genomes relative to reference genome is shown.

## Installation

MGCplotter is implemented in Python3.

**Install bioconda package:**

    conda install -c bioconda -c conda-forge mgcplotter

**Install PyPI stablepakcage:**

    pip install mgcplotter

**Use Docker (Docker Image):**

    docker pull moshi4/mgcplotter:latest

## Dependencies

- [Circos](http://circos.ca/)  
  Software package for visualizing data and information in a circular layout
- [COGclassifier](https://github.com/moshi4/COGclassifier)  
  A tool for classifying prokaryote protein sequences into COG functional category
- [MMseqs2](https://github.com/soedinglab/MMseqs2)  
  Ultra fast and sensitive sequence search and clustering suite

## Usage

### Options

    General Options:
      -r R, --ref_file R      Reference genbank file (*.gb|*.gbk|*.gbff)
      -o O, --outdir O        Output directory
      --query_files  [ ...]   Query fasta or genbank files (*.fa|*.faa|*.fasta|*.gb|*.gbk|*.gbff)
      --mmseqs_evalue         MMseqs e-value parameter (Default: 1e-05)
      --cog_evalue            COGclassifier e-value parameter (Default: 1e-02)
      -t , --thread_num       Threads number parameter (Default: MaxThreads - 1)
      -f, --force             Forcibly overwrite previous calculation result (Default: OFF)
      -v, --version           Print version information
      -h, --help              show this help message and exit

    Graph Size Options:
      --ticks_labelsize       Ticks label size (Default: 35)
      --forward_cds_r         Forward CDS track radius size (Default: 0.07)
      --reverse_cds_r         Reverse CDS track radius size (Default: 0.07)
      --rrna_r                rRNA track radius size (Default: 0.07)
      --trna_r                tRNA track radius size (Default: 0.07)
      --conserved_seq_r       Conserved seq track radius size (Default: 0.04)
      --gc_content_r          GC content track radius size (Default: 0.15)
      --gc_skew_r             GC skew track radius size (Default: 0.15)

    Graph Color Options:
      --assign_cog_color      Assign COG classification color to reference CDSs (Default: OFF)
      --cog_color_json        User-defined COG classification color json file
      --forward_cds_color     Forward CDS color (Default: 'red')
      --reverse_cds_color     Reverse CDS color (Default: 'blue')
      --rrna_color            rRNA color (Default: 'green')
      --trna_color            tRNA color (Default: 'magenta')
      --conserved_seq_color   Conserved sequence color (Default: 'chocolate')
      --gc_content_p_color    GC content color for positive value from average (Default: 'black')
      --gc_content_n_color    GC content color for negative value from average (Default: 'grey')
      --gc_skew_p_color       GC skew color for positive value (Default: 'olive')
      --gc_skew_n_color       GC skew color for negative value (Default: 'purple')

> :warning: If graph size option is set 0 value, target is not shown in figure.

## Output Contents

## Gallery

![MGCplotter_example_fig](https://github.com/moshi4/MGCplotter/blob/main/images/05_mycoplasma.png?raw=true)  

![MGCplotter_example_fig](https://github.com/moshi4/MGCplotter/blob/main/images/04_ecoli.png?raw=true)  
