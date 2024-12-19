**Background:**

This repository accompanies the manuscript *SARS-CoV-2 CoCoPUTs: an
Analysis of Sequences from Two of the Largest SARS-CoV-2 Databases
(GISAID and NCBI Virus) to understand Trends in Codon Usage and
Secondary Structure over a multi-year period*. It is intended to:

1.  Facilitate the initial parsing of genomic sequence data and metadata
    from both GISAID and NCBI, using nextclade to QC sequences

2.  Create a combined sequence resource from these two databases
    (de-duplicate the overlapping sequences between the two and then use
    all the resulting sequences for downstream analysis)

3.  Obtain one-sequence-per-variant using several different
    methodologies

4.  Calculate many different sequence composition metrics for each
    SARS-CoV-2 variant (codon usage, codon-pair usage, dinucleotide
    counts, etc) as well as free energy and secondary structure

5.  Use metadata to create time-series visualizations for the above
    metrics, weighting each variant by its abundance at given timepoints
    from early 2020 to mid-2024.

**Dependencies:**

The below libraries are required dependencies for the scripts in this
repository -- many are built in python packages or available on linux
command line. With the exception of forgi, all can be installed into a
single python3 conda environment. To run forgi without issues, it is
best to do so in a separate conda environment, following the
installation instructions on their documentation
(<https://viennarna.github.io/forgi/download.html>). Additionally, one
should download the forgi git repository
(<https://github.com/ViennaRNA/forgi>), and take note of the path of the
rnaConvert.py script within this directory.

1.  os

2.  Bio

3.  subprocess

4.  nextclade

5.  awk

6.  Dask

7.  Sqlite3

8.  Time

9.  Multiprocessing

10. Concurrent

11. Hashlib

12. dateutil

13. Collections

14. pandas

15. Levenshtein

16. Pango aliasor

17. Numpy

18. Re

19. Codonbias

20. Warnings

21. Matplotlib

22. seaborn

23. Statsmodels

24. halign

25. Forgi

26. Linearpartition

**Downloading Data:**

The following data should be downloaded:

1.  Nextclade dataset: In this study, the wuhan-hu-1 dataset was used
    (<https://github.com/nextstrain/nextclade_data/tree/release/data/nextstrain/sars-cov-2/wuhan-hu-1>

2.  SHAPE data from Manfredonia *et al*., available here (download the
    xml file):
    <https://www.incarnatolab.com/datasets/SARS_Manfredonia_2020.php>

3.  GISAID data -- must be downloaded from gisaid.org, after agreeing to
    their terms of use. You may need to reach out to gisaid support in
    order to have access to download the full sequence data and
    metadata. Specifically, after logging in, one should click
    "Downloads", followed by "FASTA" and "metadata" under "Download
    packages". Again, if these options do not appear to you, please
    reach out to GISAID support.

4.  NCBI Data -- First the datasets package must be installed -- this
    can be done via conda. The full sequence data and metadata can be
    then downloaded with: **datasets download virus genome taxon
    SARS-CoV-2 \--filename sars_cov_2.zip**

**Running python files:**

All python files are designed to be run in the order of the number in
their filename. If two files share the same order, they can be run in
any order.

-   run_nextclade_on_large_fastas_1_final.py

This script runs nextclade on the (very large) raw fasta files from both
GISAID and NCBI. This enables the later QC steps.

-   parse_gisaid_accessions_and_merge_metadata_2_final.py

This script merges relevant metadata for each GISAID sequence into the
sequence header, and also QC's the sequences based on information from
the previous nextclade step. It also separates the sequences into
individual fasta files for each PANGO variant.

-   parse_ncbi_accessions_and_merge_metadata_2_final.py

Identical to #2 but for ncbi sequences. Can be run before or after #2.

-   combine_gisaid_and_ncbi_3_final.py

This creates the combined dataset of sequences from gisaid and ncbi
(deduplicating the overlap), using a combination of the sequence data
itself and the associated metadata.

-   take_unique_sequences_and_align_4_final.py

this script reduces each fasta file to only the unique sequences (to
enable easier/faster alignment), while still preserving the count of
these sequences. It then aligns them using Halign3.

-   build_representative_consensus_ancestral_seqs_5_final.py

This script uses three different methodologies to reduce the dataset to
1-sequence-per-variant. These are elaborated in our publication.

-   align_final_sequences_and_add_ancestry_6_final.py

This performs a final alignment on just the reduced datasets. It also
adds lineage ancestry information with the Aliasor package.

-   prepare_linear_partition_7_final.py

This script prepares the experimental shape data for use on the many
different sars-cov-2 variants.

-   run_linear_partition_8_final.py

This script runs linearpartition and forgi on all of our sars-cov-2
variants, generating free energy data, secondary structure data and the
corresponding forgi strings.

-   calculate_stats_9_final.py

This script calculates all the sequence composition statistics for all
of our sars-cov-2 variants

-   plotting_10_final.py

This script plots all of our time-series analyses and most other
visualizations featured in our publication. It also creates all of the
supplemental files.
