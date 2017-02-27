# SL-Seq
This repository contains all neccesairy tools to analyse data generated with Spliced-Leader Sequencing (SL-Seq).

The tool SLSeq-GFFconvert.py converts a standard GFF file to a format that is suited for SL-seq analysis. Briefly, the scripts expands the provided annotations to include the 5' UTRs of the genes, where a significant proportion of the SL-seq reads map. The tool can be excecuted by calling 'python SLSeq-GFFconvert.py' from the command line. To see all options, type 'python SLSeq-GFFconvert.py --help'.

Full analysis scripts from raw sequence data (FASTQ) to count tables are available for both paired-end (SL-SeqPipeline_PairedEnd.sh) and single-end sequencing data (SL-SeqPipeline_SingleEnd.sh). Please read carefully the comments contained within those scripts in order to install the right versions of every tool and provide the necessary annotation files.
