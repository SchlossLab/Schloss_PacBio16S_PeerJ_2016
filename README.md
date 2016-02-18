Sequencing 16S rRNA gene fragments using the PacBio SMRT DNA sequencing system
=======

This is the repository for the manuscript "Sequencing 16S rRNA gene fragments
using the PacBio SMRT DNA sequencing system" written by Patrick D. Schloss,
Matthew L. Jenior, Charles C. Koumpouras, Sarah L. Westcott, and Sarah K.
Highlander. The raw data can be obtained from the Sequence Read Archive at NCBI
under accession  SRP051686, which are associated with BioProject PRJNA271568.
This project relies heavily on an earlier version of this project that used an
older version of the chemistry. Because that project was published as a preprint
and the current project represents a significant shift in structure, we have
created this repository. You can find the older version at
https://github.com/SchlossLab/Schloss_PacBio16S_PeerJ_2015.


Overview
--------

    project
    |- README          # the top level description of content
    |
    |- doc/            # documentation for the study
    |  |- notebook/    # preliminary analyses (dead branches of analysis)
    |  +- paper/       # manuscript(s), whether generated or not
    |
    |- data            # raw and primary data, are not changed once created
    |  |- references/  # reference files to be used in analysis
    |  |- raw_june/    # raw data from june sequencing
    |  |- raw_october/ # raw data from october sequencing
    |  |- mothur_june/ # processing of june data with mothur
    |  |- mothur_october # processing of october data with mothur
    |  +- process/     # cleaned data, will not be altered once created;
    |                  # will be committed to repo
    |
    |- code/           # any programmatic code
    |  |
    |  +- smrtanalysis/smrtcmds/bin/ #PacBio software. exe's are here
	|
	|- results         # all output from workflows and analyses
    |  |- tables/      # text version of tables to be rendered with kable in R
    |  |- figures/     # graphs, likely designated for manuscript figures
    |  +- pictures/    # diagrams, images, and other non-graph graphics
    |
    |- scratch/        # temporary files that can be safely deleted or lost
    |
    |- Schloss_PacBio16S_PeerJ_2016.Rmd	# Rmarkdown for this study
	|
	|- submission		# stuff to be submitted
	|	|- Schloss_PacBio16S_PeerJ_2016.md		# Markdown version of the *Rmd file
	|	|- Schloss_PacBio16S_PeerJ_2016.pdf		# pdf version of the *Rmd file
	|	|- various figures and tables
	|	+- rebuttal_v1.pdf						# response to reviewers' comments
    |
    +- Makefile        # executable Makefile for this study, if applicable
