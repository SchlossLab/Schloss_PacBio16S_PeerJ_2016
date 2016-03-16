REFS = data/references
FIGS = submission
TABLE = results/tables
PROC = data/process

# utility function to print various variables. For example, running the
# following at the command line:
#
#	make print-BAM
#
# will generate:
#	BAM=data/raw_june/V1V3_0001.bam data/raw_june/V1V3_0002.bam ...
print-%:
	@echo '$*=$($*)'



################################################################################
#
# Part 1: Get the references
#
# We will need several reference files to complete the analysis: the mock
# community sequences, the SILVA reference alignment, and the SILVA, RDP, and
# greengenes reference taxonomies.
#
################################################################################



# We want the latest greatest reference alignment and the SILVA reference
# alignment is the best reference alignment on the market. This version is from
# v123 and described at http://blog.mothur.org/2014/08/08/SILVA-v119-reference-files/.
# We will use the SEED and NR versions of the database, which contain 12,083 and 152,308
# bacterial sequences. This also contains the reference taxonomy. We will limit
# the databases to only include bacterial sequences.

$(REFS)/silva.seed.align :
	wget -N http://mothur.org/w/images/1/15/Silva.seed_v123.tgz
	tar xvzf Silva.seed_v123.tgz silva.seed_v123.align silva.seed_v123.tax
	mothur "#get.lineage(fasta=silva.seed_v123.align, taxonomy=silva.seed_v123.tax, taxon=Bacteria);degap.seqs(fasta=silva.seed_v123.pick.align, processors=8)"
	mv silva.seed_v123.pick.align $(REFS)/silva.seed.align
	rm Silva.seed_v123.tgz silva.seed_v123.*


$(REFS)/silva.nr.% :
	wget -N http://mothur.org/w/images/b/be/Silva.nr_v123.tgz
	tar xvzf Silva.nr_v123.tgz silva.nr_v123.align silva.nr_v123.tax
	mothur "#get.lineage(fasta=silva.nr_v123.align, taxonomy=silva.nr_v123.tax, taxon=Bacteria);degap.seqs(fasta=silva.nr_v123.pick.align, processors=8)"
	mv silva.nr_v123.pick.tax $(REFS)/silva.nr.tax
	mv silva.nr_v123.pick.ng.fasta $(REFS)/silva.nr.fasta
	rm Silva.nr_v123.tgz silva.nr_v123.*


# We also want the greengenes reference taxonomy. This version is from the
# greengenes v13_8_99 and is described at
# http://blog.mothur.org/2014/08/12/greengenes-v13_8_99-reference-files/

$(REFS)/gg_13_8_99.% :
	wget -N http://www.mothur.org/w/images/6/68/Gg_13_8_99.taxonomy.tgz
	tar xvzf Gg_13_8_99.taxonomy.tgz gg_13_8_99.fasta gg_13_8_99.gg.tax
	mv gg_13_8_99.* $(REFS)/
	rm Gg_13_8_99.taxonomy.tgz


# Next, we want the RDP reference taxonomy. The current version is v10 and we
# use a "special" pds version of the database files, which are described at
# http://blog.mothur.org/2014/10/28/RDP-v10-reference-files/

$(REFS)/trainset14_032015.% :
	wget -N http://www.mothur.org/w/images/8/88/Trainset14_032015.pds.tgz
	tar xvzf Trainset14_032015.pds.tgz trainset14_032015.pds/trainset14_032015.pds.*
	mv trainset14_032015.pds/* $(REFS)/
	rmdir trainset14_032015.pds
	rm Trainset14_032015.pds.tgz


# Now, we want to align the mock community reference sequences to our newly
# created silva.seed.align file...

$(REFS)/HMP_MOCK.% : $(REFS)/silva.seed.align
	wget --no-check-certificate -N -P $(REFS) https://raw.githubusercontent.com/SchlossLab/Kozich_MiSeqSOP_AEM_2013/master/data/references/HMP_MOCK.fasta
	mothur "#align.seqs(fasta=$(REFS)/HMP_MOCK.fasta, reference=$^)"
	rm $@.report



################################################################################
#
# Part 2: BAM to fasta/qual files
#
# The goal of this part was to take the PacBio subreads files and convert them
# into fasta and qual files. In some cases there were multiple movies per region
# per data generation event (e.g. june and  october of 2015). The project root
# directory has several folders that are useful for this makefile, namely data/
# and code/. Within data/ are raw_june/, mothur_june, raw_october, raw_june.
# These correspond to data dumps that we received from pacbio in the specified
# month. the october data was for the longer regions with longer movies.
#
################################################################################


# will replace this with whatever we pull down from SRA. this step was done to
# combie the hd5 files into subreads.bam files for the june raw data. probably
# not the way one would want to run these commands. would be better to split
# them into individual pbs scripts

generate_bam_subread_files : code/run_bax2bam.sh
	bash code/run_bax2bam.sh



# this assumes that the pbccs package is installed in code/pbcss and uses the
# default. probably not the way one would want to run these commands. would be
# better to split them into individual pbs scripts. it also assumes that all of
# the raw data are stored like: data/raw_<month>/<region>_<movie>.subreads.bam

SR_BAM = $(wildcard data/raw*/*.subreads.bam)
BAM = $(subst .subreads.bam,.bam,$(SR_BAM))

.SECONDEXPANSION:
$(BAM) : $$(subst .bam,.subreads.bam,$$@)
	code/smrtanalysis/smrtcmds/bin/ccs $@ $^ --numThreads=8



# each of these *.bam files should have a *.scrap.bam and *.ccs_report.csv file
# with them. We'll extract the fastq files from each ccs-processed bam file
# using samtools:

BAM_FASTQ = $(subst .subreads.bam,.bam.fastq,$(SR_BAM))

.SECONDEXPANSION:
$(BAM_FASTQ) : $$(subst .fastq,,$$@)
	samtools bam2fq $^ > $@



# Now we need to extract the number of passes and predicted error rate from the
# bam file:

BAM_CCS_STATS = $(subst .subreads.bam,.bam.ccs_stats,$(SR_BAM))

.SECONDEXPANSION:
$(BAM_CCS_STATS) : $$(subst .ccs_stats,,$$@)
	samtools view $^ | cut -f 1,17,16 | sed -s "s/np:i://" | sed -s "s/rq:f://" > $@



# We want to pool the *.fastq and *.ccs_stats files for each region within each
# data and we'll drop them in the appropriate data/mothur_* folder

FASTQ = $(sort $(subst raw,mothur,$(addsuffix .fastq,$(basename $(subst _0,.,$(subst .bam.fastq,,$(BAM_FASTQ)))))))

$(FASTQ) : $$(subst mothur,raw,$$(subst .fastq,_*.bam.fastq,$$@))
	cat $^ > $@


CCS_STATS = $(sort $(subst raw,mothur,$(addsuffix .ccs_stats,$(basename $(subst _0,.,$(subst .bam.fastq,,$(BAM_FASTQ)))))))

$(CCS_STATS) : $$(subst mothur,raw,$$(subst .ccs_stats,_*.bam.ccs_stats,$$@))
	cat $^ > $@



# Finally, let's extract the fasta and quality score data from the fastq files

FASTA_QUAL = $(subst fastq,fasta,$(FASTQ)) $(subst fastq,qual,$(FASTQ))
$(FASTA_QUAL) : $$(addsuffix .fastq,$$(basename $$@))
	mothur "#fastq.info(fastq=$^, pacbio=T)"




################################################################################
#
# Part 3: Processing to separate reads by library
#
# Now we're all set to run some mothur commands. Since each file is a mixture of
# our mock community and data from soil, human feces, and mouse feces, we need
# to split the fasta and qual files by barcode. Let's initially be generous and
# allow for 2 mismatches to each barcode and 4 mismatches to each primer. To
# keep things simple, we'll concatenate the three mock community fasta, quality
# score, and groups files.
#
################################################################################


SAMPLES = mock soil human mouse

SAMPLE_FASTA = $(foreach S,$(SAMPLES),$(foreach B,$(sort $(basename $(FASTA_QUAL))),$B.$S.fasta))
SAMPLE_QUAL = $(subst fasta,qual,$(SAMPLE_FASTA))
SAMPLE_GROUPS = $(subst fasta,groups,$(SAMPLE_FASTA))

$(SAMPLE_FASTA) $(SAMPLE_QUAL) $(SAMPLE_GROUPS) : $$(addsuffix .fasta,$$(basename $$(basename $$@))) $$(addsuffix .qual,$$(basename $$(basename $$@))) $(REFS)/pacbio.oligos
	$(eval RAW_FASTA = $(word 1, $^))
	$(eval RAW_QUAL = $(word 2, $^))
	$(eval OLIGOS = $(word 3, $^))
	$(eval REGION = $(subst october/,,$(subst june/,,$(patsubst data/mothur_%.fasta,%,$(RAW_FASTA)))))
	mothur "#trim.seqs(fasta=$(RAW_FASTA), qfile=$(RAW_QUAL), oligos=$(OLIGOS), checkorient=T, pdiffs=4, bdiffs=2, allfiles=T, processors=8)"
	for S in $(SAMPLES) ; do \
		cat $(patsubst %.fasta,%.$$S*.$(REGION).fasta,$(RAW_FASTA)) > $(patsubst %.fasta,%.$$S.fasta,$(RAW_FASTA)); \
		cat $(patsubst %.fasta,%.$$S*.$(REGION).qual,$(RAW_FASTA)) > $(patsubst %.fasta,%.$$S.qual,$(RAW_FASTA)); \
		cat $(patsubst %.fasta,%.$$S*.$(REGION).groups,$(RAW_FASTA)) > $(patsubst %.fasta,%.$$S.groups,$(RAW_FASTA)); \
		rm $(patsubst %.fasta,%.$$S*.V**,$(RAW_FASTA)); \
	done
	rm $(patsubst %.fasta,%.trim.*,$(RAW_FASTA));
	rm $(patsubst %.fasta,%.scrap.*,$(RAW_FASTA));
	rm $(patsubst %.fasta,%.groups,$(RAW_FASTA));




################################################################################
#
# Part 4: Processing mock community data
#
# Here we'll work with the mock community samples to get a sense of their error
# rate
#
################################################################################

MOCK_FASTA = $(addsuffix .mock.fasta,$(sort $(basename $(basename $(SAMPLE_FASTA)))))
MOCK_QUAL = $(subst fasta,qual,$(MOCK_FASTA))

# The number of mismatches to the barcodes and primers is on the header line for
# each sequence in the trim file and is extracted here to a *.mismatches file

MISMATCH = $(subst fasta,mismatches,$(MOCK_FASTA))

$(MISMATCH) : $$(subst mismatches,fasta,$$@)
	grep ">" $^ | cut -c 2- > $@



# We are now ready to calculate the error rate for the various regions using the
# mock community and the HMP_MOCK sequence data. We will use mothur to align the
# sequences to HMP_MOCK.align, determine the start and end positions of the
# alignment, and calculate the error rate.

ALIGN_SUMMARY = $(subst fasta,filter.summary,$(MOCK_FASTA))
ERROR_SUMMARY = $(subst	summary,error.summary,$(ALIGN_SUMMARY))
ERROR_QUALITY = $(subst summary,error.quality,$(ALIGN_SUMMARY))
MOCK_ANALYSIS = $(ALIGN_SUMMARY) $(ERROR_SUMMARY) $(ERROR_QUALITY)

$(MOCK_ANALYSIS) : $$(subst .filter,.fasta,$$(subst .error,,$$(basename $$@)))\
		$$(subst .filter,.qual,$$(subst .error,,$$(basename $$@)))\
		$$(REFS)/HMP_MOCK.align
	$(eval FASTA = $(word 1, $^))
	$(eval QUAL = $(word 2, $^))
	$(eval MOCK = $(word 3, $^))
	$(eval STUB = $(subst .fasta,,$(FASTA)))
	cp $(MOCK) $(STUB).HMP_MOCK.align
	mothur "#align.seqs(fasta=$(FASTA), reference=$(STUB).HMP_MOCK.align, processors=8);\
	    filter.seqs(fasta=$(STUB).align-$(STUB).HMP_MOCK.align, vertical=T);\
	    summary.seqs();\
	    seq.error(fasta=$(STUB).filter.fasta, reference=$(STUB).HMP_MOCK.filter.fasta, report=$(STUB).align.report, qfile=$(QUAL));"
	rm $(STUB).align*
	rm $(STUB).filter.fasta
	rm $(STUB).filter.error.seq
	rm $(STUB).filter.error.chimera
	rm $(STUB).filter.error.qual.*
	rm $(STUB).filter.error.seq.forward
	rm $(STUB).filter.error.seq.*
	rm $(STUB).filter.error.count
	rm $(STUB).filter.error.ref
	rm $(STUB).HMP_MOCK*



# Now we need to synthesize the various output files into a report that we can
# use to identify the best parameter settings to minimize the error rate.

MOCK_REPORT = $(subst filter.summary,report,$(ALIGN_SUMMARY))

$(MOCK_REPORT) : code/consolidate_data.R\
				$$(subst report,filter.error.summary,$$@)\
				$$(subst mock.report,ccs_stats,$$@)\
				$$(subst report,mismatches,$$@)\
				$$(subst report,filter.summary,$$@)\
				$$(subst report,qual,$$@)
	R -e "source('$<');generate_report('$@')"




REGIONS = V4 V1V3 V1V5 V1V6 V1V9 V3V5

POOL_FASTA = $(sort $(subst mothur_june,mothur_pool,$(subst mothur_october,mothur_pool,$(SAMPLE_FASTA))))

$(POOL_FASTA) : $(SAMPLE_FASTA)
	$(eval F = $(notdir $@))
	cat $(filter %$F, $(SAMPLE_FASTA)) > $@


POOL_CCS_STATS = $(sort $(subst mothur_june,mothur_pool,$(subst mothur_october,mothur_pool,$(CCS_STATS))))
$(POOL_CCS_STATS) : $(CCS_STATS)
	$(eval F = $(notdir $@))
	cat $(filter %$F, $(CCS_STATS)) > $@



PRED_ERROR_SCREEN = $(subst fasta,screen.fasta,$(POOL_FASTA))
$(PRED_ERROR_SCREEN) : $$(addsuffix .ccs_stats, $$(basename $$(basename $$(basename $$@))))\
			$$(addsuffix .fasta, $$(basename $$(basename $$@)))\
			code/pred_error_screen.R
	$(eval CCS_STATS = $(word 1, $^))
	$(eval FASTA = $(word 2, $^))
	R -e 'source("code/pred_error_screen.R");pred_error_screen("$(FASTA)", "$(CCS_STATS)")'




UNIQUE_FASTA = $(subst fasta,unique.good.filter.unique.fasta,$(PRED_ERROR_SCREEN))
UNIQUE_NAMES = $(subst unique.fasta,names,$(UNIQUE_FASTA))
PRECLUSTER_FASTA = $(subst fasta,precluster.fasta,$(UNIQUE_FASTA))
PRECLUSTER_NAME = $(subst fasta,names,$(PRECLUSTER_FASTA))

UNIQUE_FILES = $(UNIQUE_FASTA) $(UNIQUE_NAMES)
$(UNIQUE_FILES) : $$(subst unique.good.filter.names,fasta,$$(subst unique.fasta,names,$$@)) code/get_unique_pc_fasta.sh
	bash code/get_unique_pc_fasta.sh $<


UNIQUE_ERROR = $(subst fasta,error.summary,$(filter %.mock.screen.unique.good.filter.unique.fasta,$(UNIQUE_FASTA)))
$(UNIQUE_ERROR) : $$(subst error.summary,fasta,$$@) $$(subst unique.error.summary,names,$$@) $(REFS)/HMP_MOCK.fasta
	$(eval F=$(word 1,$^))

	$(eval N=$(word 2,$^))
	mothur "#seq.error(fasta=$F, name=$N, reference=$(REFS)/HMP_MOCK.fasta, processors=8, aligned=F)"




PRECLUSTER_FILES = $(PRECLUSTER_FASTA) $(PRECLUSTER_NAMES)
$(PRECLUSTER_FILES) : $$(subst unique.good.filter.unique.precluster,fasta,$$(basename $$@)) code/get_unique_pc_fasta.sh
	bash code/get_unique_pc_fasta.sh $<




PRECLUSTER_ERROR = $(subst fasta,error.summary,$(filter %.mock.screen.unique.good.filter.unique.precluster.fasta,$(PRECLUSTER_FASTA)))
$(PRECLUSTER_ERROR) : $$(subst error.summary,fasta,$$@) $$(subst error.summary,names,$$@) $(REFS)/HMP_MOCK.fasta
	$(eval F=$(word 1,$^))
	$(eval N=$(word 2,$^))
	mothur "#seq.error(fasta=$F, name=$N, reference=$(REFS)/HMP_MOCK.fasta, processors=8, aligned=F)"




CHIMERA_FASTA = $(subst fasta,pick.fasta,$(PRECLUSTER_FASTA))
CHIMERA_NAMES = $(subst names,pick.names,$(PRECLUSTER_NAMES))
CHIMERA = $(CHIMERA_FASTA) $(CHIMERA_NAMES)
$(CHIMERA) : $$(addsuffix .fasta,$$(basename $$(basename $$@))) $$(addsuffix .names,$$(basename $$(basename $$@)))
	$(eval F=$(word 1,$^))
	$(eval N=$(word 2,$^))
	mothur "#chimera.uchime(fasta=$F, name=$N);remove.seqs(fasta=current, name=current, accnos=current)"



LIST_FILES = $(subst fasta,an.list,$(CHIMERA_FASTA))
$(LIST_FILES) : $$(subst an.list,fasta,$$@) $$(subst an.list,names,$$@)
	$(eval F = $(word 1, $^))
	$(eval N = $(word 2, $^))
	$(eval S = $(basename $N))
	mothur "#dist.seqs(fasta=$F, cutoff=0.15, processors=8);cluster(name=$N);"
	rm $S.dist
	rm $S.an.sabund
	rm $S.an.rabund



PERFECT_LIST = $(subst error.summary,perfect.an.list,$(PRECLUSTER_ERROR))
$(PERFECT_LIST) : $$(subst perfect.an.list,error.summary,$$@) $$(subst perfect.an.list,fasta,$$@) $$(subst perfect.an.list,names,$$@)
	$(eval E=$(word 1,$^))
	$(eval F=$(word 2,$^))
	$(eval N=$(word 3,$^))
	$(eval S=$(basename $F))
	grep "2$$" $E | cut -f 1 > $S.perfect.accnos
	cp $F $S.perfect.fasta
	cp $N $S.perfect.names
	mothur "#remove.seqs(fasta=$S.perfect.fasta, name=$S.perfect.names, accnos=$S.perfect.accnos);dist.seqs(cutoff=0.15, processors=8);cluster(name=current)"
	mv $S.perfect.pick.an.list $S.perfect.an.list
	rm $S.perfect.pick.*abund
	rm $S.perfect.pick.dist
	rm $S.perfect.*names
	rm $S.perfect.*fasta
	rm $S.perfect.accnos



UCHIME_SOBS = $(subst list,ave-std.summary,$(LIST_FILES))
$(UCHIME_SOBS) : $$(subst ave-std.summary,list,$$@)
	mothur "#summary.single(list=$^, label=0.03, subsample=1000,calc=nseqs-sobs-coverage)"
	touch $@


NOCHIM_SOBS = $(subst list,ave-std.summary,$(PERFECT_LIST))
$(NOCHIM_SOBS) : $$(subst ave-std.summary,list,$$@)
	mothur "#summary.single(list=$^, label=0.03, subsample=1000,calc=nseqs-sobs-coverage)"
	touch $@


MOCK_REGION_ALIGN = $(foreach R,$(REGIONS),data/mothur_pool/HMP_MOCK.$R.align)
$(MOCK_REGION_ALIGN) : $(REFS)/HMP_MOCK.align code/get_mock_region.sh
	$(eval REGION = $(subst .,,$(suffix $(basename $@))))
	bash code/get_mock_region.sh $(REGION)


MOCK_REGION_FASTA = $(subst align,fasta,$(MOCK_REGION_ALIGN))
$(MOCK_REGION_FASTA) : $$(subst fasta,align,$$@)
	mothur "#degap.seqs(fasta=$^)"
	mv $(subst align,ng.fasta,$^) $@


NOERROR_SOBS = $(foreach R,$(REGIONS),data/mothur_pool/HMP_MOCK.$R.pick.phylip.an.summary)
data/mothur_pool/HMP_MOCK.%.pick.phylip.an.summary : data/mothur_pool/%.mock.screen.unique.good.filter.unique.precluster.error.summary data/mothur_pool/HMP_MOCK.%.align
	$(eval E = $(word 1, $^))
	grep "1$$" $E | cut -f 2 | sort | uniq > data/mothur_pool/$*.mock.screen.unique.good.filter.unique.precluster.ref.accnos
	mothur "#get.seqs(fasta=data/mothur_pool/HMP_MOCK.$*.align, accnos=data/mothur_pool/$*.mock.screen.unique.good.filter.unique.precluster.ref.accnos);dist.seqs(cutoff=0.15, output=lt); cluster(phylip=current); summary.single(label=0.03, calc=sobs)"
	rm data/mothur_pool/HMP_MOCK.$*.pick.phylip.an.sabund
	rm data/mothur_pool/HMP_MOCK.$*.pick.phylip.an.rabund
	rm data/mothur_pool/HMP_MOCK.$*.pick.phylip.an.list
	rm data/mothur_pool/HMP_MOCK.$*.pick.phylip.dist
	rm data/mothur_pool/HMP_MOCK.$*.pick.align
	rm data/mothur_pool/$*.mock.screen.unique.good.filter.unique.precluster.ref.accnos




RDP = $(subst fasta,pds.wang.taxonomy,$(CHIMERA_FASTA) $(MOCK_REGION_FASTA))
$(RDP) : $$(subst pds.wang.taxonomy,fasta,$$@) data/references/trainset10_082014.pds.fasta data/references/trainset10_082014.pds.tax
	mothur "#classify.seqs(fasta=$<,reference=data/references/trainset10_082014.pds.fasta, taxonomy=data/references/trainset10_082014.pds.tax, cutoff=80, processors=8);"
	rm $(subst taxonomy,tax.summary,$@)
	#keep: data/mothur_pool/*.pds.wang.taxonomy


GG = $(subst fasta,gg.wang.taxonomy,$(CHIMERA_FASTA) $(MOCK_REGION_FASTA))
$(GG) : $$(subst gg.wang.taxonomy,fasta,$$@) data/references/gg_13_8_99.fasta data/references/gg_13_8_99.gg.tax
	mothur "#classify.seqs(fasta=$<, reference=data/references/gg_13_8_99.fasta, taxonomy=data/references/gg_13_8_99.gg.tax, cutoff=80, processors=8)"
	rm $(subst taxonomy,tax.summary,$@)
	#keep: data/mothur_pool/*.gg.wang.taxonomy


SILVA = $(subst fasta,nr.wang.taxonomy,$(CHIMERA_FASTA) $(MOCK_REGION_FASTA))
$(SILVA) : $$(subst nr.wang.taxonomy,fasta,$$@) data/references/silva.nr.fasta data/references/silva.nr.tax
	mothur "#classify.seqs(fasta=$<, reference=data/references/silva.nr.fasta, taxonomy=data/references/silva.nr.tax, cutoff=80, processors=8)"
	rm $(subst taxonomy,tax.summary,$@)
	#keep: data/mothur_pool/*.nr.wang.taxonomy



# Let's pool all of the mock report files together and toss the output into the
# $(PROC) folder

$(PROC)/mock.error.report : $(MOCK_REPORT)
	$(eval FIRST = $(word 1, $^))
	head -n 1 $(FIRST) > $@
	for FILE in $^; do tail -n +2 $$FILE >> $@; done


# Let's pool the *.error.quality files into a single file to see whether there are
# relationships between quality scores and specific error types

$(PROC)/mock.quality.report : $(ERROR_QUALITY)
	R -e "source('code/pool_error_quality.R'); pool('$^')"


# Let's compress the final report files to have a smaller version that we can put
# in the repo

%.gz : %
	gzip < $^ > $@



$(PROC)/error_profile.json : code/get_error_profile.R \
								$(PROC)/mock.quality.report\
								$(PROC)/mock.error.report
	R -e "source('code/get_error_profile.R')"


$(PROC)/error_summary.tsv : code/get_error_rate_table.R\
								$(PROC)/mock.error.report\
								$(PRECLUSTER_ERROR)
	R -e "source('code/get_error_rate_table.R')"


$(PROC)/sobs_table.tsv : code/get_sobs_table.R\
							$(UCHIME_SOBS)\
							$(NOCHIM_SOBS)\
							$(NOERROR_SOBS)
	R -e "source('code/get_sobs_table.R')"


$(PROC)/taxonomy_depth_analysis.tsv : $$(filter data/mothur_pool/V%, $$(RDP) $$(GG) $$(SILVA)) code/consolidate_taxonomy.R
	R -e "source('code/consolidate_taxonomy.R')"


$(PROC)/non_random_analysis.tsv : code/get_one_off_table.R\
				$(UNIQUE_ERROR)
	R -e "source('code/get_one_off_table.R')"




$(FIGS)/figure_1.pdf : code/build_figure1.R\
				$(PROC)/mock.error.report
	R -e "source('code/build_figure1.R')"

$(FIGS)/figure_2.pdf : code/build_figure2.R\
				$(PROC)/error_summary.tsv
	R -e "source('code/build_figure2.R')"

$(FIGS)/figure_3.pdf : code/build_figure3.R\
				$(PROC)/taxonomy_depth_analysis.tsv
	R -e "source('code/build_figure3.R')"

$(FIGS)/figure_4.pdf : code/build_figure4.R\
				$(PROC)/non_random_analysis.tsv
	R -e "source('code/build_figure4.R')"

submission/table_1.% : $(TABLE)/build_table1.Rmd\
				$(REFS)/pacbio.oligos\
				$(PROC)/mock.error.report
	R -e 'library(rmarkdown); render("results/tables/build_table1.Rmd", output_file="table_1.pdf")'
	mv results/tables/table_1.pdf submission/
	mv results/tables/table_1.tex submission/

submission/table_2.% : $(TABLE)/build_table2.Rmd\
				$(PROC)/sobs_table.tsv\
				$(PROC)/mock.error.report
	R -e 'library(rmarkdown); render("results/tables/build_table2.Rmd", output_file="table_2.pdf")'
	mv results/tables/table_2.pdf submission/
	mv results/tables/table_2.tex submission/

submission/Schloss_PacBio16S_PeerJ_2016.% : \
						$(PROC)/error_profile.json\
						$(PROC)/error_summary.tsv\
						$(PROC)/taxonomy_depth_analysis.tsv\
						$(PROC)/sobs_table.tsv\
						$(PROC)/non_random_analysis.tsv\
						\
						peerj.csl\
						references.bib\
						Schloss_PacBio16S_PeerJ_2016.Rmd
	R -e 'render("Schloss_PacBio16S_PeerJ_2016.Rmd", clean=FALSE)'
	mv Schloss_PacBio16S_PeerJ_2016.knit.md $@
	rm Schloss_PacBio16S_PeerJ_2016.utf8.md
	mv Schloss_PacBio16S_PeerJ_2016.pdf submission/Schloss_PacBio16S_PeerJ_2016.pdf
	mv Schloss_PacBio16S_PeerJ_2016.tex submission/Schloss_PacBio16S_PeerJ_2016.tex

write.paper :	Schloss_PacBio16S_PeerJ_2016.Rmd\
				submission/table_1.tex\
				submission/table_1.pdf\
				submission/table_2.tex\
				submission/table_2.pdf\
				$(FIGS)/figure_1.pdf\
				$(FIGS)/figure_2.pdf\
				$(FIGS)/figure_3.pdf\
				$(FIGS)/figure_4.pdf\
				submission/Schloss_PacBio16S_PeerJ_2016.md\
				submission/Schloss_PacBio16S_PeerJ_2016.tex\
				submission/Schloss_PacBio16S_PeerJ_2016.pdf

submission/rebuttal_v1.pdf : doc/paper/rebuttal_v1.md
	pandoc --variable mainfont="Helvetica" --variable fontsize=12pt --variable geometry:margin=1in $^ --latex-engine=xelatex -o $@
