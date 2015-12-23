# utility function to print various variables. For example, running the
# following at the command line:
#
#	make print-BAM
#
# will generate:
#	BAM=data/raw_june/V1V3_0001.bam data/raw_june/V1V3_0002.bam data/raw_june/V1V3_0003.bam data/raw_june/V1V3_0004.bam data/raw_june/V1V3_0005.bam data/raw_june/V1V3_0006.bam data/raw_june/V1V5_0001.bam data/raw_june/V1V5_0002.bam data/raw_june/V1V5_0003.bam data/raw_june/V1V5_0004.bam data/raw_june/V1V5_0005.bam data/raw_june/V1V5_0006.bam data/raw_june/V1V6_0001.bam data/raw_june/V1V6_0002.bam data/raw_june/V1V6_0003.bam data/raw_june/V1V6_0004.bam data/raw_june/V1V6_0005.bam data/raw_june/V1V9_0001.bam data/raw_june/V1V9_0002.bam data/raw_june/V1V9_0003.bam data/raw_june/V1V9_0004.bam data/raw_june/V3V5_0001.bam data/raw_june/V3V5_0002.bam data/raw_june/V3V5_0003.bam data/raw_june/V3V5_0004.bam data/raw_june/V3V5_0005.bam data/raw_june/V3V5_0006.bam data/raw_june/V4_0001.bam data/raw_june/V4_0002.bam data/raw_june/V4_0003.bam data/raw_june/V4_0004.bam data/raw_june/V4_0005.bam data/raw_october/V1V5_0001.bam data/raw_october/V1V6_0001.bam data/raw_october/V1V9_0001.bam

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

REFS = data/references


# We want the latest greatest reference alignment and the SILVA reference
# alignment is the best reference alignment on the market. This version is from
# v123 and described at http://blog.mothur.org/2014/08/08/SILVA-v119-reference-files/.
# We will use the full-length version of the database, which contains 137,879
# bacterial sequences. This also contains the reference taxonomy. We will limit
# the databases to only include bacterial sequences.

$(REFS)/silva.bacteria.% :
	wget -N http://mothur.org/w/images/b/be/Silva.nr_v123.tgz
	tar xvzf Silva.nr_v123.tgz silva.nr_v123.align silva.nr_v123.tax
	mothur "#get.lineage(fasta=silva.nr_v123.align, taxonomy=silva.nr_v123.tax, taxon=Bacteria);degap.seqs(fasta=silva.nr_v123.pick.align, processors=8)"
	mv silva.nr_v123.pick.align $(REFS)/silva.bacteria.align
	mv silva.nr_v123.pick.tax $(REFS)/silva.bacteria.tax
	mv silva.nr_v123.pick.ng.fasta $(REFS)/silva.bacteria.fasta
	rm Silva.nr_v123.tgz mothur*logfile silva.nr_v123.*


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
# created silva.bacteria.fasta file...

$(REFS)/HMP_MOCK.% :
	wget --no-check-certificate -N -P $(REFS) https://raw.githubusercontent.com/SchlossLab/Kozich_MiSeqSOP_AEM_2013/master/data/references/HMP_MOCK.fasta
	mothur "#align.seqs(fasta=$(REFS)/HMP_MOCK.fasta, reference=$(REFS)/silva.bacteria.align)"



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

$(SAMPLE_FASTA) $(SAMPLE_QUAL) $(SAMPLE_GROUPS) : $$(addsuffix .fasta,$$(basename $$(basename $$@))) $$(addsuffix .qual,$$(basename $$(basename $$@))) data/references/pacbio.oligos
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
ERROR_MATRIX = $(subst summary,error.matrix,$(ALIGN_SUMMARY))
ERROR_QUALITY = $(subst summary,error.quality,$(ALIGN_SUMMARY))

$(ALIGN_SUMMARY) : $$(subst filter.summary,fasta,$$@) $$(subst filter.summary,qual,$$@) $$(REFS)/HMP_MOCK.align
	$(eval FASTA = $(word 1, $^))
	$(eval QUAL = $(word 2, $^))
	$(eval MOCK = $(word 3, $^))
	$(eval STUB = $(subst .fasta,,$(FASTA)))
	@echo $(STUB)
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

$(ERROR_SUMMARY) $(ERROR_MATRIX) $(ERROR_QUALITY) : $$(subst error,summary,$$(basename $$@))


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
