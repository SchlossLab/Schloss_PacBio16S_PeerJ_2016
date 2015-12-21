# utility function to print various variables. For example, running the following at the command line:
#	make print-BAM
#
# will generate:
#	BAM=data/raw_june/V1-V3_0001.bam data/raw_june/V1-V3_0002.bam data/raw_june/V1-V3_0003.bam data/raw_june/V1-V3_0004.bam data/raw_june/V1-V3_0005.bam data/raw_june/V1-V3_0006.bam data/raw_june/V1-V5_0001.bam data/raw_june/V1-V5_0002.bam data/raw_june/V1-V5_0003.bam data/raw_june/V1-V5_0004.bam data/raw_june/V1-V5_0005.bam data/raw_june/V1-V5_0006.bam data/raw_june/V1-V6_0001.bam data/raw_june/V1-V6_0002.bam data/raw_june/V1-V6_0003.bam data/raw_june/V1-V6_0004.bam data/raw_june/V1-V6_0005.bam data/raw_june/V1-V9_0001.bam data/raw_june/V1-V9_0002.bam data/raw_june/V1-V9_0003.bam data/raw_june/V1-V9_0004.bam data/raw_june/V3-V5_0001.bam data/raw_june/V3-V5_0002.bam data/raw_june/V3-V5_0003.bam data/raw_june/V3-V5_0004.bam data/raw_june/V3-V5_0005.bam data/raw_june/V3-V5_0006.bam data/raw_june/V4_0001.bam data/raw_june/V4_0002.bam data/raw_june/V4_0003.bam data/raw_june/V4_0004.bam data/raw_june/V4_0005.bam data/raw_october/V1-V5_0001.bam data/raw_october/V1-V6_0001.bam data/raw_october/V1-V9_0001.bam

print-%:
	@echo '$*=$($*)'



###########################################################################################################
#
# Part 1:
#
# The goal of the first part was to take the PacBio subreads files and convert them into fasta and qual
# files. In some cases there were multiple movies per region per data generation event (e.g. june and
# october of 2015). The project root directory has several folders that are useful for this makefile,
# namely data/ and code/. Within data/ are raw_june/, mothur_june, raw_october, raw_june. These correspond 
# to data dumps that we received from pacbio in the specified month. the october data was for the longer
# regions with longer movies.
#
###########################################################################################################


# will replace this with whatever we pull down from SRA. this step was done to combie the hd5 files into
# subreads.bam files for the june raw data. probably not the way one would want to run these commands. would
# be better to split them into individual pbs scripts

generate_bam_subread_files : code/run_bax2bam.sh 
	bash code/run_bax2bam.sh



# this assumes that the pbccs package is installed in code/pbcss and uses the default. probably not the way
# one would want to run these commands. would be better to split them into individual pbs scripts. it also
# assumes that all of the raw data are stored like: data/raw_<month>/<region>_<movie>.subreads.bam

SR_BAM = $(wildcard data/raw*/*.subreads.bam)
BAM = $(subst .subreads.bam,.bam,$(SR_BAM))

.SECONDEXPANSION:
$(BAM) : $$(subst .bam,.subreads.bam,$$@)
	code/smrtanalysis/smrtcmds/bin/ccs $@ $^ --numThreads=8



# each of these *.bam files should have a *.scrap.bam and *.ccs_report.csv file with them. We'll extract the
# fastq files from each ccs-processed bam file using samtools:

BAM_FASTQ = $(subst .subreads.bam,.bam.fastq,$(SR_BAM))

.SECONDEXPANSION:
$(BAM_FASTQ) : $$(subst .fastq,,$$@)
	samtools bam2fq $^ > $@



# Now we need to extract the number of passes and predicted error rate from the bam file:

BAM_CCS_STATS = $(subst .subreads.bam,.bam.ccs_stats,$(SR_BAM))

.SECONDEXPANSION:
$(BAM_CCS_STATS) : $$(subst .ccs_stats,,$$@)
	samtools view $^ | cut -f 1,17,16 | sed -s "s/np:i://" | sed -s "s/rq:f://" > $@



# Finally, we want to pool the *.fastq and *.ccs_stats files for each region within each data drop

data/mothur_%.ccs_stats : $(BAM_CCS_STATS)
	cat $(join $(join data/raw_,$*),_*.ccs_stats) > $@

data/mothur_%.fastq : $(BAM_FASTQ)
	cat $(join $(join data/raw_,$*),_*.fastq) > $@



# Let's extract the fasta and quality score data from the fastq files

data/mothur_%.fasta data/mothur_%.qual : data/mothur_%.fastq
	mothur "#fastq.info(fastq=$^, pacbio=T)"



