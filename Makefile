print-%:
	@echo '$*=$($*)'

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



