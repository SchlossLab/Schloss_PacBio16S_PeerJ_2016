# utility function to print various variables. For example, running the following at the command line:
#	make print-BAM
#
# will generate:
#	BAM=data/raw_june/V1-V3_0001.bam data/raw_june/V1-V3_0002.bam data/raw_june/V1-V3_0003.bam data/raw_june/V1-V3_0004.bam data/raw_june/V1-V3_0005.bam data/raw_june/V1-V3_0006.bam data/raw_june/V1-V5_0001.bam data/raw_june/V1-V5_0002.bam data/raw_june/V1-V5_0003.bam data/raw_june/V1-V5_0004.bam data/raw_june/V1-V5_0005.bam data/raw_june/V1-V5_0006.bam data/raw_june/V1-V6_0001.bam data/raw_june/V1-V6_0002.bam data/raw_june/V1-V6_0003.bam data/raw_june/V1-V6_0004.bam data/raw_june/V1-V6_0005.bam data/raw_june/V1-V9_0001.bam data/raw_june/V1-V9_0002.bam data/raw_june/V1-V9_0003.bam data/raw_june/V1-V9_0004.bam data/raw_june/V3-V5_0001.bam data/raw_june/V3-V5_0002.bam data/raw_june/V3-V5_0003.bam data/raw_june/V3-V5_0004.bam data/raw_june/V3-V5_0005.bam data/raw_june/V3-V5_0006.bam data/raw_june/V4_0001.bam data/raw_june/V4_0002.bam data/raw_june/V4_0003.bam data/raw_june/V4_0004.bam data/raw_june/V4_0005.bam data/raw_october/V1-V5_0001.bam data/raw_october/V1-V6_0001.bam data/raw_october/V1-V9_0001.bam

print-%:
	@echo '$*=$($*)'



###########################################################################################################
#
# Part 1: Get the references
#
# We will need several reference files to complete the analysis: the mock community sequences, the SILVA
# reference alignment, and the SILVA, RDP, and greengenes reference taxonomies. 
#
###########################################################################################################

REFS = data/references


# We want the latest greatest reference alignment and the SILVA reference alignment is the best reference
# alignment on the market. This version is from v123 and described at
# http://blog.mothur.org/2014/08/08/SILVA-v119-reference-files/.  We will use the full-length version of
# the database, which contains 137,879 bacterial sequences. This also contains the reference taxonomy. We
# will limit the databases to only include bacterial sequences.

$(REFS)/silva.bacteria.% :
	wget -N http://mothur.org/w/images/b/be/Silva.nr_v123.tgz
	tar xvzf Silva.nr_v123.tgz silva.nr_v123.align silva.nr_v123.tax
	mothur "#get.lineage(fasta=silva.nr_v123.align, taxonomy=silva.nr_v123.tax, taxon=Bacteria);degap.seqs(fasta=silva.nr_v123.pick.align, processors=8)"
	mv silva.nr_v123.pick.align $(REFS)/silva.bacteria.align
	mv silva.nr_v123.pick.tax $(REFS)/silva.bacteria.tax
	mv silva.nr_v123.pick.ng.fasta $(REFS)/silva.bacteria.fasta
	rm Silva.nr_v123.tgz mothur*logfile silva.nr_v123.*


# We also want the greengenes reference taxonomy. This version is from the  greengenes v13_8_99 and is
# described at http://blog.mothur.org/2014/08/12/greengenes-v13_8_99-reference-files/

$(REFS)/gg_13_8_99.% :
	wget -N http://www.mothur.org/w/images/6/68/Gg_13_8_99.taxonomy.tgz
	tar xvzf Gg_13_8_99.taxonomy.tgz gg_13_8_99.fasta gg_13_8_99.gg.tax
	mv gg_13_8_99.* $(REFS)/
	rm Gg_13_8_99.taxonomy.tgz


# Next, we want the RDP reference taxonomy. The current version is v10 and we use a "special" pds version
# of the database files, which are described at http://blog.mothur.org/2014/10/28/RDP-v10-reference-files/

$(REFS)/trainset14_032015.% :
	wget -N http://www.mothur.org/w/images/8/88/Trainset14_032015.pds.tgz
	tar xvzf Trainset14_032015.pds.tgz trainset14_032015.pds/trainset14_032015.pds.*
	mv trainset14_032015.pds/* $(REFS)/
	rmdir trainset14_032015.pds
	rm Trainset14_032015.pds.tgz


# Finally, we want to align the mock community reference sequences to our newly created
# silva.bacteria.fasta file...

$(REFS)/HMP_MOCK.% :
	wget --no-check-certificate -N -P $(REFS) https://raw.githubusercontent.com/SchlossLab/Kozich_MiSeqSOP_AEM_2013/master/data/references/HMP_MOCK.fasta
	mothur "#align.seqs(fasta=$(REFS)/HMP_MOCK.fasta, reference=$(REFS)/silva.bacteria.align)"




###########################################################################################################
#
# Part 2: BAM to fasta/qual files
#
# The goal of this part was to take the PacBio subreads files and convert them into fasta and qual files.
# In some cases there were multiple movies per region per data generation event (e.g. june and  october of
# 2015). The project root directory has several folders that are useful for this makefile, namely data/
# and code/. Within data/ are raw_june/, mothur_june, raw_october, raw_june. These correspond to data dumps
# that we received from pacbio in the specified month. the october data was for the longer regions with
# longer movies.
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



