# will replace this with whatever we pull down from SRA. this step was done to combie the hd5 files into
# subreads.bam files for the june raw data. probably not the way one would want to run these commands. would
# be better to split them into individual pbs scripts

generate_bam_subread_files : code/run_bax2bam.sh 
	bash code/run_bax2bam.sh


# this assumes that the pbccs package is installed in code/pbcss and uses the default. probably not the way
# one would want to run these commands. would be better to split them into individual pbs scripts

generate_bam_files : code/run_ccs.sh
	bash code/run_ccs.sh

# this should result in: 
# 
# data/raw_june/V1-V3/V1-V3_0001.bam  data/raw_june/V1-V5/V1-V5_0006.bam  data/raw_june/V3-V5/V3-V5_0002.bam
# data/raw_june/V1-V3/V1-V3_0002.bam  data/raw_june/V1-V6/V1-V6_0001.bam  data/raw_june/V3-V5/V3-V5_0003.bam
# data/raw_june/V1-V3/V1-V3_0003.bam  data/raw_june/V1-V6/V1-V6_0002.bam  data/raw_june/V3-V5/V3-V5_0004.bam
# data/raw_june/V1-V3/V1-V3_0004.bam  data/raw_june/V1-V6/V1-V6_0003.bam  data/raw_june/V3-V5/V3-V5_0005.bam
# data/raw_june/V1-V3/V1-V3_0005.bam  data/raw_june/V1-V6/V1-V6_0004.bam  data/raw_june/V3-V5/V3-V5_0006.bam
# data/raw_june/V1-V3/V1-V3_0006.bam  data/raw_june/V1-V6/V1-V6_0005.bam  data/raw_june/V4/V4_0001.bam
# data/raw_june/V1-V5/V1-V5_0001.bam  data/raw_june/V1-V9/V1-V9_0001.bam  data/raw_june/V4/V4_0002.bam
# data/raw_june/V1-V5/V1-V5_0002.bam  data/raw_june/V1-V9/V1-V9_0002.bam  data/raw_june/V4/V4_0003.bam
# data/raw_june/V1-V5/V1-V5_0003.bam  data/raw_june/V1-V9/V1-V9_0003.bam  data/raw_june/V4/V4_0004.bam
# data/raw_june/V1-V5/V1-V5_0004.bam  data/raw_june/V1-V9/V1-V9_0004.bam  data/raw_june/V4/V4_0005.bam
# data/raw_june/V1-V5/V1-V5_0005.bam  data/raw_june/V3-V5/V3-V5_0001.bam
# 
# data/raw_october/V1-V5/V1-V5_0001.bam  data/raw_october/V1-V6/V1-V6_0001.bam
# data/raw_october/V1-V9/V1-V9_0001.bam
#
# Also, each of these *.bam files should have a *.scrap.bam and *.ccs_report.csv file with them.




