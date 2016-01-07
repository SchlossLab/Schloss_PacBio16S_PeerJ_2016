REGION=$1

START=0
END=0


# | Region | Start | End   | Length | diffs |
# |--------|-------|-------|--------|-------|
#x| v13    | 1046  | 13125 | 489    |    4  | 
#x| v15    | 1046  | 27659 | 879    |    8  | 
#x| v16    | 1046  | 34113 | 1028   |   10  |
#x| v19    | 1046  | 43116 | 1458   |   14  |
#x| v35    | 6428  | 27659 | 545    |    5  |
#x| v4     | 13862 | 23444 | 253    |    2  |
# |--------|-------|-------|--------|-------|

if [ $REGION == V1V9 ]
then
	START=1046
	END=43116
elif [ $REGION == V1V3 ]
then
	START=1046
	END=13125
elif [ $REGION == V1V5 ]
then
    	START=1046
        END=27659
elif [ $REGION == V1V6 ]
then
    	START=1046
        END=34113
elif [ $REGION == V3V5 ]
then
    	START=6428
        END=27659
elif [ $REGION == V4 ]
then
    	START=13862
        END=23444
fi

cp data/references/HMP_MOCK.align data/references/HMP_MOCK.$REGION.align
mothur "#pcr.seqs(fasta=data/references/HMP_MOCK.$REGION.align, start=$START, end=$END, outputdir=data/mothur_pool);
	filter.seqs(vertical=T, trump=.)"

mv data/mothur_pool/HMP_MOCK.$REGION.pcr.filter.fasta data/mothur_pool/HMP_MOCK.$REGION.align
rm data/mothur_pool/HMP_MOCK.$REGION.pcr.align
