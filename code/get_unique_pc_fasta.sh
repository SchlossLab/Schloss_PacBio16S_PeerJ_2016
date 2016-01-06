FASTA=$1

STUB=$(echo $FASTA | sed -e 's/.fasta//')
echo $STUB

REGION=$(echo $FASTA | sed -re 's/.*\/(V[^.]*)\..*/\1/')
echo $REGION

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
	DIFFS=14
elif [ $REGION == V1V3 ]
then
	START=1046
	END=13125
	DIFFS=4
elif [ $REGION == V1V5 ]
then
    	START=1046
        END=27659
	DIFFS=8
elif [ $REGION == V1V6 ]
then
    	START=1046
        END=34113
	DIFFS=10
elif [ $REGION == V3V5 ]
then
    	START=6428
        END=27659
	DIFFS=5
elif [ $REGION == V4 ]
then
    	START=13862
        END=23444
	DIFFS=2
fi

mothur "#unique.seqs(fasta=$FASTA);
	align.seqs(fasta=current, reference=data/references/silva.seed.align, processors=8);
	screen.seqs(fasta=current, name=current, start=$START, end=$END, maxhomop=8);
	filter.seqs(vertical=T, trump=.);
	unique.seqs(fasta=current, name=current);
	pre.cluster(fasta=current, name=current, diffs=$DIFFS)"

#$STUB.unique.good.filter.names
#$STUB.unique.good.filter.unique.fasta
#$STUB.unique.good.filter.unique.precluster.names
#$STUB.unique.good.filter.unique.precluster.fasta

rm $STUB.unique.good.filter.unique.precluster.map
rm $STUB.unique.good.filter.fasta
rm $STUB.unique.bad.accnos
rm $STUB.good.names
rm $STUB.unique.good.align
rm $STUB.unique.align.report
rm $STUB.unique.align
rm $STUB.unique.flip.accnos
rm $STUB.names
rm $STUB.unique.fasta
