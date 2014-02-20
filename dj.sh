make_temp(){
	mktemp 2>/dev/null || mktemp -t $0
}

make_retainable_introns(){
## [    ]—-intron--[   ]
## [                   ]
	intersectBed -a $1 -b $2 -wa -wb \
	| awk -v OFS="\t" '{ split($5,a,",");
		if($8 == a[1] && $9 == a[2]){ ## match boundary
			print $1,$2,$3,$4,0,$6;
		}
	}' | sort -u  
	## unique retained intron candidates
}

bed12_to_junction(){
    perl -ne '
    BEGIN { 
        my %junction=(); ## interaction between exons (splicing events)
    }
    chomp;
    my @a=split/\t/,$_;
    my $chr=$a[0]; my $start=$a[1]; my $end=$a[2]; my $score = $a[4]; my $strand =$a[5];
    my $n=$a[9]; my @sizes = split /,/,$a[10]; my @starts = split /,/,$a[11];
    if($n < 2){ next;}
    for(my $i=0;$i < $n-1; $i++){
        my $ss5 = $start+$starts[$i]+$sizes[$i]-1; 
        my $ss3 = $start+$starts[$i+1]+1;
        my $k = $chr."\t".$ss5."\t".$ss3;
        $junction{$k} += $score;
    }
    END { 
        foreach my $k ( keys %junction){
            print $k,"\t",$junction{$k},"\n";
        }
    }
    ' 
}

gtf_to_bed12(){

tmp=`cat << EOF
#!/usr/bin/env python
"""
This source code was obtained from MATs program
(http://intron.healthcare.uiowa.edu/MATS/): 
processGTF.SAMs.py  file
"""

import sys

chunk = 1000
geneGroup = {}
genes = {}
supple = {}
cds={}
for line in sys.stdin: ## for each line
  ele = line.strip().split('\t');
  chr = ele[0];
  type = ele[2]; ## exon, intron, CDS, start_codon, stop_codon..
  sC = ele[3]; ## start coord, 1-base
  eC = ele[4]; ## end coord, 1-base
  group = range(int(sC)/chunk, int(eC)/chunk + 1); ## groups this line could belong to
  group = list(set(group));  ## remove duplicate groups
  strand = ele[6];
  desc = ele[8].split(';');
  gID=['','']; txID=['','']; ## init..
  for dEle in desc: ## for each element of description
    if len(dEle.strip())==0:
      continue; ## probably the last description
    dName = dEle.strip().split(' ')[0];
    dVal = dEle.strip().split(' ')[1];
    if dName.upper() == 'GENE_ID': ## it is a description for gene_id
      gID = [dName,dVal];
    elif dName.upper() == 'TRANSCRIPT_ID': ## it is a description for transcript_id
      txID = [dName, dVal];

  if gID[0].upper()!='GENE_ID' or txID[0].upper() != 'TRANSCRIPT_ID': ## wrong one..
    print("This line does not have correct description for gID or txID: %s, %s" % (gID, txID));
    print("Incorrect description: %s" % ele);
    continue; ## process next line

  for i in group: ## for each possible group
    if i in geneGroup: ## this group already exist
      geneGroup[i].append(gID[1]); ## duplicate geneIDs will get removed after the outer for loop
    else: ## first time accesing this group
      geneGroup[i] = [gID[1]];

  if type=='exon':  ## process exon
    if gID[1] in genes: # already processed this gID
      if txID[1] in genes[gID[1]]: ## this transcript is added already
        genes[gID[1]][txID[1]].append([int(sC), int(eC)]); ## add exon to the existing Tx
      else: ## first time processing this Tx
        genes[gID[1]][txID[1]] = [[int(sC), int(eC)]]; ## add first exon
    else:  ## new gene ID
      genes[gID[1]] = {};
      genes[gID[1]][txID[1]] = [[int(sC), int(eC)]]; ## add first exon
      supple[gID[1]] = [gID[1], chr, strand]; ## geneID, chromosom and strand
  if type=='CDS': ## coding region
    if gID[1] in cds: # already processed this gID
      if txID[1] in cds[gID[1]]: ## this transcript is added already
        cds[gID[1]][txID[1]].append([int(sC), int(eC)]); ## add CDS to the existing Tx
      else: ## first time processing this Tx
        cds[gID[1]][txID[1]] = [[int(sC), int(eC)]]; ## add first CDS
    else:  ## new gene ID
      cds[gID[1]] = {};
      cds[gID[1]][txID[1]] = [[int(sC), int(eC)]]; ## add first exon

## get unique gene lists per group 
for gg in geneGroup: ## for all groups in geneGroup
  geneGroup[gg] = list(set(geneGroup[gg]));

nGene=len(genes); ## number of genes in genes dict
nTx=0; ## number of transcripts
oneTx=0; ## number of one-tx genes
nExon = 0; ## number of exons
oneExon=0; ## number of one-exon transcripts
oneTxOneExon=0;

for id in genes: ## for each gene
  nTx += len(genes[id]); 
  if len(genes[id])==1:
    oneTx += 1; ## one-transcript gene
  for tx in genes[id]: ## for each tx
    nExon += len(genes[id][tx]);
    #print id, tx, genes[id][tx]
    chrom = supple[id][1]
    strand = supple[id][2]

	## 1base to 0base
    ends = map(lambda x: x[1], genes[id][tx])
    starts = map(lambda x: x[0]-1, genes[id][tx]) # 0base
    start = min(starts)
    end = max(ends)
    rstarts0 = map(lambda x: x - start,starts)
    sizes0 = map(lambda x: x[1]-x[0]+1, genes[id][tx]) # 0base 

	## sort by rstarts oannes
    order = zip(*sorted((e,i) for i,e in enumerate(rstarts0)))[1]
    sizes = map(lambda x: sizes0[x], order);
    rstarts = map(lambda x: rstarts0[x], order);

#chr1	12140	12177	HISEQ:69:C2675ACXX:5:1103:8808:4762/1	0	+	12140	12177	255,0,0	1	37	0
    print '\t'.join(map(str, (chrom,start,end,id.replace('"',"")+'::'+tx.replace('"',""),\
		0,strand, start,end,"0,0,0", len(sizes),','.join(map(str,sizes)),','.join(map(str,rstarts)))))
    if len(genes[id][tx])==1: ## one exon tx
      oneExon += 1;
      if len(genes[id])==1: ## one tx gene
        oneTxOneExon+=1;
sys.stderr.write("There are %d distinct gene ID in the gtf file\n" % nGene);
sys.stderr.write("There are %d distinct transcript ID in the gtf file\n" % nTx);
sys.stderr.write("There are %d one-transcript genes in the gtf file\n" % oneTx);
sys.stderr.write("There are %d exons in the gtf file\n" % nExon);
sys.stderr.write("There are %d one-exon transcripts in the gtf file\n" % oneExon);
sys.stderr.write("There are %d one-transcript genes with only one exon in the transcript\n" % oneTxOneExon);
EOF
`
	python -c "$tmp"
}

bed12_to_introns(){
	##  [es   ]s----e[    ee]
	awk -v OFS="\t" '{
		## take introns
		split($11,sizes,",");
		split($12,starts,",");
		for(i=2;i<= $10;i++){
			## intron
			s = $2 + starts[i-1]+sizes[i-1];	
			e = $2 + starts[i];
			ls = $2 + starts[i-1];
			le = ls + sizes[i-1];
			rs = $2 + starts[i];
			re = rs + sizes[i];
			split($4,gene,"::");
			## genebase
			print $1,s,e,gene[1],ls "," re,$6;
		}	
	}' | sort -u  
}
bed12_to_exons(){
	## [s1  ]e1----[s   ]e----[s2   ]e2
	awk -v OFS="\t" '{
		## take introns
		split($11,sizes,",");
		split($12,starts,",");
		for(i=1;i<= $10;i++){
			s = $2 + starts[i];
			e = $2 + starts[i]+sizes[i];	
			s1 = -1; e1 = -1;
			s2 = -1; e2 = -1;
			if( i > 1){ s1 = $2 + starts[i-1]; e1 = s1 + sizes[i-1]; }
			if( i < $10){ s2 = $2 + starts[i+1]; e2 = s2 + sizes[i+1]; }
			split($4,gene,"::");
			## gene base
			print $1,s,e,gene[1],s1 "," e1 "," s2 "," e2,$6;
		}	
	}' | sort -u  
}




each_chrom(){
	## lambda function accepts chrom and size parameters
	## lambda(){ 
	##   # .. handle $1 $2 
	## } 
	chrom_size_file=$1; lambda_func=$2;
	tmp=( `cat $chrom_size_file` )
	for (( i=0; i< ${#tmp[@]}; i+=2 ))
	do
		chrom=${tmp[$i]}; size=${tmp[$i+1]};
		$lambda_func $chrom $size
	done
}

bed12_to_junction(){
	## IN:     [  ]-[  ] : x 
    ##        [   ] [ ]  : y
    ## OUT:       [ ]    : x + y
	perl -ne '
	BEGIN { 
		my %junction=(); ## interaction between exons (splicing events)
	}
	chomp;
	my @a=split/\t/,$_;
	my $chr=$a[0]; my $start=$a[1]; my $end=$a[2]; my $score = $a[4]; my $strand =$a[5];
	my $n=$a[9]; my @sizes = split /,/,$a[10]; my @starts = split /,/,$a[11];
	if($n < 2){ next;}
	for(my $i=0;$i < $n-1; $i++){
		my $ss5 = $start+$starts[$i]+$sizes[$i]-1; 
		my $ss3 = $start+$starts[$i+1]+1;
		my $k = $chr."\t".$ss5."\t".$ss3;
		#$junction{$k} += $score; ## when you use junction.bed directly
		$junction{$k} ++;
	}
	END { 
		foreach my $k ( keys %junction){
			print $k,"\t",$junction{$k},"\n";
		}
	}
	' 
}

bam_to_junctions(){
	## split bam by chrom to reduce memory usage
	bam=$1; chromsize=$2; quality=$3;
	lambda(){
		samtools view -bq $quality $bam $1 | bamToBed -bed12 | awk '$10 > 1' | bed12_to_junction
	}
	each_chrom $chromsize lambda
}
#bam_to_junctions $1 $2 255
count_exon_junction_events(){
    ## input: exon and junction_count
    ## output: exon skipping and junction counts 
    EXON=$1;JUNCTION=$2;
	intersectBed -a $EXON -b $JUNCTION -wa -wb  \
    | awk -v OFS="\t" '{
        ## format:  chr1 249152520   249152523  gene ls,le,rs,re strand chr1  249152519   249153117   12
		## [ls  le]----[   ]---[rs   re]
        s=$8;e=$9;c=$10; ## junction start, end, count
		split($5,a,",");
		l = a[2]-1; r = a[3]+1;  # le -1,  rs+1 
        x=0;y=0;z=0;
        if(s==l && e==$2+1){ ##[  ]^[   ] : left splicing  
            x=c;
        }else if(s==$3-1 && e == r ){ ##   [    ]^[  ]: right splicing 
            y=c;
        }else if(s==l && e == r){ ## [  ]/  [  ] \[  ]skipping
            z=c;
        }
		if( x+y+z > 0){
			print $1,$2,$3,$4,l+1 "," r-1,$6,x,z,y; ## left splicing, skipping, right splicing
		}
    }' | sort -k1,1 -k2,3n -k4,5 | groupBy -g 1,2,3,4,5,6 -c 7,8,9 -o sum,sum,sum
}

count_nonjunction_events(){
	#   --  --  ---  : left unspliced, within, right unspliced
    #   -----------  : within = 0, left unspliced == right unspliced 
	#    [       ]
	bed6=$1; bam=$2; chromsize=$3; quality=$4
	lambda(){
		samtools view -bq $quality $bam $1| bamToBed -bed12 \
		| intersectBed -a $bed6 -b stdin -wa -wb \
		| awk -v OFS="\t" '{ 
			L = 0; C = 0; R = 0;
			if($8 < $2 && $9-1 > $2){ L += 1;}  # left unspliced
			if($8 < $3-1 && $9 > $3){ R += 1;}  # right unspliced 
			if($8 > $2 && $9 < $3){ C += 1;}    # within 
			if( L + C + R > 0){
				print $1,$2,$3,$4,$5,$6,L,C,R;
			}
		}' | sort -k1,1 -k2,3n -k4,6 | groupBy -g 1,2,3,4,5,6 -c 7,8,9 -o sum,sum,sum  
	}
	each_chrom $chromsize lambda
}
#count_nonjunction_events $1 $2 $3 $4

count_intron_junction_events(){
	intron=$1; junction=$2;
	intersectBed -a $intron -b $junction -wa -wb  | awk -v OFS="\t" '{
		##     /       \
		## [   ]-------[    ]
		if($2-1 == $8 && $3+1 == $9){ 
			print $1,$2,$3,$4,$5,$6,$10;
		}
	}'
}
		
count_intron_events(){
	intron=$1;bam=$2;junction=$3;chromsize=$4;quality=$5;
	tmp1=`make_temp`
	tmp2=`make_temp`
	count_intron_junction_events  $intron $junction  > $tmp1
	count_nonjunction_events $intron $bam $chromsize $quality > $tmp2
	intersectBed -a $tmp2 -b $tmp1 -wa -wb -f 1 -r -s | awk -v OFS="\t" '{
		print $1,$2,$3,$4,$5,$6,$7,$8,$9,$16;	
	}'
}
#count_intron_events Events/rintrons.bed ../Tophat/Wt1/accepted_hits.bam Events/Wt1/jc.bed Data/chrom.size 255
