#!/bin/bash 

usage="
Author	: Hyunmin Kim (hyun.kim@ucdenver.edu)

Usage	: cat <bed> | $0 <cmd>

<cmd>	:
 	-3 	: take 3\'
 	-5 	: take 5\'
	-w	: switch strand 
	-l <int>	: left flank
	-r <int> 	: right flank
	-s  : strand specific (left and right became 5' and 3' in the strand)
	Note: options will be applied in the order
"
if [ -t 0 ] || [ $# -lt 1 ] || [ $1 == "-h" ]; then 
	echo "$usage"
	exit 0
fi
p3='if($a[5] eq "+"){ $a[1]=$a[2]-1} $a[2]=$a[1]+1;'
p5='if($a[5] eq "-"){ $a[1]=$a[2]-1} $a[2]=$a[1]+1;'
sw='if($a[5] eq "+"){ $a[5]="-"}else{ $a[5]="+";}'
s='$strand=1;'
ll='if($strand && $a[5] eq "-"){$a[2] += $X; }else{$a[1] -= $X; }'
rr='if($strand && $a[5] eq "-"){$a[1] -= $X; }else{$a[2] += $X; }'

#tmp=$( mktemp -t $0 )
#arr=$(echo $cmds | tr "," "\n")
tmp='my $strand=0;'
strand=0
while getopts "sw35l:r:" arg; do
	case $arg in 
		w) tmp="$tmp $sw";;
		3) tmp="$tmp $p3";;
		5) tmp="$tmp $p5";;
		s) tmp="$tmp $s";;
		l) tmp="$tmp ;my \$X=${OPTARG};$ll";; 
		r) tmp="$tmp ;my \$X=${OPTARG};$rr";; 
	esac
done
shift $(( OPTIND - 1 ));

# build script
tmp='chomp;my @a=split/\t/,$_;'$tmp'print join( "\t",@a),"\n";'
# run
#echo $tmp
perl -ne "$tmp"

