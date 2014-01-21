#!/usr/bin/perl -w
use strict;
my $usage = " usage: $0 -i <bed> -a <genome directory>";
if(scalar @ARGV ne 4){ print $usage,"\n";exit(-1);}

my $fin;
my $fh;
my $dir;
my %genome = ();

while(@ARGV){
	my $e = shift @ARGV;
	if($e eq "-i"){ $fin = shift @ARGV;}
	elsif($e eq "-a"){ $dir =shift @ARGV;}
	else{ print $usage,"\n";exit(-1);}
}
if($fin eq "stdin"){ $fh = *STDIN;
}else{ open($fh, $fin) or die "$fin $!";}

if(substr($dir,-1,1) ne "\/"){ $dir.= "\/";}
if(-d $dir){
	opendir(D,$dir) or die "$!";
      	my @ff = readdir(D);
	foreach my $f (@ff){
		if($f=~ /(\w+).fa/){
			$genome{$1} = $dir.$f;
		}
	}

	close(D);
}else{
	print $usage,"\nno such dir $dir\n";exit(-1);
}
#foreach my $chr (keys %genome){ print $chr," ",$genome{$chr},"\n"; }


my %data = ();
while(<$fh>){
	chomp;$_=~s/\r//g;
	next if(substr($_,0,1) eq "#");
	my ($chr) = split /\s/,$_;
	push @{$data{$chr}},$_;
}
if($fin ne "stdin"){ close($fh);}

foreach my $chr (keys %data){
	my $seq = "";
	my $f = $genome{$chr};
	next unless defined $f;
	open (F, $f) or die "$! $f";
	while(<F>){
		chomp;
		next if (substr($_,0,1) eq ">");
		$seq .= $_;
	}
	close(F);
	foreach my $e (@{$data{$chr}}){
		#my ($chr, $start, $end, $name, $score, $strand,@other) = split /\s/,$e;
		my ($chr,$start,$end,$name,$score,$strand,$thickStart,$thickEnd,$itemRgb,$blockCount,$blockSizes,$blockStarts) = split /\s/, $e;
		my $subseq = "";
		if(defined $blockSizes && $blockCount> 0){
			my @ss = split /,/,$blockStarts;
			my @ll = split /,/,$blockSizes;
			for(my $i=0;$i<$blockCount;$i++){
				$subseq .= substr($seq,$start+$ss[$i]+1,$ll[$i]);
			}
		}else{
			$subseq = substr($seq,$start,$end-$start);	
		}	
		my $head = $e; $head =~ s/\t/:/g;
		print ">",$head,"\n";#$chr,$start,$end,$name,$score,$strand\n";
		#print substr($seq,$start,$end-$start),"\n";
		print $subseq,"\n";
	}
}

