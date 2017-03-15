#!/usr/bin/perl

$files = shift;
$outdir = shift;

$outfile = shift;
$minSupp = shift;

$progSetup = shift;
$targetBAM = shift;
$clean = shift;

#$oncotator = shift;
#$oncotator_bin = shift;
#$oncotator_datasource = shift;
#$oncotator_transcript = shift;
#$oncotator_build = shift;


open IN,"<$progSetup" or die "Cant open $file\n";
while(<IN>){
	chomp;
	next if(!$_ || m/^\#/);

	@arr = split /[\t\s\=\:]+/;

	if($arr[0] eq "TUMOR_NORMAL_PAIRS"){
		push @{$SETUP{$arr[0]}},$arr[1];
	}
	else{
		$SETUP{$arr[0]}=$arr[1];
	}
}
close IN;

$oncotator = $SETUP{'ONCOTATOR_ANNOT_LOCAL'};		#/usr/local/bin/oncotator
$oncotator_bin = $SETUP{'ONCOTATOR_BIN'};		#/usr/local/bin/oncotator
$oncotator_datasource = $SETUP{'ONCOTATOR_DATASOURCE'};	#/home/lolab/datasource/oncotator_v1_ds_Jan262015/
$oncotator_transcript = $SETUP{'ONCOTATOR_TRANSCRIPT'};	#/home/lolab/datasource/oncotator_v1_ds_Jan262015/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt
$oncotator_build = $SETUP{'ONCOTATOR_BUILD'};
$combineDinuc = $SETUP{'COMBINE_DINUC'};

$commonDBSNP = $SETUP{'COMMON_DBSNP'};
$clinicalDBSNP = $SETUP{'CLINICAL_DBSNP'};
$exacSites = $SETUP{'EXAC_COUNT_GEQ_TWO'};
$ownGermline = $SETUP{'GERMLINE_OUR_OWN'};

@files = split /\,/,$files;
$total = scalar(@files);
@arr1;
%DETAIL=();
%COUNT=();
%COUNT2=();

foreach $file (@files){
	#print "$file\n";
	#exit;
	$method = (split /\./,(split /\/+/,$file)[-3])[0];

	open IN11,"<$file";
	$mutect=0;
	if($file =~ m/MUTECT/){
		$mutect=1;
		$header = <IN11>;
		$header = <IN11>;
		$ref = 3;
		$mut = 4;
	}
	else{
		$header = <IN11>;
		chomp($header);
		$ref = 2;
		$mut = 3;
	}

	while(<IN11>){
		chomp;
		@arr1 = split /\t/;
		if($method eq "SomaticIndelDetector"){
			@mut = split /[\/\>]+/,$arr1[2];
			$DETAIL{"$arr1[0]\t$arr1[1]\t$mut[0]\t$mut[1]"}{$method}="$arr1[-3]:$arr1[-2]\|$arr1[-1]";	
		}
		else{
			$DETAIL{"$arr1[0]\t$arr1[1]\t$arr1[$ref]\t$arr1[$mut]"}{$method}="$arr1[-2]\|$arr1[-1]";
		}
		$COUNT{$method}++;
	}
	close IN11;
}

if(! -e "$outdir/$outfile"||$clean>0){
	$methods = join "\t",sort(keys %COUNT);
	open OUT11,">$outdir/$outfile";
	print OUT11 "contig\tposition\tref_allele\talt_allele\tCountMethod\t$methods\n";
	foreach $pos (sort keys %DETAIL){
		@methods = sort keys %{$DETAIL{$pos}};
		$COUNT2{(join("\_",@methods))}++;

		next if(($ct=(scalar @methods))/$total<=$minSupp);
		print OUT11 "$pos\t$ct";

		foreach $method (sort keys %COUNT){
			if(exists $DETAIL{$pos}{$method}){
				print OUT11 "\t$DETAIL{$pos}{$method}";
			}
			else{
				print OUT11 "\tNA";
			}
		}
		print OUT11 "\n";		
	}
	close OUT11;
}

$outfile2=$outfile;
$outfile3=$outfile;
$outfile4=$outfile;
$outprefix=$outfile;
$outfile2=~ s/\.xls$/\.all\.oncotator\.xls/;
$outfile3=~ s/\.xls$/\.gene\.oncotator\.xls/;
$outfile4=~ s/\.xls$/\.RNA\.oncotator\.xls/;
$outprefix=~ s/\.xls$//;

@tmp0=(split /\./,$outfile);
$currSample = "$tmp0[0]\.$tmp0[1]"; ##sample name and "indel"/"snv"

print "perl $oncotator $outdir/$outfile $outdir/$outprefix $oncotator_bin $oncotator_datasource $oncotator_transcript $oncotator_build $commonDBSNP $clinicalDBSNP $exacSites $ownGermline $clean\n";
#exit;
if(!-e "$outdir/$outfile2"||$clean>0){
	system "perl $oncotator $outdir/$outfile $outdir/$outprefix $oncotator_bin $oncotator_datasource $oncotator_transcript $oncotator_build $commonDBSNP $clinicalDBSNP $exacSites $ownGermline $clean";

	system "mkdir -p $outdir/OncotatorTmp/";
	system "mv $outdir/$currSample\.\*\.tmp $outdir/OncotatorTmp/";
	system "mv $outdir/$currSample\.\*\.maflite $outdir/OncotatorTmp/";
}
#exit;



exit if($outfile3 !~ m/snv/);

## get further annotation for the dinuc and trinuc mutation in gene annotation
$outfileConsec=$outfile;
$outfileConsec=~ s/\.xls$/\.gene\.consecMut\.xls/;

if(-e "$outdir/$outfile2"  && (! -e "$outdir/$outfileConsec"||$clean>0)){
	open IN11,"<$outdir/$outfile2";
	$_=<IN11>;
	$prevLine = <IN11>;
	chomp($prevLine);
	@arrPrev = split /\t/,$prevLine;
	%LINES=();
	%CHR=();
	%START=();
	%END=();
	%REF=();
	%MUT=();

	#print "@arrPrev\n";
	#exit;

	while(<IN11>){
		chomp;
		@arr = split /\t/;
		
		if($arr[4] eq $arrPrev[4] && $arr[8] ne "NA" && $arrPrev[8] ne "NA"){
			$AA1 = $arr[8];
			$AA1 =~ s/[A-Z\*]+$//g;
			$AA0 = $arrPrev[8];
			$AA0 =~ s/[A-Z\*]+$//g;

			if($AA1 eq $AA0){
				#print "@arrPrev\n@arr\n";
				#exit;

				@tmp1 = split /\>/,$arr[2];
				if(! exists $REF{"$arr[4]\t$AA1"}){
					@tmp0 = split /\>/,$arrPrev[2];
					$REF{"$arr[4]\t$AA1"}=$tmp0[0];					
					$MUT{"$arr[4]\t$AA1"}=$tmp0[1];

					push @{$LINES{"$arr[4]\t$AA1"}},(join "\t",@arrPrev[13..16]);
					$START{"$arr[4]\t$AA1"}=$arrPrev[1];
					$CHR{"$arr[4]\t$AA1"}=$arrPrev[0];
				}

				$REF{"$arr[4]\t$AA1"}.=$tmp1[0];
				$MUT{"$arr[4]\t$AA1"}.=$tmp1[1];
				push @{$LINES{"$arr[4]\t$AA1"}},(join "\t",@arr[13..16]);
				$END{"$arr[4]\t$AA1"}=$arr[1];
			}
		}	

		@arrPrev = @arr;
	}
	close IN11;


	open OUT11,">$outdir/$outfileConsec";
	print OUT11 "contig\tposition\tref_allele\talt_allele\tCountMethod\t$methods\n";

	foreach $gene (sort keys %REF){
		@tmp = split /\t/,$gene;		
		$maxLen = length($REF{$gene});

		$lines = `samtools view $targetBAM $CHR{$gene}\:$START{$gene}\-$END{$gene}`;
		chomp($lines);
		
		@lines = split /\n/,$lines;
		$refCount=0;
		$mutCount=0;
		foreach $line (@lines){
			@arr2 = split /\t/,$line;
			$sign1 = $sign2 = $signature = &extractPos($START{$gene},$maxLen,$arr2[3],$arr2[9],$arr2[5]);
			next if(!$signature);			
			
			if($sign1 eq $REF{$gene}){
				$refCount++;
			}
			elsif($sign2 eq $MUT{$gene}){
				$mutCount++;
			}

			#print "$line\n$signature\n\n";
			#exit;
		}
		
		#print "$refCount\t$mutCount\n";
		
		#if($mutCount > 2 && $mutCount/($mutCount+$refCount)>=0.05){
		#if($mutCount > 3 && $mutCount/($mutCount+$refCount)>=0.05){
		if($mutCount > 3){
			$line = ${$LINES{$gene}}[0];
			@line = split /\t/,$line;

			for($i=1;$i<=$#line;$i++){
				$line[$i] =~ s/\(\d+\,\d+\)/\($refCount\,$mutCount\)/;
				$line[$i] =~ s/\:[ACGT]{1}\/[ACGT]{1}\:/\:$REF{$gene}\/$MUT{$gene}\:/;
				$line[$i] =~ s/\:[ACGT]{1}\/[ACGT]{1}\:/\:$REF{$gene}\/$REF{$gene}\:/;
			}
						
			$line = join "\t",@line;
			#print "$line\n";
			#exit;

			#print "$CHR{$gene}\t$START{$gene}\t$END{$gene}\t$REF{$gene}\t$MUT{$gene}\t$line\n";
			print OUT11 "$CHR{$gene}\t$START{$gene}\t$REF{$gene}\t$MUT{$gene}\t$line\n";
		}
	}
	close OUT11;
}

$outfileConsec2=$outfileConsec;
$outfileConsec3=$outfileConsec;
$outfileConsec4=$outfileConsec;
$outprefixConsec=$outfileConsec;
$outfileConsec2=~ s/\.xls$/\.all\.oncotator\.xls/;
$outfileConsec3=~ s/\.xls$/\.gene\.oncotator\.xls/;
$outfileConsec4=~ s/\.xls$/\.RNA\.oncotator\.xls/;
$outprefixConsec=~ s/\.xls$//;

print "perl $oncotator $outdir/$outfile $outdir/$outprefix $oncotator_bin $oncotator_datasource $oncotator_transcript $oncotator_build $commonDBSNP $clinicalDBSNP $exacSites $ownGermline $clean\n";
if(!-e "$outdir/$outfileConsec2"||$clean>0){
	
	system "perl $oncotator $outdir/$outfileConsec $outdir/$outprefixConsec $oncotator_bin $oncotator_datasource $oncotator_transcript $oncotator_build $commonDBSNP $clinicalDBSNP $exacSites $ownGermline $clean";
	system "perl $combineDinuc $outdir/$outfile3 $outdir/$outfileConsec3 $clean";

	system "rm $outdir/$outfileConsec2";
	system "rm $outdir/$outfileConsec4";

	system "mkdir -p $outdir/OncotatorTmp/";
	system "mv $outdir/$currSample\.\*\.tmp $outdir/OncotatorTmp/";
	system "mv $outdir/$currSample\.\*\.maflite $outdir/OncotatorTmp/";
	system "mv $outdir/$currSample\.\*\.noDinuc $outdir/OncotatorTmp/";
	system "mv $outdir/$currSample\*\.consecMut\.\* $outdir/OncotatorTmp/";
}



sub extractPos{
	local $pos = shift;
	local $maxLen = shift;
	local $alignSt = shift;
	local $seq = shift;
	local $cigarOp = shift;
	local $nucl = shift;
	
	local @op = ();
	local $i;
	local $j=0;
	local $k;

	@op = $cigarOp =~ m/(\d+)([MIDNSHPX\=])/g;
	local @status=();
	local $len = length($seq);
	for($i=0;$i<=$#op;$i+=2){
		#next if($op[$i+1] eq "D");
		for($k=0;$j<$len && $k<$op[$i];$j++,$k++){
			$status[$j] = $op[$i+1];
		}
	}
	#print "@status\n";
	#exit;	
	$localpos = $pos - $alignSt;
	@seq = split //,$seq;

	$j=0;
	$k=0;
	if($status[0] eq 'S'){
		while($status[$j] eq 'S'){
			$j++;
			$k++;
		}
	}

	#exit;	
	for($i=0;$i<$localpos && $j<=$#status && $k <=$#seq;$j++){
		if($status[$j] eq 'M'){
			$i++;
			$k++;
		}
		elsif($status[$j] eq 'D'){
			$i++;
		}
		elsif($status[$j] eq 'I'||$status[$j] eq 'S'){
			$k++;
		}
	}
	if($i==$localpos){
		$localpos=$k;
	}
	else{
		$localpos=-1;
	}
	
	#return "" if($localpos<0 || $status[$j] eq 'D' || $seq[$localpos] ne $nucl);
	#return "" if($localpos<0 || $status[$j] eq 'D');	
	return "" if($localpos<0||$localpos+$maxLen-1>$#seq);	

	@sign = ();
	$offset=0;
	$insert=0;
	for($r=$localpos;$r<=$#seq && $r<$localpos+$maxLen && $offset<$maxLen;$r++,$offset++){
		if($status[$j+$offset] eq "M" && $seq[$r] ne 'N'){
			push @sign,$seq[$r];
		}
		if($status[$j+$offset] eq "I" && $seq[$r] ne 'N'){
			push @sign,$seq[$r];
			$insert=1;
		}
		elsif($status[$j+$offset] eq "D"){
			push @sign,"\*";
			$r--;
		}		
	}
	
	#print "@sign\n";
	#exit;
	$sign = (join "",@sign);
	if($insert){
		$sign = lc($sign);
	}
	return $sign;
}

