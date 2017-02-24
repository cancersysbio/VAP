package bwaMapping;

use strict;

#
# BWA Mapping
#

sub bwaPairMapping {

  my ($class, $bwaBin, $samtoolsBin, $ReadGroup, $threads, $bwaindex, $outBam, $readfiles1, $readfiles2) = @_;

  my $cmd = "$bwaBin mem -r 1.2 -t $threads -R \'$ReadGroup\' $bwaindex \'\<zcat $readfiles1\' \'\<zcat $readfiles2\' | $samtoolsBin view -bS - >$outBam";

  return $cmd;

}


sub bwaSmartMapping {   #smart paring, for interleaved fastqs

  my ($class, $bwaBin, $samtoolsBin, $ReadGroup, $threads, $bwaindex, $outBam, $readfiles1) = @_;

  my $cmd = "$bwaBin mem -p -r 1.2 -t $threads -R \'$ReadGroup\' $bwaindex \'\<zcat $readfiles1\' | $samtoolsBin view -bS - >$outBam";

  return $cmd;

}


sub bwaSingleMapping {

  my ($class, $bwaBin, $samtoolsBin, $ReadGroup, $threads, $bwaindex, $outBam, $readfiles1) = @_;

  my $cmd = "$bwaBin mem -r 1.2 -t $threads -R \'$ReadGroup\' $bwaindex \'\<zcat $readfiles1\' | $samtoolsBin view -bS - >$outBam";

  return $cmd;

}


sub bowtieMappingSnv {

  my ($class, $bowtieBin, $bowtie2index, $inputFa, $outSam, $threads) = @_;

  my $cmd = "$bowtieBin -k 22 -p $threads -f --no-unal -D 15 -R 2 -N 1 -L 20 -i S,1,0.75 --score-min L,-2,-0.3 -x $bowtie2index $inputFa >$outSam";

  return $cmd;
}


sub bamSort {

  my ($class, $samtoolsBin, $threads, $bamTmp, $outBam, $inBam) = @_;

  my $cmd = "$samtoolsBin sort -@ $threads -O bam -T $bamTmp -o $outBam $inBam";

  return $cmd;

}


sub bamIndex {

  my ($class, $samtoolsBin, $inBam) = @_;

  my $cmd = "$samtoolsBin index $inBam";

  return $cmd;
}


sub samToBam {

  my ($class, $samtoolsBin, $inSam, $outBam) = @_;

  my $cmd = "$samtoolsBin view -Sb $inSam -o $outBam";

  return $cmd;
}


sub recalMD {

  my ($class, $samtoolsBin, $inBam, $GFASTA, $outBam) = @_;

  my $cmd = "$samtoolsBin calmd -brE $inBam $GFASTA >$outBam";

  return $cmd;

}


sub indelRealignment1 {

  my ($class, $gatkBin, $inBam, $gfasta, $knownindel1, $knownindel2, $CHR, $outList, $threads, $mem) = @_;

  my $cmd = "java -Xmx$mem -jar $gatkBin -T RealignerTargetCreator --allow_potentially_misencoded_quality_scores -R $gfasta -I $inBam -known $knownindel1 -known $knownindel2 -L $CHR -nt $threads -o $outList";
  if ($CHR eq 'ALL') {
    $cmd = "java -Xmx$mem -jar $gatkBin -T RealignerTargetCreator --allow_potentially_misencoded_quality_scores -R $gfasta -I $inBam -known $knownindel1 -known $knownindel2 -nt $threads -o $outList";
  }

  return $cmd;

}


sub indelRealignment2 {

  my ($class, $gatkBin, $inBam, $gfasta, $targetList, $knownindel1, $knownindel2, $CHR, $outBam, $threads, $mem) = @_;

  my $cmd = "java -Xmx$mem -jar $gatkBin -T IndelRealigner --allow_potentially_misencoded_quality_scores -R $gfasta -I $inBam -targetIntervals $targetList -known $knownindel1 -known $knownindel2 -compress 5 -L $CHR -o $outBam";
  if ($CHR eq 'ALL') {
    $cmd = "java -Xmx$mem -jar $gatkBin -T IndelRealigner --allow_potentially_misencoded_quality_scores -R $gfasta -I $inBam -targetIntervals $targetList -known $knownindel1 -known $knownindel2 -compress 5 -o $outBam";
  }

  return $cmd;

}

sub MarkDuplicates {

  my ($class, $MarkDuplicatesBin, $inBam, $outBam, $metric, $mem) = @_;

  my $cmd = "java -Xmx$mem -jar $MarkDuplicatesBin I=$inBam O=$outBam M=$metric REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true";

  return $cmd;

}


sub BaseRecalibration {

  my ($class, $gatkBin, $inBam, $gfasta, $DBSNP, $knownindel1, $knownindel2, $outTable, $threads, $mem) = @_;

  #--fix_misencoded_quality_scores

  my $cmd = "java -Xmx$mem -jar $gatkBin -T BaseRecalibrator -R $gfasta -I $inBam -knownSites $DBSNP -knownSites $knownindel1 -knownSites $knownindel2 -nct $threads -o $outTable";

  return $cmd;

}

sub BaseRecalibrationPrint {

  my ($class, $gatkBin, $inBam, $gfasta, $inTable, $outBam, $threads, $mem) = @_;

  my $cmd = "java -Xmx$mem -jar $gatkBin -T PrintReads -R $gfasta -I $inBam -BQSR $inTable -DIQ --emit_original_quals -nct $threads -o $outBam";

  return $cmd;

}

1;


