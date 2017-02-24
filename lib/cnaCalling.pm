package cnaCalling;

use strict;

#
# cna Calling
#

sub runTitan {

  my ($class, $RscriptBin, $titanRBin, $PATH, $sampleName, $alleleCount, $tumorWig, $normalWig, $gcWig, $mapWig, $plp, $plpe, $nc, $ncm, $sym, $exons) = @_;

  my $cmd = "$RscriptBin $titanRBin $PATH $sampleName $alleleCount $tumorWig $normalWig $gcWig $mapWig $plp $plpe $nc $ncm $sym $exons";

  return $cmd;

}


1;


