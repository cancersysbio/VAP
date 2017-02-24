use strict;
use Data::Dumper;

my $bed = shift;

my @chrs;
for (1..22){
  push(@chrs, $_);
}
push(@chrs, 'X');
push(@chrs, 'Y');
print STDERR Dumper(\@chrs);


open IN, "$bed";
my %bed;
while ( <IN> ){
  chomp;
  my ($chr, $start, $end, $tags, $starts) = split /\t/;
  push (@{$bed{$chr}}, $starts);
}
close IN;

print STDERR "all bed line stored\n";


foreach my $chr (@chrs){
  print "fixedStep chrom=$chr start=1 step=1000 span=1000\n";
  foreach my $starts (@{$bed{$chr}}){
    print "$starts\n";
  }
}
