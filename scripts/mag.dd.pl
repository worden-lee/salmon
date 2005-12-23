#!/usr/bin/perl
# this script transforms a series of complex numbers to just magnitude
#  i.e input must be 2 columns
my $infile=$ARGV[0];
my $outfile; ($outfile=$infile) =~ s/\.out/.mag.out/;
print "from $infile to $outfile\n";
open IN, $infile || die "couldn't open $infile";
open OUT, ">$outfile" || die "couldn't open $outfile";
my @y0;
my $mean;
while (<IN>)
{
#  next if /^\#/;
  if (/^\(([^,]+),([^,]+)\)/ || /^(\S*) (\S*)/)
  {
    my $mag = sqrt($1*$1 + $2*$2);
    #my $phase = atan2($1,$2);
    print OUT "$mag\n";
  }
}
