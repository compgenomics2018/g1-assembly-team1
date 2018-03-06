#!/usr/bin/perl


use strict;
use warnings;

# warning: the script was partially hardcoded.
# the hardcoded part is not shown here but in the precusor file
my $pool = $ARGV[0];
my $sourceFile = $ARGV[1];
open (INPUT, $pool)
        or die "Could not open the file: $pool. Exit! $!";

my $tempLine1 = <INPUT>; # skip the first line
my %hash1;
my $counter = 0;
while ($tempLine1 = <INPUT>)
{
        my @tempArr1 = (split '\s+', $tempLine1);
        if ($tempArr1[0] =~ /^SRR/) # if the first column starts with SRR
        {
                $counter += 1;
                # whenever there's a new key, define it with address of arrays
                if (defined $hash1{$tempArr1[2]})
                {
                        $hash1{$tempArr1[2]} = "$hash1{$tempArr1[2]}\n\t$tempArr1[0]";
                }
                else
                {
                        $hash1{$tempArr1[2]} = $tempArr1[0];
                }
        }
}
my $size = scalar keys %hash1;
#print "Test1 $size Counter: $counter\n";
close INPUT;

open (INPUT, $sourceFile)
        or die "Could not open the file: $sourceFile. Exit! $!";

my $tempLine2;
while ($tempLine2 = <INPUT>)
{
        chomp $tempLine2;
        if ( defined $hash1{$tempLine2} )
        {
#       print "$hash1{$tempLine2} \n";
        print "$tempLine2 $hash1{$tempLine2} \n";
        }
        else
        {
                print "missing(need to find mannually): $tempLine2 \n";
        }

}


close INPUT;

