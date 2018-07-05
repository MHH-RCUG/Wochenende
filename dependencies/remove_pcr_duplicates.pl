#!perl -w

use warnings;
use strict;

my $input=$ARGV[0];
my $output=$ARGV[1];

my %reads;
my $original_reads = 0;
my $removed_reads = 0;
my $written_reads = 0;

open(INFILE, $input) or die "Can't open ${input}\n";
open(OUTFILE, ">".$output) or die "Can't open ${output}\n";

while(<INFILE>) {
    my $header_a = $_;
    my $read = <INFILE>;
    my $header_b = <INFILE>;
    my $quals = <INFILE>;

    $original_reads++;

    my $r = substr($read, 0, 100);
    #print $read;
    #print $r, "\n";
    
    if (defined $reads{$r}) {
        $removed_reads++;
    } else {
        $reads{$r} = 1;
        print OUTFILE $header_a;
        print OUTFILE $read;
        print OUTFILE $header_b;
        print OUTFILE $quals;
        $written_reads++;
    }   
}

close(OUTFILE);
close(INFILE);

my $pc = 0;

if ($removed_reads > 0) {
    $pc = (100 * $removed_reads) / $original_reads;
}

open(STATFILE, ">".$output.".stats") or die "Can't open stat file\n";
printf STATFILE "%s\t%d\t%d\t%d\t%.2f\n", $input, $original_reads, $removed_reads, $written_reads, $pc;
close(STATFILE);
