#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my ($tbl_file, $prefix, $max_evalue, $outdir);
$prefix = "aps4_ncrna";
$max_evalue = 0.01;
$outdir = ".";

GetOptions(
    "tbl=s"        => \$tbl_file,
    "prefix=s"     => \$prefix,
    "max-evalue=f" => \$max_evalue,
    "outdir=s"     => \$outdir,
) or die "Error in command line arguments\n";

die "Usage: cmscan_to_gff3.pl --tbl <tbl> [--prefix <p>] [--max-evalue <e>] [--outdir <d>]\n" unless $tbl_file;

my $fh;
if ($tbl_file =~ /\.gz$/) {
    open $fh, '-|', 'gzip', '-dc', $tbl_file or die "Cannot decompress $tbl_file: $!\n";
} else {
    open $fh, '<', $tbl_file or die "Cannot open $tbl_file: $!\n";
}

my @records;
while (my $line = <$fh>) {
    chomp $line;
    next if $line =~ /^\s*#/ or $line =~ /^\s*$/;
    my @parts = split /\s+/, $line;
    next if @parts < 17;

    my $target_name = $parts[0];
    my $target_acc  = $parts[1];
    my $query_name  = $parts[2];
    my $seq_from    = $parts[7];
    my $seq_to      = $parts[8];
    my $strand      = $parts[9];
    my $score       = $parts[14];
    my $evalue      = $parts[15];

    next if $evalue > $max_evalue;

    my $start = $seq_from < $seq_to ? $seq_from : $seq_to;
    my $end   = $seq_from > $seq_to ? $seq_from : $seq_to;
    my $desc  = @parts > 17 ? join(" ", @parts[17..$#parts]) : $target_name;

    push @records, {
        chrom  => $query_name,
        start  => int($start),
        end    => int($end),
        strand => $strand,
        score  => $score,
        evalue => $evalue,
        target => $target_name,
        acc    => $target_acc,
        desc   => $desc,
    };
}
close $fh;

@records = sort { $a->{chrom} cmp $b->{chrom} || $a->{start} <=> $b->{start} } @records;

my $gff3_path = "$outdir/$prefix.gff3";
my $bed_path  = "$outdir/${prefix}_exclusion.bed";

open my $gff, '>', $gff3_path or die "Cannot write to $gff3_path: $!\n";
open my $bed, '>', $bed_path  or die "Cannot write to $bed_path: $!\n";

print $gff "##gff-version 3\n";

my $idx = 1;
for my $r (@records) {
    my $feat_id = "ncRNA_${prefix}_$idx";
    print $gff "$r->{chrom}\tInfernal\tncRNA\t$r->{start}\t$r->{end}\t$r->{score}\t$r->{strand}\t.\tID=$feat_id;Name=$r->{target};Accession=$r->{acc};Evalue=$r->{evalue};Note=$r->{desc}\n";
    my $bed_start = $r->{start} - 1;
    my $int_score = int($r->{score});
    print $bed "$r->{chrom}\t$bed_start\t$r->{end}\t$feat_id\t$int_score\t$r->{strand}\n";
    $idx++;
}
close $gff;
close $bed;

printf STDERR "[cmscan_to_gff3] Emitted %d ncRNA features to %s and %s\n", scalar(@records), $gff3_path, $bed_path;
