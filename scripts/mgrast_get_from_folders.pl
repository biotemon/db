#!/usr/bin/perl
#This is version 2021MAR12
use strict;
use warnings;
use 5.010;
use Data::Dumper qw(Dumper);
use List::MoreUtils qw(first_index);

#This script aims to generate a counts file by using the list of folders with results from mgrast in cvs format.

my ($com1, $com2, $res1, $res2, $path);
my (@files1, @files2, $file1);
my ($genuspath, @genus_vec);
my ($i, $j, @counts_vec, @samples_index);
my ($outfile1, $outfile2, $fh);

$path = $ARGV[0];
$com1 = "ls ".$path;
$res1 = `$com1`;
chomp ($res1);
@files1 = split("\n", $res1);

$j = 1;

$outfile1 = 'mgrast_raw_counts.txt';
open($fh, '>', $outfile1) or die "Could not open file '$outfile1' $!";

foreach $file1 (@files1){
    $genuspath = $path.$file1."/genus.csv";
    $com2 = "cat ".$genuspath;
    $res2 = `$com2`;
    chomp ($res2);
    @files2 = split("\n", $res2);
    @genus_vec = split("\t", $files2[0]);
    @counts_vec = split("\t", $files2[1]);
    for ($i=0; $i < scalar(@genus_vec); $i++){
        print $fh $genus_vec[$i]."\t".$j."\t".$counts_vec[$i]."\n";
    }
    push @samples_index, $j;
    $j++;
}
close $fh;

$outfile2 = 'mgrast_index_samples.txt';
open($fh, '>', $outfile2) or die "Could not open file '$outfile2' $!";

for ($i=0; $i < scalar(@samples_index); $i++){
print $fh $files1[$i]."\t".$samples_index[$i]."\n";
}

close $fh;
