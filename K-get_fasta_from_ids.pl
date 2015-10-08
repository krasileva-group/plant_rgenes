#!/usr/bin/perl
# Author        : Ksenia
# Date          : Wed Aug  8 20:42:38 UTC 2007
# Description   : Given an identifier, find a sequence in the database file


# Feed a file containg ids and a database to the script

use strict;
use warnings;
use Getopt::Long;

my $usage = "Usage: perl K-get_fasta_from_ids.pl <options>\n
Necessary options for the script to run:
-i|--ids location of the file with sequence ids
-f|--fasta location of the fasta file
-o|--output location of the output file
";
my ($fasta, $ids_file, $output_file);

GetOptions(
       'f|fasta:s'        => \$fasta,
       'i|ids:s'    => \$ids_file,
       'o|output:s'    => \$output_file,
       );
die $usage unless ( defined($fasta) and defined($ids_file) );

my %ids;
my @info;

open (my $FILEHANDLE, "<", $ids_file) or die "cannot open this file";

my $FILEOUT;

if (defined($output_file)) {

    open ($FILEOUT, ">", $output_file);
}

while (my $ids_line = <$FILEHANDLE>){

chomp $ids_line;

@info = split (/\s/, $ids_line);
$info[0] =~ s/\>//;
$info[0] =~ s/,//;
$ids{$info[0]} = 1;

}

close $FILEHANDLE;

my %sequence;
my %header; 
my @query_info;
my $query;

open (FILE1, "<", $fasta) or die "cannot open this file";

while ( my $line = <FILE1> ){
    
    chomp $line; 

if ($line=~ m/\>/)	{

@query_info = split (/\s/, $line);

$query = shift(@query_info);

$query =~ s/\>//;
$query =~ s/,//;

#for ORF Predictor
#$query =~ s/ .+//;
#use for Unigene only
$query =~ s/gnl\|UG\|//;

$header{$query} = join(" ", @query_info);

}

elsif($sequence{$query}) {

$sequence{$query} = $sequence{$query} . $line;

}
else{
$sequence{$query} = $line;

}
}
close FILE1;


my $key;
my $value;

while (($key, $value) = each(%sequence)){

if ($ids{$key}){

    if (defined($output_file)){
	print $FILEOUT ">", $key, " ", $header{$key}, "\n", $sequence{$key}, "\n";
    }

    else {
	print ">", $key, " ", $header{$key}, "\n", $sequence{$key}, "\n";
    }
}

else {

#print ">", $key, " ", $header{$key}, "\n", $sequence{$key}, "\n";

}

}


__END__
