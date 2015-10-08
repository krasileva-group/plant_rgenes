#!/usr/bin/perl -w

# Author: Ksenia V Krasileva

# This script parses the output of K-parse_Pfam_domains-v3.1.pl
# Domains from two files -1 and -2 are analysed and the overlap between the lists is printed
# Usage case:
# This script can be used to calculate overlap between NLR-fusions domains and domains of proteins targeted by effectors
# Set any domains to skip at line 112

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $today=`date +'%m%d%Y'`;

chomp $today;

my ($file_1, $file_2, $output_dir);
my $evalue_cutoff=1e-3;
    
GetOptions(
       '1|file1:s'        => \$file_1,
       'e|evalue:s'        => \$evalue_cutoff,
       '2|file2:s'        => \$file_2,
       'o|output_dir:s'   => \$output_dir, 
);  

my $usage = "usage: perl script.pl <options>
-1|--file1 first list of genes with parsed domains
-2|--file2 second list of genes with parsed domains  
-e|--evalue evalue cutoff for determining domain fusions [default 1e-3]
-o|--output output directory
";
die $usage unless ( (defined $file_1) and (defined $file_2) );

my (%domains_1, %domains_2);
my (%domains_list_1, %domains_list_2);

my $basefile_1=`basename $file_1 .parsed`;
chomp $basefile_1;
my $basefile_2=`basename $file_2 .parsed`;
chomp $basefile_2;

my $output_1 = $output_dir . "/" . $basefile_1 . "_" . $basefile_2 . ".overlap";
my $output_2 = $output_dir . "/" . $basefile_2 . "_" . $basefile_1 . ".overlap";
    
open (my $FILE1, "<", $file_1) or die $usage;
    
    # sample input
    # Traes_4DS_92264B9F4.1LRR_4(start=58, stop=98, evalue=0.0021)~DARPP-32(start=145, stop=241, evalue=0.015)
	
    while (my $line = <$FILE1> ){
    
	chomp $line; 

	my ($seqid, $domainstring)=split("\t", $line);
	    
	    $domains_1{$seqid}=$domainstring; #record all sequences and their associated domains
	    
	    my @domains = split('~', $domainstring);

	    foreach (@domains){

		my $domain=$_;

		my ($domainid, $attributes)=split(/\(/, $domain);
		$attributes=~ s/\)//; #remove right hand bracket
		$attributes=~ s/\s//g;
		my ($start_strng, $stop_string, $evalue_string)=split(",", $attributes);
		my ($foo, $evalue)=split("=", $evalue_string);

		if ($evalue < $evalue_cutoff){
		    
		    unless ($domainid=~ m/NB-ARC|LRR|AAA|AAA.+|TIR|RPW8/ ){ #skip known domains
			
			push @{$domains_list_1{$domainid}}, $seqid; #make a list of sequences for each domain
		    }
		}
	    }

	}

close $FILE1;

open (my $FILE2, "<", $file_2) or die $usage;

while (my $line = <$FILE2> ){

    chomp $line;

    my ($seqid, $domainstring)=split("\t", $line);

    $domains_2{$seqid}=$domainstring; #record all NLRs for a species                                                                                  

    my @domains = split('~', $domainstring);

    foreach (@domains){

	my $domain=$_;

	my ($domainid, $attributes)=split(/\(/, $domain);
        $attributes=~ s/\)//; #remove right hand bracket                                                                                              
        $attributes=~ s/\s//g;
        my ($start_strng, $stop_string, $evalue_string)=split(",", $attributes);
        my ($foo, $evalue)=split("=", $evalue_string);

        if ($evalue < $evalue_cutoff){

	    unless ($domainid=~ m/NB-ARC|LRR|AAA|AAA.+|TIR|RPW8/ ){ #skip known domains  

	    push @{$domains_list_2{$domainid}}, $seqid; #make a list of sequences for each domain                                                 

	    }
        }
    }

}



# create output files and a summary table

open SUMMARYFILE, ">", "overlap_domains_summary" . $basefile_1 . "_" . $basefile_2 . "-" . $today . ".tsv"; 
open FILEOUT1, ">", $output_1;
open FILEOUT2, ">", $output_2;

my $overlap_counter=0;

foreach my $domain (sort keys (%domains_list_1) ){


    if ( exists $domains_list_2{$domain}){

	print $domain, "\n";

	$overlap_counter++;

	print SUMMARYFILE $domain, "\n";

	foreach(@{$domains_list_1{$domain}}){

	    print FILEOUT1 $_, "\t", $domains_1{$_}, "\n";
		}

        foreach(@{$domains_list_2{$domain}}){

	    print FILEOUT2 $_, "\t", $domains_2{$_}, "\n";
	    
	}

    }
}

close SUMMARYFILE;

    

sub uniq {

    my %seen;
    grep !$seen{$_}++, @_;

}

__END__
