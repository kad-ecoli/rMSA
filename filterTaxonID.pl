#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';

my $rootdir=dirname(abs_path(__FILE__));
my $bindir ="$rootdir/bin";
my $dbdir  ="$rootdir/database";
my $db1tsv ="$dbdir/rnacentral.tsv";
my $db2    ="$dbdir/nt";
my $max_split_seqs=5000;

my $docstring=<<EOF
filterTaxonID.pl seq.afa seq.afa.tsv seq.filter.afa
    filter fasta format alignment so that on the first hit for
    each species is kept

Input:
    seq.afa - fasta file, header should start with NCBI NT accession
              or RNAcentral accession

Output:
    seq.afa.tsv    - mapping file from fasta header to taxon ID
    seq.filter.afa - filtered fasta file 
EOF
;

if (@ARGV<2)
{
    print $docstring;
    exit(1);
}

my $infile =$ARGV[0];
my $tsvfile=$ARGV[1];
my $outfile=$ARGV[2];

if (!-s "$db1tsv")
{
    print "ERROR! No mapping file $db1tsv\n";
    exit(1);
}
if (!-s "$db2")
{
    print "ERROR! No such file $db2\n";
    exit(1);
}
if (!-s "$db2.nal" && !-s "$db2.ndb")
{
    print "ERROR! $db2 not in blastn format. Please run\n";
    print "$dbdir/script/makeblastdb -in $db2 -parse_seqids -hash_index -dbtype nucl\n";
    exit(1);
}

&getTaxonID($infile,$tsvfile) if (!-s "$tsvfile");
&filterTaxonID($infile,$tsvfile,$outfile);

exit();

sub filterTaxonID
{
    my ($infile,$tsvfile,$outfile)=@_;
    my %accept_taxon_dict;
    my %reject_header_dict;
    my %accept_header_dict;
    foreach my $line(`tail -n+2 $tsvfile`)
    {
        if ($line=~/\S+\t(\S+)\t(\S+)/)
	{
	    my $header="$1";
	    my $taxonIDs="$2";
	    my $accept_taxonIDs="";
	    foreach my $taxonID(split(',',$taxonIDs))
	    {
	        next if (defined $accept_taxon_dict{$taxonID});
		$accept_taxonIDs.=",$taxonID";
		$accept_taxon_dict{$taxonID}=$header;
	    }
	    if ($accept_taxonIDs eq "")
	    {
	        $reject_header_dict{$header}=$taxonIDs;
            }
            else
	    {
	        #print "$taxonIDs\t$accept_taxonIDs\n";
	        $accept_header_dict{$header}=substr($accept_taxonIDs,1);
            }

	}
    }
    my $txt="";
    foreach my $line(`$bindir/fasta2pfam $infile`)
    {
	if ($line=~/(\S+)\t(\S+)/)
	{
	    my $header="$1";
	    my $sequence="$2";
	    next if (defined $reject_header_dict{$header});
	    my $taxonIDs="0";
	    $taxonIDs=$accept_header_dict{$header} if (defined $accept_header_dict{$header});
	    $txt.=">$header\t$taxonIDs\n$sequence\n";
	}
    }
    open(FP,">$outfile");
    print FP $txt;
    close(FP);
}

sub getTaxonID
{
    my ($infile,$tsvfile)=@_;

    my $tmpfile=$tsvfile.".tmp";

    my @header_list=();
    my @accession_list=();
    my $accession="";
    my @rc_list=();
    my @nt_list=();
    my $header;
    my $cat="cat";
    $cat="zcat" if ($infile=~/.gz/);
    foreach $header(`$cat $infile|grep '^>' |sed 's/>//g'|cut -f1`)
    {
        chomp($header);
        $accession="";
        if ($header=~/^(URS[A-Z0-9]+)/) # RNAcentral
        {
            $accession="$1";
            push(@rc_list,($accession)) if (! grep( /^$accession$/, @rc_list));
        }
        elsif ($header=~/^([_A-Z0-9]+)[.]/)
        {
            $accession="$1";
            push(@nt_list,($accession)) if (! grep( /^$accession$/, @nt_list));
        }
        elsif ($header=~/^([A-Z0-9]{4}_[A-Za-z0-9]{1,4})/) # PDB chain
        {
            $accession="$1";
            push(@nt_list,($accession)) if (! grep( /^$accession$/, @nt_list));
        }
        if (length $accession)
        {
            push(@header_list,($header));
            push(@accession_list,($accession));
        }
        else
        {
            print "skip unmappable entry >".$header."\n";
        }
    }

    printf "mapping %d RNAcentral accession(s)\n", scalar @rc_list;
    my %taxon_dict;
    %taxon_dict = map { $_ => "" } @accession_list;
    my $taxonIDs;
    foreach my $line(`$cat $infile| grep -ohP ">URS[A-Z0-9]+" | sed 's/>//g' | $bindir/getRNAcentralTaxonID - $db1tsv -`)
    {
        chomp($line);
        if ($line=~/(\S+)\t(\S+)$/)
        {
            $accession="$1";
            $taxonIDs="$2";
            $taxon_dict{$accession}=$taxonIDs;
        }
    }
 
    printf "mapping %d NCBI nucleotide (NT) accession(s)\n", scalar @nt_list;
    for (my $j=0;$j<scalar @nt_list;$j+=$max_split_seqs)
    {
        my $txt="";
        foreach (my $i=$j;$i<=$j+$max_split_seqs;$i++)
        {
            $txt.="$nt_list[$i]\n";
        }
        open(FP,">$tmpfile");
        print FP $txt;
        close(FP);
        foreach my $line(`$bindir/blastdbcmd -db $db2 -entry_batch $tmpfile -outfmt '%a %T'`)
        {
            if ($line=~/(\S+)\s(\S+)/)
            {
                $accession="$1";
                $taxonIDs="$2";
                $accession="$1" if ($accession=~/(\S+)[.]\d+/);
                $taxon_dict{$accession}=$taxonIDs;
            }
        }
    }
    system("rm $tmpfile");

    printf "writing mapping file for %d hits\n", scalar @header_list;
    my $txt="#accession\thit\ttaxonID\n";
    my $names;
    for (my $i=0;$i<scalar @header_list;$i++)
    {
        $header=$header_list[$i];
        $accession=$accession_list[$i];
        $taxonIDs=$taxon_dict{$accession};
        if (length $taxonIDs==0)
        {
            print "failed to map >$header\n";
        }
        $txt.="$accession\t$header\t$taxonIDs\n";
    }
    open(FP,">$tsvfile");
    print FP $txt;
    close(FP);
}
