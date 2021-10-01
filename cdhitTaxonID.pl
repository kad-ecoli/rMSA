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
cdhitTaxonID.pl seq.fasta db.gz db.gz.tsv db.filter
    filter fasta format alignment so that sequences are non-redundant
    within each species

Input:
    seq.fasta - query fasta. used to calculate length of throw_away_sequences
    db.gz     - fasta file to filter, header should start with NCBI NT
                accession or RNAcentral accession

Output:
    db.gz.tsv - mapping file from fasta header to taxon ID
    db.filter - filtered fasta file 
EOF
;

if (@ARGV<2)
{
    print $docstring;
    exit(1);
}

my $queryfile=$ARGV[0];
my $infile   =$ARGV[1];
my $tsvfile  =$ARGV[2];
my $outfile  =$ARGV[3];
my $max_aln_seqs     =200000; # max number of blastn alignmnets to parse

if (!-s "$db1tsv" && !-s "$tsvfile")
{
    print "ERROR! No mapping file $db1tsv\n";
    exit(1);
}
if (!-s "$db2" && !-s "$tsvfile")
{
    print "ERROR! No such file $db2\n";
    exit(1);
}
if (!-s "$db2.nal" && !-s "$db2.ndb" && !-s "$tsvfile")
{
    print "ERROR! $db2 not in blastn format. Please run\n";
    print "$dbdir/script/makeblastdb -in $db2 -parse_seqids -hash_index -dbtype nucl\n";
    exit(1);
}
if (!-s "$queryfile")
{
    print "ERROR! No such file $queryfile\n";
    exit(1);
}
my $sequence=`$bindir/fastaOneLine $queryfile|tail -1`;
chomp($sequence);
my $Lch=length $sequence;

&getTaxonID($infile,$tsvfile) if (!-s "$tsvfile");
# TODO: perhaps a cd-hit-est -c 1.0 and cd-hit-est-2d -c 1.0 regardless of species anyway, esp when the number of sequences are large
&cdhitTaxonID($Lch,$infile,$tsvfile,$outfile);

exit();

sub cdhitTaxonID
{
    my ($Lch,$infile,$tsvfile,$outfile)=@_;
    
    my $throw_away_sequences=int(0.4*$Lch);
    $throw_away_sequences=9 if ($throw_away_sequences<10);
    print "throw_away_sequences<=$throw_away_sequences\n";

    my %taxon2fasta_dict;
    $taxon2fasta_dict{"0"}="";
    my %taxon_count_dict;
    $taxon_count_dict{"0"}=0;
    my %header2taxon_dict;
    my @taxonIDs_list=("0");
    foreach my $line(`tail -n+2 $tsvfile`)
    {
        if ($line=~/\S+\t(\S+)\t(\S+)/)
        {
            my $header="$1";
            my $taxonIDs="$2";
            $header2taxon_dict{$header}=$taxonIDs;
            push(@taxonIDs_list,($taxonIDs)) if (!defined $taxon2fasta_dict{$taxonIDs});
            $taxon2fasta_dict{$taxonIDs}="";
        }
    }
    my $cat="cat";
    $cat="zcat" if ($infile=~/.gz/);
    foreach my $line(`$cat $infile | $bindir/fasta2pfam -`)
    {
        if ($line=~/(\S+)\t(\S+)/)
        {
            my $header="$1";
            my $sequence="$2";
            next if (length $sequence<=$throw_away_sequences);
            my $taxonIDs="0";
            $taxonIDs=$header2taxon_dict{$header} if (defined $header2taxon_dict{$header});
            $taxon2fasta_dict{$taxonIDs}.=">$header\n$sequence\n";
            $taxon_count_dict{$taxonIDs}=0 if (!defined $taxon_count_dict{$taxonIDs});
            $taxon_count_dict{$taxonIDs}+=1;
        }
    }
    my $min_group_size=2; # smallest group size to triger a cd-hit
    # 2: 26min 
    # 3: 18min 
    # 4: 14min 
    # 5: 11min
    foreach my $c((1.00,0.99))
    {
        my $total_count=0;
        foreach my $taxonIDs(@taxonIDs_list)
        {
            if ($taxon_count_dict{$taxonIDs}<$min_group_size)
            {
                $total_count+=$taxon_count_dict{$taxonIDs};
                next;
            }
            print "$taxonIDs\t$taxon_count_dict{$taxonIDs}\n";
            open(FP,">$tsvfile.tmp");
            print FP "$taxon2fasta_dict{$taxonIDs}";
            close(FP);
            # always use -T 1 because each cluster is usually small.
            # multi-threading therefore does not make much sense.
            system("$bindir/cd-hit-est -T 1 -i $tsvfile.tmp -c $c -o $tsvfile.cdhit -l $throw_away_sequences -M 5000 >/dev/null");
            $taxon2fasta_dict{$taxonIDs}=`cat $tsvfile.cdhit`;
            $taxon_count_dict{$taxonIDs}=`grep '^>' $tsvfile.cdhit|wc -l`+0;
            $total_count+=$taxon_count_dict{$taxonIDs};
            system("rm $tsvfile.tmp $tsvfile.cdhit $tsvfile.cdhit.clstr");
        }
        print "$total_count sequences at -c $c\n";
        last if ($total_count<=$max_aln_seqs);
        $min_group_size--;
        $min_group_size=2 if ($min_group_size<2);
    }
    my $txt="";
    foreach my $taxonIDs(@taxonIDs_list)
    {
        $txt.=$taxon2fasta_dict{$taxonIDs};
    }
    open(FP,">$outfile");
    print FP $txt;
    close(FP);
}

sub getTaxonID
{
    my ($infile,$tsvfile)=@_;

    my $tmpfile="$tsvfile.tmp";

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
