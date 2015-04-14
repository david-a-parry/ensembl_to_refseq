#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use FindBin;
#use lib "$FindBin::Bin/bioperl-1.2.3";
use lib "$FindBin::Bin/lib/modules";
use Bio::EnsEMBL::Registry;

my @symbols;
my @gene_ids;
my @transcript_ids;
#my @refseq_ids;
my $species = 'Human';
my $help;
GetOptions(
            "symbol=s{,}"       => \@symbols,
            "gene=s{,}"         => \@gene_ids,
            "transcript=s{,}"   => \@transcript_ids,
#            "refseq=s{,}"       => \@refseq_ids,
            "organism=s"        => \$species,
            "help"              => \$help,
            ) or usage("Syntax error");
usage() if $help;
if (not @gene_ids 
    and not @transcript_ids 
    and not @symbols 
    #and not @refseq_ids
){
    usage("At least one --symbol, --gene or --transcript must be specified.") ;
}

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',    # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

my @genes = ();
my @transcripts = ();

if (@gene_ids or @symbols){
    my $gene_adaptor  = $registry->get_adaptor( $species, 'Core', 'Gene' );
    foreach my $s (@symbols){
        my $g = $gene_adaptor->fetch_by_display_label($s);
        if (not defined $g){
            print STDERR "$s not found as gene symbol for $species.\n";
        }else{
            push @genes, $g; 
        }
    }
    foreach my $g_id (@gene_ids){
        my $g = $gene_adaptor->fetch_by_stable_id($g_id);
        if (not defined $g){
            print STDERR "$g_id not found as gene ID for $species.\n";
        }else{
            push @genes, $g; 
        }
    }
}

if (@transcript_ids){
    my $transcript_adaptor =
      $registry->get_adaptor( $species, 'Core', 'Transcript' );
    foreach my $t (@transcript_ids){
        my $g = $transcript_adaptor->fetch_by_stable_id($t);
        if (not defined $g){
            print STDERR "$t not found as transcript ID for $species.\n";
        }else{
            push @transcripts, $g; 
        }
    }
}

#if (@refseq_ids){
#    my $db_entry_adaptor =
#        $registry->get_adaptor( $species, 'Core', 'DBEntry' );
#    foreach my $r (@refseq_ids){
#        my $db_entry = $db_entry_adaptor->fetch_by_db_accession('RefSeq_mRNA', $r);   
#    }
#}


exit if not @genes and not @transcripts; 
print join("\t", qw( 
                SYMBOL 
                ENSEMBL_GENE_ID 
                ENSEMBL_TRANSCRIPT 
                ENSEMBL_PROTEIN 
                REFSEQ_TRANSCRIPT
                REFSEQ_PROTEIN
                    )
            ) ."\n"; 

foreach my $g (@genes){
    foreach my $transcript ( @{ $g->get_all_Transcripts() } ) {
        print_transcript($transcript, $g);
    }
}
foreach my $transcript (@transcripts){
    print_transcript($transcript);
}


sub print_transcript{
    my $transcript = shift;
    my $g = shift; 
    my $ensp_id = "-";
    my $refp_id = "-";
    my $refn_id = "-";
    if (not defined $g){
        $g = $transcript->get_Gene();
    }
    my @n_ids = ();
    foreach my $dbe ( @{ $transcript->get_all_DBEntries()} ){
        if ($dbe->dbname() eq 'RefSeq_mRNA'){
            push @n_ids, $dbe->display_id();
        }
    }
    $refn_id = join(",", @n_ids) if @n_ids; 

    # Watch out: pseudogenes have no translation
    if ( defined $transcript->translation() ) {
        my $translation = $transcript->translation();
        $ensp_id = $translation->stable_id();
        my @p_ids = ();
        foreach my $dbe ( @{ $translation->get_all_DBEntries() }){
            if ($dbe->dbname() eq 'RefSeq_peptide'){
                push @p_ids, $dbe->display_id();
            }
        }
        $refp_id = join(",", @p_ids) if @p_ids; 
    }
    print join("\t", (  $g->external_name,
                        $g->stable_id(),
                        $transcript->stable_id(),
                        $ensp_id,
                        $refn_id,
                        $refp_id, 
                        )
                    ). "\n"; 
}

sub usage{
    my $msg = shift; 
    if ($msg){
        print "\t$msg\n\n";
    }
    print <<EOT;

    Cross-references Ensembl and RefSeq transcript and protein IDs.

    $0 [-s <gene symbols>] [-g <ensembl gene IDs>] [-t <ensembl transcript IDs>]

    Arguments: 
    
    -s    --symbol      <one or more gene symbols to search>
    -g    --gene        <one or more Ensembl gene IDs to search>
    -t    --transcript  <one or more Ensembl transcript IDs to search>
    -o    --organism    <name of the organism to search. Default = Human>
    -h    --help        <show this help message>

EOT
;
exit 1 if $msg;
exit;
}
