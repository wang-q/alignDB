package AlignDB::GUI::Frame;
use strict;
use warnings;

use Config::Tiny;
use FindBin;

#----------------------------#
# read-out configs
#----------------------------#
sub read_config {
    my $self = shift;

    my $Config = Config::Tiny->new;
    $Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

    # Database init values
    $self->set_value( "server",   $Config->{database}->{server} );
    $self->set_value( "port",     $Config->{database}->{port} );
    $self->set_value( "username", $Config->{database}->{username} );
    $self->set_value( "password", $Config->{database}->{password} );

    # target, query init values
    $self->set_value( "target_taxon_id",
        $Config->{taxon}->{target_taxon_id} );
    $self->set_value( "target_name",    $Config->{taxon}->{target_name} );
    $self->set_value( "query_taxon_id", $Config->{taxon}->{query_taxon_id} );
    $self->set_value( "query_name",     $Config->{taxon}->{query_name} );

    $self->set_value( "db_name",      $Config->{database}->{db} );
    $self->set_value( "ensembl",      $Config->{database}->{ensembl} );
    $self->set_value( "gene_ensembl", $Config->{database}->{ensembl} );

    # generate
    $self->set_value( "axt_dir",       $Config->{taxon}->{axt_dir} );
    $self->set_value( "insert_dG",     $Config->{generate}->{insert_dG} );
    $self->set_value( "axt_threshold", $Config->{generate}->{axt_threshold} );
    $self->set_value( "parallel",      $Config->{generate}->{parallel} );

    # insert GC
    $self->set_value( "insert_gc",      $Config->{gc}->{insert_gc} );
    $self->set_value( "insert_segment", $Config->{gc}->{insert_segment} );

    # insert gene
    $self->set_value( "insert_genesw",   $Config->{gene}->{insert_genesw} );
    $self->set_value( "insert_exonsw",   $Config->{gene}->{insert_exonsw} );
    $self->set_value( "insert_codingsw", $Config->{gene}->{insert_codingsw} );

    # update feature
    $self->set_value( "process_align",  $Config->{feature}->{align} );
    $self->set_value( "process_indel",  $Config->{feature}->{indel} );
    $self->set_value( "process_isw",    $Config->{feature}->{isw} );
    $self->set_value( "process_snp",    $Config->{feature}->{snp} );
    $self->set_value( "process_window", $Config->{feature}->{window} );

    # common stat parameter
    $self->set_value( "common_run",        $Config->{stat}->{run} );
    $self->set_value( "common_threshold",  $Config->{stat}->{sum_threshold} );
    $self->set_value( "common_jc",         $Config->{stat}->{jc_correction} );
    $self->set_value( "common_time_stamp", $Config->{stat}->{time_stamp} );
    $self->set_value( "common_add_index_sheet",
        $Config->{stat}->{add_index_sheet} );
    
    # gc stat parameter
    $self->set_value( "gc_run",        $Config->{stat}->{run} );
    $self->set_value( "gc_threshold",  $Config->{stat}->{sum_threshold} );
    $self->set_value( "gc_jc",         $Config->{stat}->{jc_correction} );
    $self->set_value( "gc_time_stamp", $Config->{stat}->{time_stamp} );
    $self->set_value( "gc_add_index_sheet",
        $Config->{stat}->{add_index_sheet} );
    
    # gene stat parameter
    $self->set_value( "gene_run",        $Config->{stat}->{run} );
    $self->set_value( "gene_jc",         $Config->{stat}->{jc_correction} );
    $self->set_value( "gene_time_stamp", $Config->{stat}->{time_stamp} );
    $self->set_value( "gene_add_index_sheet",
        $Config->{stat}->{add_index_sheet} );
    
    # three stat parameter
    $self->set_value( "three_run",        $Config->{stat}->{run} );
    $self->set_value( "three_threshold",  $Config->{stat}->{sum_threshold} );
    $self->set_value( "three_jc",         $Config->{stat}->{jc_correction} );
    $self->set_value( "three_time_stamp", $Config->{stat}->{time_stamp} );
    $self->set_value( "three_add_index_sheet",
        $Config->{stat}->{add_index_sheet} );

    return;
}

1;
