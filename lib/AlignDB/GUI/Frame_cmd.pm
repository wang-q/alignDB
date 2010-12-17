package AlignDB::GUI::Frame;
use strict;
use warnings;

use Proc::Background;

#----------------------------#
# execute system commands
#----------------------------#
sub exec_cmd {
    my $self = shift;
    my $cmd  = shift;

    print "\n", "=" x 12, "CMD", "=" x 15, "\n";
    print $cmd , "\n";
    print "=" x 30, "\n";

    $cmd .= ' 1>&2';    # redirect STDOUT to STDERR so we can see the outputs

    my $proc = Proc::Background->new($cmd);
    $proc->{cmd} = $cmd;

    my $processes = $self->processes;
    if ( ref $processes eq 'ARRAY' ) {
        push @$processes, $proc;
    }
    else {
        $processes = [$proc];
    }
    $self->processes($processes);
    
    #$proc->wait;

    return;
}

#----------------------------#
# run *.pl events
#----------------------------#
sub event_init_alignDB {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $db_name  = $self->get_value("db_name");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");
    my $cmd
        = "perl $FindBin::Bin/../init/init_alignDB.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$db_name"
        . " -u=$username"
        . " --password=$password";
    $self->exec_cmd($cmd);
    return;
}

sub event_ref_outgroup {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $first_db  = $self->get_value("first_db");
    my $second_db = $self->get_value("second_db");
    my $goal_db   = $self->get_value("goal_db");

    my $first    = $self->get_value("first");
    my $second   = $self->get_value("second");
    my $outgroup = $self->get_value("outgroup");

    my $cmd
        = "perl $FindBin::Bin/../extra/ref_outgroup.pl"
        . " --server=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " --first_db=$first_db"
        . " --second_db=$second_db"
        . " --goal_db=$goal_db"
        . " --target=$first"
        . " --query=$second"
        . " --outgroup=$outgroup"
        . " --chr_id_runlist=1-1000";
    $self->exec_cmd($cmd);
    return;
}

sub event_gen_alignDB {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $db_name  = $self->get_value("db_name");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $target_taxon_id = $self->get_value("target_taxon_id");
    my $target_name     = $self->get_value("target_name");
    my $query_taxon_id  = $self->get_value("query_taxon_id");
    my $query_name      = $self->get_value("query_name");
    my $axt_dir         = $self->get_value("axt_dir");
    my $axt_threshold   = $self->get_value("axt_threshold");
    my $insert_dG       = $self->get_value("insert_dG") || 0;
    my $parallel        = $self->get_value("parallel");

    my $cmd
        = "perl $FindBin::Bin/../init/gen_alignDB.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$db_name"
        . " -u=$username"
        . " --password=$password"
        . " -t=\"$target_taxon_id,$target_name\""
        . " -q=\"$query_taxon_id,$query_name\""
        . " -a=$axt_dir"
        . " -l=$axt_threshold"
        . " --insert_dG=$insert_dG"
        . " --parallel=$parallel";

    $self->exec_cmd($cmd);
    return;
}

sub event_insert_gc {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $db_name  = $self->get_value("db_name");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $insert_gc      = $self->get_value("insert_gc")      || 0;
    my $insert_segment = $self->get_value("insert_segment") || 0;
    my $parallel       = $self->get_value("parallel");

    my $cmd
        = "perl $FindBin::Bin/../init/insert_gc.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$db_name"
        . " -u=$username"
        . " --password=$password"
        . " --insert_gc=$insert_gc"
        . " --insert_segment=$insert_segment"
        . " --parallel=$parallel";

    $self->exec_cmd($cmd);
    return;
}

sub event_insert_gene {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $db_name  = $self->get_value("db_name");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $ensembl = $self->get_value("gene_ensembl");

    my $insert_genesw   = $self->get_value("insert_genesw")   || 0;
    my $insert_exonsw   = $self->get_value("insert_exonsw")   || 0;
    my $insert_codingsw = $self->get_value("insert_codingsw") || 0;

    my $cmd
        = "perl $FindBin::Bin/../gene/insert_gene.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$db_name"
        . " -u=$username"
        . " --password=$password"
        . " -e=$ensembl"
        . " --insert_genesw=$insert_genesw"
        . " --insert_exonsw=$insert_exonsw"
        . " --insert_codingsw=$insert_codingsw";

    $self->exec_cmd($cmd);
    return;
}

sub event_update_feature {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $db_name  = $self->get_value("db_name");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $ensembl = $self->get_value("ensembl");

    my $process_align  = $self->get_value("process_align")  || 0;
    my $process_indel  = $self->get_value("process_indel")  || 0;
    my $process_isw    = $self->get_value("process_isw")    || 0;
    my $process_snp    = $self->get_value("process_snp")    || 0;
    my $process_window = $self->get_value("process_window") || 0;
    my $parallel       = $self->get_value("parallel");

    my $cmd
        = "perl $FindBin::Bin/../init/update_feature.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$db_name"
        . " -u=$username"
        . " --password=$password"
        . " -e=$ensembl"
        . " --process_align=$process_align"
        . " --process_indel=$process_indel"
        . " --process_isw=$process_isw"
        . " --process_snp=$process_snp"
        . " --process_window=$process_window"
        . " --parallel=$parallel";

    $self->exec_cmd($cmd);
    return;
}

sub event_update_indel_slippage {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $db_name  = $self->get_value("db_name");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $cmd
        = "perl $FindBin::Bin/../init/update_indel_slippage.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$db_name"
        . " -u=$username"
        . " --password=$password";

    $self->exec_cmd($cmd);
    return;
}

sub event_update_isw_indel_id {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $db_name  = $self->get_value("db_name");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $cmd
        = "perl $FindBin::Bin/../init/update_isw_indel_id.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$db_name"
        . " -u=$username"
        . " --password=$password";

    $self->exec_cmd($cmd);
    return;
}

sub event_update_segment {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $db_name  = $self->get_value("db_name");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $cmd
        = "perl $FindBin::Bin/../init/update_segment.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$db_name"
        . " -u=$username"
        . " --password=$password";

    $self->exec_cmd($cmd);
    return;
}

sub event_update_indel_occured {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");
    my $aim_db   = $self->get_value("aim_db");
    my $ref_db   = $self->get_value("ref_db");

    my $cmd
        = "perl $FindBin::Bin/../extra/update_indel_occured.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " --aim_db=$aim_db"
        . " --ref_db=$ref_db";

    $self->exec_cmd($cmd);
    return;
}

sub event_update_isw_dxr {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");
    my $goal_db  = $self->get_value("goal_db");

    my $cmd
        = "perl $FindBin::Bin/../extra/update_isw_dxr.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$goal_db"
        . " -u=$username"
        . " --password=$password";

    $self->exec_cmd($cmd);
    return;
}

sub event_update_snp_cpg {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $db_name  = $self->get_value("db_name");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $cmd
        = "perl $FindBin::Bin/../extra/update_snp_cpg.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$db_name"
        . " -u=$username"
        . " --password=$password";

    $self->exec_cmd($cmd);
    return;
}

sub event_apply_sql {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $db_name  = $self->get_value("db_name");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $sql_file = $self->get_value("sql_file");

    my $cmd
        = "perl $FindBin::Bin/../util/apply_sql.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$db_name"
        . " -u=$username"
        . " --password=$password"
        . " -f=$sql_file";

    $self->exec_cmd($cmd);
    return;
}

sub event_query_sql {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $db_name  = $self->get_value("db_name");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $query_output      = $self->get_value("query_output");
    my $query_output_type = $self->query_output_type;
    my $sql_query         = $self->get_value("query_sql");
    $sql_query =~ s/\s+/ /g;

    my $cmd
        = "perl $FindBin::Bin/../util/query_sql.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$db_name"
        . " -u=$username"
        . " --password=$password"
        . " -o=$query_output"
        . " -t=$query_output_type"
        . qq{ --query="$sql_query"};

    $self->exec_cmd($cmd);
    return;
}

sub event_common_stat {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $db_name  = $self->get_value("db_name");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $output    = $self->get_value("common_stat_file");
    my $run       = $self->get_value("common_run");
    my $threshold = $self->get_value("common_threshold");

    my $cmd
        = "perl $FindBin::Bin/../stat/common_stat_factory.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$db_name"
        . " -u=$username"
        . " --password=$password"
        . " -o=$output"
        . " -r=$run"
        . " -t=$threshold";

    $self->exec_cmd($cmd);
    return;
}

sub event_common_chart {
    my $self  = shift;
    my $event = shift;

    my $common_stat_file = $self->get_value("common_stat_file");
    return if $^O ne "MSWin32";
    return if !$common_stat_file;

    my $jc_correction   = $self->get_value("common_jc")              || 0;
    my $time_stamp      = $self->get_value("common_time_stamp")      || 0;
    my $add_index_sheet = $self->get_value("common_add_index_sheet") || 0;

    my $cmd
        = "perl $FindBin::Bin/../stat/common_chart_factory.pl"
        . " -i=$common_stat_file"
        . " -j=$jc_correction"
        . " -t=$time_stamp"
        . " -a=$add_index_sheet";

    $self->exec_cmd($cmd);
    return;
}

sub event_gc_stat {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $db_name  = $self->get_value("db_name");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $output    = $self->get_value("gc_stat_file");
    my $run       = $self->get_value("gc_run");
    my $threshold = $self->get_value("gc_threshold");

    my $cmd
        = "perl $FindBin::Bin/../stat/gc_stat_factory.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$db_name"
        . " -u=$username"
        . " --password=$password"
        . " -o=$output"
        . " -r=$run"
        . " -t=$threshold";

    $self->exec_cmd($cmd);
    return;
}

sub event_gc_chart {
    my $self  = shift;
    my $event = shift;

    my $gc_stat_file = $self->get_value("gc_stat_file");
    return if $^O ne "MSWin32";
    return if !$gc_stat_file;

    my $jc_correction   = $self->get_value("gc_jc")              || 0;
    my $time_stamp      = $self->get_value("gc_time_stamp")      || 0;
    my $add_index_sheet = $self->get_value("gc_add_index_sheet") || 0;

    my $cmd
        = "perl $FindBin::Bin/../stat/gc_chart_factory.pl"
        . " -i=$gc_stat_file"
        . " -j=$jc_correction"
        . " -t=$time_stamp"
        . " -a=$add_index_sheet";

    $self->exec_cmd($cmd);
    return;
}

sub event_gene_stat {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $db_name  = $self->get_value("db_name");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $output = $self->get_value("gene_stat_file");
    my $run    = $self->get_value("gene_run");

    my $cmd
        = "perl $FindBin::Bin/../stat/gene_stat_factory.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$db_name"
        . " -u=$username"
        . " --password=$password"
        . " -o=$output"
        . " -r=$run";

    $self->exec_cmd($cmd);
    return;
}

sub event_gene_chart {
    my $self  = shift;
    my $event = shift;

    my $gene_stat_file = $self->get_value("gene_stat_file");
    return if $^O ne "MSWin32";
    return if !$gene_stat_file;

    my $jc_correction   = $self->get_value("gene_jc")              || 0;
    my $time_stamp      = $self->get_value("gene_time_stamp")      || 0;
    my $add_index_sheet = $self->get_value("gene_add_index_sheet") || 0;

    my $cmd
        = "perl $FindBin::Bin/../stat/gene_chart_factory.pl"
        . " -i=$gene_stat_file"
        . " -j=$jc_correction"
        . " -t=$time_stamp"
        . " -a=$add_index_sheet";

    $self->exec_cmd($cmd);
    return;
}

sub event_three_stat {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");
    my $goal_db  = $self->get_value("goal_db");

    my $output    = $self->get_value("three_stat_file");
    my $run       = $self->get_value("three_run");
    my $threshold = $self->get_value("three_threshold");

    my $cmd
        = "perl $FindBin::Bin/../stat/three_stat_factory.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$goal_db"
        . " -u=$username"
        . " --password=$password"
        . " -o=$output"
        . " -r=$run"
        . " -t=$threshold";

    $self->exec_cmd($cmd);
    return;
}

sub event_three_chart {
    my $self  = shift;
    my $event = shift;

    my $three_stat_file = $self->get_value("three_stat_file");
    return if $^O ne "MSWin32";
    return if !$three_stat_file;

    my $jc_correction   = $self->get_value("three_jc")              || 0;
    my $time_stamp      = $self->get_value("three_time_stamp")      || 0;
    my $add_index_sheet = $self->get_value("three_add_index_sheet") || 0;

    my $cmd
        = "perl $FindBin::Bin/../stat/three_chart_factory.pl"
        . " -i=$three_stat_file"
        . " -j=$jc_correction"
        . " -t=$time_stamp"
        . " -a=$add_index_sheet";

    $self->exec_cmd($cmd);
    return;
}

1;
