use strict;
use warnings;
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

my @aliases;
my $host = 'localhost';
my $user = 'alignDB';
my $pass = 'alignDB';
my $port = 3306;

{    # human
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Homo sapiens',
        -group   => 'core',
        -dbname  => 'human_54',
    );

    @aliases = ( 'Homo_sapiens', 'H_sapiens', 'human', 'human_54');

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Homo sapiens',
        -alias   => \@aliases
    );
}

{    # yeast
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Saccharomyces cerevisiae',
        -group   => 'core',
        -dbname  => 'yeast_65',
    );

    @aliases
        = ( 'Saccharomyces_cerevisiae', 'S_cerevisiae', 'yeast', 'S288C', 'yeast_65');

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Saccharomyces cerevisiae',
        -alias   => \@aliases
    );
}

{    # mouse
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Mus musculus',
        -group   => 'core',
        -dbname  => 'mouse_65',
    );

    @aliases = ( 'M_musculus', 'Mus_musculus', 'mouse', 'mouse_65');

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Mus musculus',
        -alias   => \@aliases
    );
}

{    # arabidppsis
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Arabidopsis thaliana',
        -group   => 'core',
        -dbname  => 'ath_65',
    );

    @aliases
        = ( 'A_thaliana', 'Arabidopsis_thaliana', 'arabidppsis', 'ath', 'ath_65');

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Arabidopsis thaliana',
        -alias   => \@aliases
    );
}

{    # compara
    Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Compara',
        -dbname  => 'compara_65',
    );

    @aliases = ( 'ensembl_compara', 'compara', 'compara_65', );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Compara',
        -alias   => \@aliases
    );
}

1;
