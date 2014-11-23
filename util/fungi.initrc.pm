{    # Aspergillus clavatus
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Aspergillus clavatus',
        -group   => 'core',
        -dbname  => 'aspergillus_clavatus_core_10_63_1',
    );

    my @aliases = (
        'acla',
        'Aspergillus_clavatus',
        'Aclavatus',
        'A_clavatus',
        'acla_core_63',
        'acla_63',
        'Acla',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Aspergillus clavatus',
        -alias   => \@aliases,
    );
}

{    # Aspergillus flavus
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Aspergillus flavus',
        -group   => 'core',
        -dbname  => 'aspergillus_flavus_core_10_63_1',
    );

    my @aliases = (
        'afla',
        'Aspergillus_flavus',
        'Aflavus',
        'A_flavus',
        'afla_core_63',
        'afla_63',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Aspergillus flavus',
        -alias   => \@aliases,
    );
}

{    # Aspergillus fumigatus
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Aspergillus fumigatus',
        -group   => 'core',
        -dbname  => 'aspergillus_fumigatus_core_10_63_2',
    );

    my @aliases = (
        'afum',
        'Aspergillus_fumigatus',
        'Afumigatus',
        'A_fumigatus',
        'afum_core_63',
        'afum_63',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Aspergillus fumigatus',
        -alias   => \@aliases,
    );
}

{    # Aspergillus nidulans
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Aspergillus nidulans',
        -group   => 'core',
        -dbname  => 'aspergillus_nidulans_core_10_63_6',
    );

    my @aliases = (
        'anid',
        'Aspergillus_nidulans',
        'Anidulans',
        'A_nidulans',
        'anid_core_63',
        'anid_63',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Aspergillus nidulans',
        -alias   => \@aliases,
    );
}

{    # Aspergillus niger
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Aspergillus niger',
        -group   => 'core',
        -dbname  => 'aspergillus_niger_core_10_63_1',
    );

    my @aliases = (
        'anig',
        'Aspergillus_niger',
        'Aniger',
        'A_niger',
        'anig_core_63',
        'anig_63',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Aspergillus niger',
        -alias   => \@aliases,
    );
}

{    # Aspergillus oryzae
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Aspergillus oryzae',
        -group   => 'core',
        -dbname  => 'aspergillus_oryzae_core_10_63_2',
    );

    my @aliases = (
        'aory',
        'Aspergillus_oryzae',
        'Aoryzae',
        'A_oryzae',
        'aory_core_63',
        'aory_63',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Aspergillus oryzae',
        -alias   => \@aliases,
    );
}

{    # Aspergillus terreus
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Aspergillus terreus',
        -group   => 'core',
        -dbname  => 'aspergillus_terreus_core_10_63_1',
    );

    my @aliases = (
        'ater',
        'Aspergillus_terreus',
        'Aterreus',
        'A_terreus',
        'ater_core_63',
        'ater_63',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Aspergillus terreus',
        -alias   => \@aliases,
    );
}

{    # Fungi
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Fungi',
        -group   => 'compara',
        -dbname  => 'ensembl_compara_fungi_10_63',
    );

    my @aliases = (
        'fungi',
        'Fungi',
        'fungi_compara_63',
        'compara_fungi',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Fungi',
        -alias   => \@aliases,
    );
}

{    # Neosartorya fischeri
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Neosartorya fischeri',
        -group   => 'core',
        -dbname  => 'neosartorya_fischeri_core_10_63_1',
    );

    my @aliases = (
        'nfis',
        'Neosartorya_fischeri',
        'Nfischeri',
        'N_fischeri',
        'nfis_core_63',
        'nfis_63',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Neosartorya fischeri',
        -alias   => \@aliases,
    );
}

{    # Neurospora crassa
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Neurospora crassa',
        -group   => 'core',
        -dbname  => 'neurospora_crassa_core_10_63_1',
    );

    my @aliases = (
        'ncra',
        'Neurospora_crassa',
        'Ncrassa',
        'N_crassa',
        'ncra_core_63',
        'ncra_63',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Neurospora crassa',
        -alias   => \@aliases,
    );
}

{    # Saccharomyces cerevisiae
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Saccharomyces cerevisiae',
        -group   => 'core',
        -dbname  => 'saccharomyces_cerevisiae_core_10_63_3',
    );

    my @aliases = (
        'scer',
        'Saccharomyces_cerevisiae',
        'Scerevisiae',
        'S_cerevisiae',
        'scer_core_63',
        'scer_63',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Saccharomyces cerevisiae',
        -alias   => \@aliases,
    );
}

{    # Saccharomyces cerevisiae
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Saccharomyces cerevisiae',
        -group   => 'funcgen',
        -dbname  => 'saccharomyces_cerevisiae_funcgen_10_63_3',
    );

    my @aliases = (
        'scer',
        'Saccharomyces_cerevisiae',
        'Scerevisiae',
        'S_cerevisiae',
        'scer_funcgen_63',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Saccharomyces cerevisiae',
        -alias   => \@aliases,
    );
}

{    # Saccharomyces cerevisiae
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Saccharomyces cerevisiae',
        -group   => 'otherfeatures',
        -dbname  => 'saccharomyces_cerevisiae_otherfeatures_10_63_3',
    );

    my @aliases = (
        'scer',
        'Saccharomyces_cerevisiae',
        'Scerevisiae',
        'S_cerevisiae',
        'scer_otherfeatures_63',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Saccharomyces cerevisiae',
        -alias   => \@aliases,
    );
}

{    # Saccharomyces cerevisiae
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Saccharomyces cerevisiae',
        -group   => 'variation',
        -dbname  => 'saccharomyces_cerevisiae_variation_10_63_3',
    );

    my @aliases = (
        'scer',
        'Saccharomyces_cerevisiae',
        'Scerevisiae',
        'S_cerevisiae',
        'scer_variation_63',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Saccharomyces cerevisiae',
        -alias   => \@aliases,
    );
}

{    # Schizosaccharomyces pombe
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Schizosaccharomyces pombe',
        -group   => 'core',
        -dbname  => 'schizosaccharomyces_pombe_core_10_63_1',
    );

    my @aliases = (
        'spom',
        'Schizosaccharomyces_pombe',
        'Spombe',
        'S_pombe',
        'spom_core_63',
        'spom_63',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Schizosaccharomyces pombe',
        -alias   => \@aliases,
    );
}

