/*==============================================================*/
/* DBMS name:      MySQL 4.0                                    */
/* Created on:     12/10/16 23:59:34                            */
/*==============================================================*/

/*==============================================================*/
/* Table: align                                                 */
/*==============================================================*/
create table align
(
   align_id                       int                            not null AUTO_INCREMENT,
   align_length                   int,
   align_comparables              int,
   align_identities               int,
   align_differences              int,
   align_gaps                     int,
   align_ns                       int,
   align_error                    int,
   align_pi                       double,
   align_indels                   int,
   align_gc                       double,
   align_average_gc               double,
   align_comparable_runlist       text,
   align_indel_runlist            text,
   align_coding                   double,
   align_repeats                  double,
   align_coding_runlist           text,
   align_repeats_runlist          text,
   primary key (align_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Table: chromosome                                            */
/*==============================================================*/
create table chromosome
(
   chr_id                         int                            not null AUTO_INCREMENT,
   common_name                    char(64),
   taxon_id                       int,
   chr_name                       char(64),
   chr_length                     int,
   primary key (chr_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Table: exon                                                  */
/*==============================================================*/
create table exon
(
   exon_id                        int                            not null AUTO_INCREMENT,
   prev_exon_id                   int                            not null,
   window_id                      int                            not null,
   gene_id                        int                            not null,
   exon_stable_id                 char(64)                       not null,
   exon_strand                    char(1),
   exon_phase                     int,
   exon_end_phase                 int,
   exon_frame                     int,
   exon_rank                      int,
   exon_is_full                   int,
   exon_tl_runlist                text,
   exon_seq                       text,
   exon_peptide                   text,
   primary key (exon_id, prev_exon_id),
   key AK_exon_stable_id (exon_stable_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: window_exon_FK                                        */
/*==============================================================*/
create index window_exon_FK on exon
(
   window_id
);

/*==============================================================*/
/* Index: gene_exon_FK                                          */
/*==============================================================*/
create index gene_exon_FK on exon
(
   gene_id
);

/*==============================================================*/
/* Table: extreme                                               */
/*==============================================================*/
create table extreme
(
   extreme_id                     int                            not null AUTO_INCREMENT,
   prev_extreme_id                int                            not null,
   window_id                      int                            not null,
   extreme_type                   char(8),
   extreme_left_amplitude         double,
   extreme_right_amplitude        double,
   extreme_left_wave_length       int,
   extreme_right_wave_length      int,
   primary key (extreme_id, prev_extreme_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: window_extreme_FK                                     */
/*==============================================================*/
create index window_extreme_FK on extreme
(
   window_id
);

/*==============================================================*/
/* Table: gene                                                  */
/*==============================================================*/
create table gene
(
   gene_id                        int                            not null AUTO_INCREMENT,
   window_id                      int                            not null,
   gene_stable_id                 char(64)                       not null,
   gene_external_name             char(64),
   gene_biotype                   char(64),
   gene_strand                    char(1),
   gene_is_full                   int,
   gene_is_known                  int,
   gene_multitrans                int,
   gene_multiexons                int,
   gene_tc_runlist                text,
   gene_tl_runlist                text,
   gene_description               text,
   primary key (gene_id),
   key AK_gene_stable_id (gene_stable_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: window_gene_FK                                        */
/*==============================================================*/
create index window_gene_FK on gene
(
   window_id
);

/*==============================================================*/
/* Table: gsw                                                   */
/*==============================================================*/
create table gsw
(
   gsw_id                         int                            not null AUTO_INCREMENT,
   extreme_id                     int                            not null,
   prev_extreme_id                int                            not null,
   window_id                      int                            not null,
   gsw_type                       char(8),
   gsw_distance                   int,
   gsw_distance_crest             int,
   gsw_wave_length                int,
   gsw_amplitude                  double,
   gsw_trough_gc                  double,
   gsw_crest_gc                   double,
   gsw_gradient                   double,
   gsw_cv                         double,
   primary key (gsw_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: window_gsw_FK                                         */
/*==============================================================*/
create index window_gsw_FK on gsw
(
   window_id
);

/*==============================================================*/
/* Index: extreme_gsw_FK                                        */
/*==============================================================*/
create index extreme_gsw_FK on gsw
(
   extreme_id,
   prev_extreme_id
);

/*==============================================================*/
/* Table: indel                                                 */
/*==============================================================*/
create table indel
(
   indel_id                       int                            not null AUTO_INCREMENT,
   prev_indel_id                  int                            not null,
   align_id                       int                            not null,
   indel_start                    int,
   indel_end                      int,
   indel_length                   int,
   indel_seq                      text,
   indel_all_seqs                 text,
   indel_ref_seq                  text,
   left_extand                    int,
   right_extand                   int,
   indel_gc                       double,
   indel_freq                     int,
   indel_occured                  text,
   indel_type                     char(8),
   indel_slippage                 double,
   indel_coding                   double,
   indel_repeats                  double,
   primary key (indel_id, prev_indel_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: align_indel_FK                                        */
/*==============================================================*/
create index align_indel_FK on indel
(
   align_id
);

/*==============================================================*/
/* Table: isw                                                   */
/*==============================================================*/
create table isw
(
   isw_id                         int                            not null AUTO_INCREMENT,
   indel_id                       int                            not null,
   prev_indel_id                  int                            not null,
   align_id                       int                            not null,
   isw_indel_id                   int,
   isw_start                      int,
   isw_end                        int,
   isw_length                     int,
   isw_type                       char(8),
   isw_distance                   int,
   isw_density                    int,
   isw_differences                int,
   isw_pi                         double,
   isw_gc                         double,
   isw_cv                         double,
   isw_d_indel                    double,
   isw_d_noindel                  double,
   isw_d_complex                  double,
   isw_d_bii                      double,
   isw_d_bnn                      double,
   isw_d_indel2                   double,
   isw_d_noindel2                 double,
   isw_d_complex2                 double,
   isw_d_bii2                     double,
   isw_d_bnn2                     double,
   isw_d_indel3                   double,
   isw_d_noindel3                 double,
   isw_d_complex3                 double,
   isw_d_bii3                     double,
   isw_d_bnn3                     double,
   isw_coding                     double,
   isw_repeats                    double,
   isw_cpg_pi                     double,
   primary key (isw_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: indel_isw_FK                                          */
/*==============================================================*/
create index indel_isw_FK on isw
(
   indel_id,
   prev_indel_id
);

/*==============================================================*/
/* Index: align_isw_FK                                          */
/*==============================================================*/
create index align_isw_FK on isw
(
   align_id
);

/*==============================================================*/
/* Table: meta                                                  */
/*==============================================================*/
create table meta
(
   meta_id                        int                            not null AUTO_INCREMENT,
   meta_key                       text,
   meta_value                     text,
   primary key (meta_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Table: ofg                                                   */
/*==============================================================*/
create table ofg
(
   ofg_id                         int                            not null AUTO_INCREMENT,
   window_id                      int                            not null,
   ofg_tag                        char(64),
   ofg_type                       char(64),
   primary key (ofg_id)
)comment = 'Other feature of genome'
ENGINE = MyISAM;

/*==============================================================*/
/* Index: window_ofg_FK                                         */
/*==============================================================*/
create index window_ofg_FK on ofg
(
   window_id
);

/*==============================================================*/
/* Table: ofgsw                                                 */
/*==============================================================*/
create table ofgsw
(
   ofgsw_id                       int                            not null AUTO_INCREMENT,
   ofg_id                         int                            not null,
   window_id                      int                            not null,
   ofgsw_type                     char(8),
   ofgsw_distance                 int,
   ofgsw_cv                       double,
   primary key (ofgsw_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: ofg_ofgsw_FK                                          */
/*==============================================================*/
create index ofg_ofgsw_FK on ofgsw
(
   ofg_id
);

/*==============================================================*/
/* Index: window_ofgsw_FK                                       */
/*==============================================================*/
create index window_ofgsw_FK on ofgsw
(
   window_id
);

/*==============================================================*/
/* Table: segment                                               */
/*==============================================================*/
create table segment
(
   segment_id                     int                            not null AUTO_INCREMENT,
   window_id                      int                            not null,
   segment_type                   char(8),
   segment_gc_mean                double,
   segment_gc_std                 double,
   segment_gc_cv                  double,
   segment_gc_mdcw                double,
   primary key (segment_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: window_segment_FK                                     */
/*==============================================================*/
create index window_segment_FK on segment
(
   window_id
);

/*==============================================================*/
/* Table: sequence                                              */
/*==============================================================*/
create table sequence
(
   seq_id                         int                            not null AUTO_INCREMENT,
   align_id                       int                            not null,
   chr_id                         int                            not null,
   seq_role                       char(8)                        not null,
   seq_position                   int,
   common_name                    char(64),
   chr_name                       char(64),
   chr_start                      int,
   chr_end                        int,
   chr_strand                     char(1),
   seq_length                     int,
   seq_gc                         double,
   seq_runlist                    text,
   seq_seq                        longtext,
   primary key (seq_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: chromosome_sequence_FK                                */
/*==============================================================*/
create index chromosome_sequence_FK on sequence
(
   chr_id
);

/*==============================================================*/
/* Index: align_sequence_FK                                     */
/*==============================================================*/
create index align_sequence_FK on sequence
(
   align_id
);

/*==============================================================*/
/* Table: snp                                                   */
/*==============================================================*/
create table snp
(
   snp_id                         int                            not null AUTO_INCREMENT,
   isw_id                         int,
   align_id                       int                            not null,
   snp_pos                        int,
   target_base                    char(1),
   query_base                     char(1),
   all_bases                      text,
   ref_base                       char(1),
   mutant_to                      char(64),
   snp_freq                       int,
   snp_occured                    text,
   snp_coding                     double,
   snp_repeats                    double,
   snp_cpg                        double,
   primary key (snp_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: align_snp_FK                                          */
/*==============================================================*/
create index align_snp_FK on snp
(
   align_id
);

/*==============================================================*/
/* Index: isw_snp_FK                                            */
/*==============================================================*/
create index isw_snp_FK on snp
(
   isw_id
);

/*==============================================================*/
/* Table: window                                                */
/*==============================================================*/
create table window
(
   window_id                      int                            not null AUTO_INCREMENT,
   align_id                       int                            not null,
   window_start                   int,
   window_end                     int,
   window_length                  int,
   window_runlist                 text,
   window_comparables             int,
   window_identities              int,
   window_differences             int,
   window_indel                   int,
   window_pi                      double,
   window_gc                      double,
   window_coding                  double,
   window_repeats                 double,
   primary key (window_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: align_window_FK                                       */
/*==============================================================*/
create index align_window_FK on window
(
   align_id
);


#----------------------------#
# chromosome
#----------------------------#
CREATE INDEX chr_common_name ON chromosome
(
   common_name, chr_name
);

CREATE INDEX taxon_id_idx ON chromosome
(
   taxon_id
);

CREATE INDEX chr_chr_name ON chromosome
(
    chr_name
);

#----------------------------#
# sequence
#----------------------------#
CREATE INDEX seq_chr_common_name ON sequence
(
    common_name, chr_name, chr_start, chr_end
);

CREATE INDEX seq_chr_name ON sequence
(
    chr_name, chr_start, chr_end
);

CREATE INDEX seq_chr_start ON sequence
(
    chr_start, chr_end
);

CREATE INDEX seq_role ON sequence
(
    seq_role
);

CREATE INDEX seq_position ON sequence
(
    seq_position
);

#----------------------------#
# indel
#----------------------------#
CREATE INDEX indel_type ON indel
(
    indel_type
);

#----------------------------#
# isw
#----------------------------#
CREATE INDEX indel_isw_id_FK ON isw
(
    isw_indel_id
);

CREATE INDEX isw_distance ON isw
(
    isw_distance
);

CREATE INDEX isw_density ON isw
(
    isw_density
);

CREATE INDEX isw_type ON isw
(
    isw_type
);

