/*==============================================================*/
/* DBMS name:      MySQL 4.0                                    */
/* Created on:     7/29/2012 5:51:14 AM                         */
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
   align_target_gc                double,
   align_average_gc               double,
   align_comparable_runlist       text,
   align_indel_runlist            text,
   align_coding                   double,
   align_repeats                  double,
   align_te                       double,
   align_paralog                  double,
   align_coding_runlist           text,
   align_repeats_runlist          text,
   align_te_runlist               text,
   align_paralog_runlist          text,
   primary key (align_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Table: chromosome                                            */
/*==============================================================*/
create table chromosome
(
   chr_id                         int                            not null AUTO_INCREMENT,
   taxon_id                       int,
   chr_name                       text,
   chr_length                     int,
   primary key (chr_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: taxon_chromosome_FK                                   */
/*==============================================================*/
create index taxon_chromosome_FK on chromosome
(
   taxon_id
);

/*==============================================================*/
/* Table: codingsw                                              */
/*==============================================================*/
create table codingsw
(
   codingsw_id                    int                            not null AUTO_INCREMENT,
   exon_id                        int,
   prev_exon_id                   int,
   window_id                      int,
   codingsw_type                  char(1),
   codingsw_distance              int,
   primary key (codingsw_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: exon_codingsw_FK                                      */
/*==============================================================*/
create index exon_codingsw_FK on codingsw
(
   exon_id,
   prev_exon_id
);

/*==============================================================*/
/* Index: window_codingsw_FK                                    */
/*==============================================================*/
create index window_codingsw_FK on codingsw
(
   window_id
);

/*==============================================================*/
/* Table: exon                                                  */
/*==============================================================*/
create table exon
(
   exon_id                        int                            not null AUTO_INCREMENT,
   prev_exon_id                   int                            not null,
   gene_id                        int,
   window_id                      int,
   exon_stable_id                 char(64)                       not null,
   exon_strand                    char(1),
   exon_phase                     int,
   exon_end_phase                 int,
   exon_frame                     int,
   exon_rank                      int,
   exon_is_full                   int,
   exon_tl_runlist                text,
   exon_seq                       longtext,
   exon_peptide                   longtext,
   primary key (exon_id, prev_exon_id),
   key AK_exon_stable_id (exon_stable_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: gene_exon_FK                                          */
/*==============================================================*/
create index gene_exon_FK on exon
(
   gene_id
);

/*==============================================================*/
/* Index: window_exon_FK                                        */
/*==============================================================*/
create index window_exon_FK on exon
(
   window_id
);

/*==============================================================*/
/* Table: exonsw                                                */
/*==============================================================*/
create table exonsw
(
   exonsw_id                      int                            not null AUTO_INCREMENT,
   window_id                      int,
   exon_id                        int,
   prev_exon_id                   int,
   exonsw_type                    char(1),
   exonsw_distance                int,
   exonsw_density                 int,
   primary key (exonsw_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: exon_exonsw_FK                                        */
/*==============================================================*/
create index exon_exonsw_FK on exonsw
(
   exon_id,
   prev_exon_id
);

/*==============================================================*/
/* Index: window_exonsw_FK                                      */
/*==============================================================*/
create index window_exonsw_FK on exonsw
(
   window_id
);

/*==============================================================*/
/* Table: extreme                                               */
/*==============================================================*/
create table extreme
(
   extreme_id                     int                            not null AUTO_INCREMENT,
   prev_extreme_id                int                            not null,
   window_id                      int,
   extreme_type                   char(1),
   extreme_left_amplitude         double,
   extreme_right_amplitude        double,
   extreme_left_wave_length       double,
   extreme_right_wave_length      double,
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
   window_id                      int,
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
   gene_go                        char(64),
   gene_feature4                  double,
   gene_feature5                  double,
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
/* Table: genesw                                                */
/*==============================================================*/
create table genesw
(
   genesw_id                      int                            not null AUTO_INCREMENT,
   gene_id                        int,
   window_id                      int,
   genesw_type                    char(1),
   genesw_distance                int,
   primary key (genesw_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: gene_genesw_FK                                        */
/*==============================================================*/
create index gene_genesw_FK on genesw
(
   gene_id
);

/*==============================================================*/
/* Index: window_genesw_FK                                      */
/*==============================================================*/
create index window_genesw_FK on genesw
(
   window_id
);

/*==============================================================*/
/* Table: gsw                                                   */
/*==============================================================*/
create table gsw
(
   gsw_id                         int                            not null AUTO_INCREMENT,
   window_id                      int,
   extreme_id                     int,
   prev_extreme_id                int,
   gsw_type                       char(1),
   gsw_distance                   int,
   gsw_density                    int,
   gsw_amplitude                  int,
   primary key (gsw_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: extreme_gsw_FK                                        */
/*==============================================================*/
create index extreme_gsw_FK on gsw
(
   extreme_id,
   prev_extreme_id
);

/*==============================================================*/
/* Index: window_gsw_FK                                         */
/*==============================================================*/
create index window_gsw_FK on gsw
(
   window_id
);

/*==============================================================*/
/* Table: indel                                                 */
/*==============================================================*/
create table indel
(
   indel_id                       int                            not null AUTO_INCREMENT,
   prev_indel_id                  int                            not null,
   align_id                       int,
   indel_start                    int,
   indel_end                      int,
   indel_length                   int,
   indel_seq                      longtext,
   left_extand                    int,
   right_extand                   int,
   indel_gc                       double,
   indel_freq                     int,
   indel_occured                  char(128),
   indel_type                     char(1),
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
   indel_id                       int,
   prev_indel_id                  int,
   isw_indel_id                   int,
   isw_start                      int,
   isw_end                        int,
   isw_length                     int,
   isw_type                       char(1),
   isw_distance                   int,
   isw_density                    int,
   isw_differences                int,
   isw_pi                         double,
   isw_target_gc                  double,
   isw_average_gc                 double,
   isw_d_indel                    double,
   isw_d_noindel                  double,
   isw_d_bii                      double,
   isw_d_bnn                      double,
   isw_d_complex                  double,
   isw_d_bnn2                     double,
   isw_d_complex2                 double,
   isw_d_bii2                     double,
   isw_d_noindel2                 double,
   isw_d_indel2                   double,
   isw_d_indel3                   double,
   isw_d_noindel3                 double,
   isw_d_bii3                     double,
   isw_d_bnn3                     double,
   isw_d_complex3                 double,
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
   window_id                      int,
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
   ofg_id                         int,
   window_id                      int,
   ofgsw_type                     char(1),
   ofgsw_distance                 int,
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
/* Table: query                                                 */
/*==============================================================*/
create table query
(
   query_id                       int                            not null AUTO_INCREMENT,
   seq_id                         int,
   query_strand                   char(1),
   query_position                 int,
   query_tag                      char(64),
   query_type                     char(64),
   primary key (query_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: sequence_query_FK                                     */
/*==============================================================*/
create index sequence_query_FK on query
(
   seq_id
);

/*==============================================================*/
/* Table: reference                                             */
/*==============================================================*/
create table reference
(
   ref_id                         int                            not null AUTO_INCREMENT,
   seq_id                         int,
   ref_raw_seq                    longtext,
   ref_complex_indel              text,
   primary key (ref_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: sequence_reference_FK                                 */
/*==============================================================*/
create index sequence_reference_FK on reference
(
   seq_id
);

/*==============================================================*/
/* Table: segment                                               */
/*==============================================================*/
create table segment
(
   segment_id                     int                            not null AUTO_INCREMENT,
   window_id                      int,
   segment_type                   char(1),
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
   chr_id                         int,
   align_id                       int,
   chr_start                      int,
   chr_end                        int,
   chr_strand                     char(1),
   seq_length                     int,
   seq_seq                        longtext,
   seq_gc                         double,
   seq_runlist                    text,
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
   align_id                       int,
   snp_pos                        int,
   target_base                    char(1),
   all_bases                      char(128),
   ref_base                       char(1),
   mutant_to                      char(128),
   snp_freq                       int,
   snp_occured                    char(128),
   snp_coding                     double,
   snp_repeats                    double,
   snp_cpg                        double,
   primary key (snp_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: isw_snp_FK                                            */
/*==============================================================*/
create index isw_snp_FK on snp
(
   isw_id
);

/*==============================================================*/
/* Index: align_snp_FK                                          */
/*==============================================================*/
create index align_snp_FK on snp
(
   align_id
);

/*==============================================================*/
/* Table: target                                                */
/*==============================================================*/
create table target
(
   target_id                      int                            not null AUTO_INCREMENT,
   seq_id                         int,
   primary key (target_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Index: sequence_target_FK                                    */
/*==============================================================*/
create index sequence_target_FK on target
(
   seq_id
);

/*==============================================================*/
/* Table: taxon                                                 */
/*==============================================================*/
create table taxon
(
   taxon_id                       int                            not null AUTO_INCREMENT,
   genus                          text,
   species                        text,
   sub_species                    text,
   common_name                    text,
   classification                 text,
   primary key (taxon_id)
)
ENGINE = MyISAM;

/*==============================================================*/
/* Table: window                                                */
/*==============================================================*/
create table window
(
   window_id                      int                            not null AUTO_INCREMENT,
   align_id                       int,
   window_start                   int,
   window_end                     int,
   window_length                  int,
   window_runlist                 text,
   window_comparables             int,
   window_identities              int,
   window_differences             int,
   window_indel                   int,
   window_pi                      double,
   window_target_gc               double,
   window_average_gc              double,
   window_coding                  double,
   window_repeats                 double,
   window_ns_indel                double,
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

