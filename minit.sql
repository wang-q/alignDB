/*==============================================================*/
/* DBMS name:      MySQL 4.0                                    */
/* Created on:     2009/1/8 21:51:52                            */
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
   primary key (align_id)
)
type = MyISAM;

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
type = MyISAM;

/*==============================================================*/
/* Index: taxon_chromosome_FK                                   */
/*==============================================================*/
create index taxon_chromosome_FK on chromosome
(
   taxon_id
);

/*==============================================================*/
/* Table: indel                                                 */
/*==============================================================*/
create table indel
(
   indel_id                       int                            not null AUTO_INCREMENT,
   foregoing_indel_id             int                            not null,
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
   primary key (indel_id, foregoing_indel_id)
)
type = MyISAM;

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
   foregoing_indel_id             int,
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
type = MyISAM;

/*==============================================================*/
/* Index: indel_isw_FK                                          */
/*==============================================================*/
create index indel_isw_FK on isw
(
   indel_id,
   foregoing_indel_id
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
type = MyISAM;

/*==============================================================*/
/* Table: query                                                 */
/*==============================================================*/
create table query
(
   query_id                       int                            not null AUTO_INCREMENT,
   seq_id                         int,
   query_strand                   char(1),
   query_position                 int,
   primary key (query_id)
)
type = MyISAM;

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
type = MyISAM;

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
type = MyISAM;

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
type = MyISAM;

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
type = MyISAM;

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
type = MyISAM;

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
type = MyISAM;

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
   window_indels                  int,
   window_pi                      double,
   window_target_gc               double,
   window_average_gc              double,
   window_coding                  double,
   window_repeats                 double,
   primary key (window_id)
)
type = MyISAM;

/*==============================================================*/
/* Index: align_window_FK                                       */
/*==============================================================*/
create index align_window_FK on window
(
   align_id
);

