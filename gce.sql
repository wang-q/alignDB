# nature08
create table gce
(
   gce_id                          int        not null AUTO_INCREMENT,
   gce_chr                         text,
   gce_start                       int,
   gce_end                         int,
   gce_runlist                     text,
   gce_tetrad                      text,
   gce_spores                      text,
   gce_type                        text,
   primary key (gce_id)
)
type = MyISAM;

create table hotspot
(
   hotspot_id                          int        not null AUTO_INCREMENT,
   hotspot_chr                         text,
   hotspot_start                       int,
   hotspot_end                         int,
   hotspot_runlist                     text,
   hotspot_type                        text,
   primary key (hotspot_id)
)
type = MyISAM;

