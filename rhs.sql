# natgen05
create table window
(
    window_id                       int        not null AUTO_INCREMENT,
    window_chr                      text,
    window_start                    int,
    window_end                      int,
    window_runlist                  text,
    primary key (window_id)
)
ENGINE = MyISAM;

# natgen05
create table hotspot 
(
   hotspot_id                       int        not null AUTO_INCREMENT,
   window_id                        int,
   hotspot_chr                      text,
   hotspot_start                    int,
   hotspot_end                      int,
   hotspot_runlist                  text,
   hotspot_type                     text,
   primary key (hotspot_id)
)
ENGINE = MyISAM;

# science05
create table spot
(
   spot_id                          int        not null AUTO_INCREMENT,
   spot_chr                         text,
   spot_start                       int,
   spot_end                         int,
   spot_runlist                     text,
   spot_type                        text,
   primary key (spot_id)
)
ENGINE = MyISAM;
    
create index window_hotspot_FK on hotspot
(
   window_id
);

