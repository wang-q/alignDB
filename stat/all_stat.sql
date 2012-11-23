#----------------------------------------------------------#
# SUMMARY TABLE                                     
#----------------------------------------------------------#
# Total_length
SELECT t.genus Genus, t.species Speciest, t.sub_species Sub_species, SUM(s.seq_length) Total_length
FROM sequence s, chromosome c, taxon t
WHERE s.chr_id = c.chr_id
AND c.taxon_id = t.taxon_id            
GROUP BY t.taxon_id

# Alignment_length
SELECT CONCAT(vs.target_name, ' vs. ', vs.query_name), SUM(a.align_length) Alignment_length
FROM align a


# No_of_windows
SELECT CONCAT(vs.target_name, ' vs. ', vs.query_name), COUNT(w.isw_id) No_of_windows
FROM indel i, isw w, align a
WHERE a.align_id = i.align_id
AND i.indel_id = w.indel_id

# No_of_SNPs
SELECT CONCAT(vs.target_name, ' vs. ', vs.query_name), SUM(a.differences) No_of_SNPs
FROM align a

# No_of_Indels
SELECT CONCAT(vs.target_name, ' vs. ', vs.query_name), COUNT(i.indel_id) No_of_Indels
FROM indel i, align a
WHERE a.align_id = i.align_id

# Total Pi
SELECT CONCAT(vs.target_name, ' vs. ', vs.query_name), (SUM(a.differences) * 1.0  / SUM(a.comparable_bases)) Pi
FROM align a

# Total Pi 2
SELECT CONCAT(vs.target_name, ' vs. ', vs.query_name), AVG(a.pi) Pi
FROM align a

# Total Pi 3
SELECT original.Name, -0.75 * log2(1- (4.0/3.0)* original.Pi) JC_Pi
FROM
    (SELECT CONCAT(vs.target_name, ' vs. ', vs.query_name) Name, (SUM(a.differences) * 1.0 / SUM(a.comparable_bases) ) Pi
    FROM align a) original

# AVG GC%
SELECT concat(vs.target_name, ' vs. ', vs.query_name),
       sum(a.align_length * a.align_average_gc) /
       sum(a.align_length)
FROM align a

#----------------------------------------------------------#
# INDIVIDUAL STAT                                     
#----------------------------------------------------------#
# base change
SELECT CONCAT(target_base, "->", query_base) base_change, COUNT(snp_id) snp_number
FROM snp
GROUP BY CONCAT(target_base, query_base)

#----------------------------------------------------------#
# SNP -- indel                                        
#----------------------------------------------------------#
# SNP -- indel regression
SELECT a.align_id align_id, a.align_length align_length, a.pi pi,
       a.differences snp_number,
       i.indel_number indel_number,
       a.gaps indel_sum,
       (a.differences / a.align_length * 1000) snp_density,
       (i.Indel_number / a.align_length * 1000) indel_density,
       (a.gaps / a.align_length * 1000) indel_sum_density
FROM (SELECT align_id, COUNT(indel_id) indel_number
      FROM indel i
      GROUP BY align_id) i,
     align a
WHERE i.align_id = a.align_id

# SNP -- indel ratio, group by align_length
SELECT FLOOR(a.align_length / 10000) align_length,
       AVG(a.differences / i.indel_number) `AVG_SNP/Indel`,
       AVG(a.differences / a.align_length * 1000) `AVG_SNP/kb`,
       AVG(i.indel_number / a.align_length * 1000) `AVG_Indel/kb`,
       AVG(a.align_length) `AVG_align_length`,
       AVG(a.pi) `AVG_pi`,
       COUNT(*) COUNT
FROM (SELECT align_id, COUNT(indel_id) indel_number
      FROM indel i
      GROUP BY align_id) i,
     align a
WHERE i.align_id = a.align_id
GROUP BY FLOOR(a.align_length / 10000)

# SNP -- indel ratio, group by align_length
SELECT AVG(a.differences / i.indel_number) `AVG_SNP/Indel`,
       AVG(a.differences / a.align_length * 1000) `AVG_SNP/kb`,
       AVG(i.indel_number / a.align_length * 1000) `AVG_Indel/kb`,
       AVG(a.align_length) `AVG_align_length`,
       AVG(a.pi) `AVG_pi`,
       COUNT(*) COUNT
FROM (SELECT align_id, COUNT(indel_id) indel_number
      FROM indel i
      GROUP BY align_id) i,
     align a
WHERE i.align_id = a.align_id
AND FLOOR(a.align_length / 10000) >= 11

# create temporary table
CREATE TABLE pi_group (p_id INT NOT NULL AUTO_INCREMENT, PRIMARY KEY (p_id))
	ENGINE=InnoDB
	SELECT a.pi `pi`,
	       a.differences `snp`,
	       i.indel_number `indel`,
	       a.gaps `gaps`,
	       a.align_length `align_length`
	FROM (SELECT align_id, COUNT(indel_id) indel_number
	      FROM indel i
	      GROUP BY align_id) i,
	     align a
	WHERE i.align_id = a.align_id
	ORDER BY pi

SELECT CEIL(p.p_id / s.spacing) `Group`,
       AVG(p.pi) `AVG_pi`,
       AVG(p.snp / p.indel) `AVG_SNP/Indel`,
       AVG(p.snp / p.gaps) `AVG_SNP/gaps`,
       COUNT(*) COUNT,
       AVG(p.align_length) `AVG_align_length`,
       AVG(p.snp / p.align_length * 1000) `AVG_SNP/kb`,
       AVG(p.indel / p.align_length * 1000) `AVG_Indel/kb`,
       AVG(p.gaps / p.align_length * 1000) `AVG_gaps/kb`
FROM pi_group p,
     (SELECT IF(MAX(p_id) <= 10000, CEIL(MAX(p_id) / 10), 1000) spacing FROM pi_group) s
GROUP BY CEIL(p.p_id / s.spacing)

DROP TABLE IF EXISTS pi_group

# Processed chr
SELECT c.chr_id, t.genus, t.species, t.sub_species, c.chr_name
FROM target ta, sequence s, chromosome c, taxon t
WHERE ta.seq_id = s.seq_id
AND s.chr_id = c.chr_id
AND c.taxon_id = t.taxon_id
GROUP BY chr_id

# indel_insert count
SELECT indel_insert, count(indel_insert) FROM indel i
GROUP BY indel_insert

#----------------------------------------------------------#
# CHROMOSOME                                          
#----------------------------------------------------------#
# indel_length distribution
SELECT i.indel_length, COUNT(i.indel_length) indel_number, SUM(i.indel_length) indel_sum
FROM indel i, (SELECT distinct i.indel_id indel_id
               FROM indel i, target t, query q, sequence s1, sequence s2
               WHERE i.align_id = t.align_id
               AND i.align_id = q.align_id
               AND t.seq_id = s1.seq_id
               AND q.seq_id = s2.seq_id
               AND s1.chr_id IN (102, 103)
               AND s2.chr_id IN (235, 234)
              ) chr
WHERE i.indel_id = chr.indel_id
GROUP BY i.indel_length

# stand-alone indel length distribution
SELECT i.indel_length, COUNT(i.indel_length) indel_number, SUM(i.indel_length) indel_sum
FROM indel i, (SELECT distinct i.indel_id indel_id
               FROM indel i, target t, query q, sequence s1, sequence s2
               WHERE i.align_id = t.align_id
               AND i.align_id = q.align_id
               AND t.seq_id = s1.seq_id
               AND q.seq_id = s2.seq_id
               AND s1.chr_id IN (102, 103)
               AND s2.chr_id IN (235, 234)
              ) chr
WHERE i.indel_id = chr.indel_id
AND i.left_extand >= 100
AND i.right_extand >= 100
GROUP BY i.indel_length

# base change
SELECT CONCAT(s.target_base, "->", s.query_base) base_change, COUNT(s.snp_id) SNP_number
FROM snp s, (SELECT distinct s.snp_id snp_id
               FROM snp s, target t, query q, sequence seq1, sequence seq2
               WHERE s.align_id = t.align_id
               AND s.align_id = q.align_id
               AND t.seq_id = seq1.seq_id
               AND q.seq_id = seq2.seq_id
               AND seq1.chr_id IN (102, 103)
               AND seq2.chr_id IN (235, 234)
              ) chr
WHERE s.snp_id = chr.snp_id
GROUP BY CONCAT(s.target_base, s.query_base)

# isw average in chromosome
SELECT AVG(isw_pi), COUNT(isw_pi), STD(isw_pi) 
FROM isw, (SELECT distinct i.indel_id indel_id
           FROM indel i, target t, query q, sequence s1, sequence s2
           WHERE i.align_id = t.align_id
           AND i.align_id = q.align_id
           AND t.seq_id = s1.seq_id
           AND q.seq_id = s2.seq_id
           AND s1.chr_id IN (102, 103)
           AND s2.chr_id IN (235, 234)
          ) indel
WHERE isw.indel_id = indel.indel_id

# density effect in chromosome
SELECT isw_density density,
       AVG(isw_pi) AVG_pi,
       COUNT(isw_pi) COUNT,
       STD(isw_pi) STD_pi
FROM isw, (SELECT distinct i.indel_id indel_id
           FROM indel i, target t, query q, sequence s1, sequence s2
           WHERE i.align_id = t.align_id
           AND i.align_id = q.align_id
           AND t.seq_id = s1.seq_id
           AND q.seq_id = s2.seq_id
           AND s1.chr_id IN (102, 103)
           AND s2.chr_id IN (235, 234)
          ) indel
WHERE isw.indel_id = indel.indel_id
GROUP BY isw_density

# distance effect in chromosome
SELECT isw_distance distance, 
       AVG(isw_pi) AVG_pi, 
       COUNT(isw_pi) COUNT, 
       STD(isw_pi) STD_pi
FROM isw, (SELECT distinct i.indel_id indel_id
           FROM indel i, target t, query q, sequence s1, sequence s2
           WHERE i.align_id = t.align_id
           AND i.align_id = q.align_id
           AND t.seq_id = s1.seq_id
           AND q.seq_id = s2.seq_id
           AND s1.chr_id IN (102, 103)
           AND s2.chr_id IN (235, 234)
          ) indel
WHERE isw.indel_id = indel.indel_id
GROUP BY isw_distance


#----------------------------------------------------------#
# GC% & SNP                                            
#----------------------------------------------------------#
# density effect
SELECT isw_density density, 
       AVG(isw_average_gc) AVG_gc, 
       COUNT(isw_average_gc) COUNT, 
       STD(isw_average_gc) STD_gc
FROM isw i
GROUP BY isw_density

# base change
SELECT AVG(i.isw_distance) AVG_density, COUNT(s.snp_id) snp_number
FROM snp s, isw i
WHERE s.isw_id = i.isw_id
AND CONCAT(target_base, query_base) IN ("AC", "CA")
AND isw_density IN (-1, 0)

# density effect in chromosome
SELECT isw_density density,
       AVG(isw_pi) AVG_pi,
       COUNT(isw_pi) COUNT,
       STD(isw_pi) STD_pi
FROM isw, (SELECT distinct i.indel_id indel_id
           FROM indel i
           WHERE i.indel_length = 2
           AND i.indel_seq IN ("TA", "TA")
          ) indel
WHERE isw.indel_id = indel.indel_id
GROUP BY isw_density

# density effect
SELECT isw_density density, 
       AVG((isw_target_dG + isw_query_dG) / 2) AVG_dG, 
       COUNT((isw_target_dG + isw_query_dG) / 2) COUNT, 
       STD((isw_target_dG + isw_query_dG) / 2) STD_dG
FROM isw i
WHERE isw_length = 100
GROUP BY isw_density

# density effect
SELECT isw_density density, 
       AVG(isw_average_gc) AVG_gc, 
       COUNT(isw_average_gc) COUNT, 
       STD(isw_average_gc) STD_gc
FROM isw i
GROUP BY isw_density

SELECT isw_distance d1, 
       isw_density d2, 
       isw_pi pi, 
       isw_average_gc gc
FROM isw i

SELECT a.align_id align_id, 
       a.align_length align_length,
       count(i.indel_id) indel,
       a.gaps gaps,
       a.pi align_pi, 
       a.align_average_gc align_gc
FROM align a, indel i
WHERE a.align_id = i.align_id
GROUP BY a.align_id

#----------------------------------------------------------#
# GC wave                                        
#----------------------------------------------------------#
# crest_trough length & indel
SELECT e.extreme_type TYPE,
       COUNT(e.window_id) COUNT, 
       AVG(w.window_length) AVG_length,
       SUM(w.window_length) SUM_length,
       SUM(w.window_length) / a.sum_length * 100 `LENGTH%`,
       SUM(w.window_indel) indel,
       SUM(w.window_indel) / i.sum_indel * 100 `INDEL%`,
       SUM(w.window_indel) / SUM(w.window_length) * 100 FACTOR
FROM extreme e, window w, (SELECT SUM(comparable_bases) sum_length
                           FROM align) a,
                          (SELECT COUNT(indel_id) sum_indel
                           FROM indel) i 
WHERE w.window_id = e.window_id
GROUP BY e.extreme_type

# gsw length & indel
SELECT g.gsw_type TYPE, 
       COUNT(g.window_id) COUNT, 
       AVG(w.window_length) AVG_length, 
       SUM(w.window_length) SUM_length, 
       SUM(w.window_length) / a.sum_length * 100 `LENGTH%`,
       SUM(w.window_indel) indel,
       SUM(w.window_indel) / i.sum_indel * 100 `INDEL%`,
       SUM(w.window_indel) / SUM(w.window_length) * 100 FACTOR
FROM gsw g, window w, (SELECT SUM(comparable_bases) sum_length
                       FROM align) a,
                      (SELECT COUNT(indel_id) sum_indel
                       FROM indel) i 
WHERE w.window_id = g.window_id
GROUP BY g.gsw_type

# comparable length & indel
SELECT 'All' TYPE, 
       COUNT(DISTINCT a.align_id) COUNT,
       AVG(a.comparable_bases) AVG_length, 
       SUM(a.comparable_bases) SUM_length,
       SUM(a.comparable_bases) / a2.sum_length * 100 `LENGTH%`,
       SUM(i.indel) indel,
       SUM(i.indel) / i2.sum_indel * 100 `INDEL%`,
       SUM(w.window_indel) / SUM(w.window_length) * 100 FACTOR
FROM align a, (SELECT a.align_id,
                      COUNT(i.indel_id) indel
               FROM align a, indel i
               WHERE a.align_id = i.align_id
               GROUP BY a.align_id) i, 
              (SELECT SUM(comparable_bases) sum_length
               FROM align) a2,
              (SELECT COUNT(indel_id) sum_indel
               FROM indel) i2 
WHERE i.align_id = a.align_id

# distance effect
SELECT g.gsw_distance distance, 
       AVG(w.window_indel) AVG_indel, 
       COUNT(w.window_indel) COUNT, 
       STD(w.window_indel) STD_indel
FROM gsw g, window w
WHERE g.window_id = w.window_id
GROUP BY gsw_distance

# distance_count
SELECT gsw_distance distance, COUNT(gsw_id) COUNT
FROM gsw g
GROUP BY gsw_distance

# combined distance effect
SELECT AVG(gsw_distance) AVG_distance,
       AVG(w.window_indel) AVG_indel,
       COUNT(w.window_indel) COUNT,
       STD(w.window_indel) STD_indel
FROM gsw g, window w
WHERE g.window_id = w.window_id
AND gsw_distance IN (1, 2, 3)

# density effect
SELECT g.gsw_density density, 
       AVG(w.window_indel) AVG_indel, 
       COUNT(w.window_indel) COUNT, 
       STD(w.window_indel) STD_indel
FROM gsw g, window w
WHERE g.window_id = w.window_id
GROUP BY gsw_density

# density_count
SELECT gsw_density density, COUNT(gsw_id) COUNT
FROM gsw g
GROUP BY gsw_density

# combined density effect
SELECT AVG(gsw_density) AVG_density,
       AVG(w.window_indel) AVG_indel,
       COUNT(w.window_indel) COUNT,
       STD(w.window_indel) STD_indel
FROM gsw g, window w
WHERE g.window_id = w.window_id
AND gsw_density IN (1, 2, 3)

# dd_group
SELECT CONCAT(gsw_type, gsw_distance) distance,
       AVG(w.window_indel) AVG_indel,
       COUNT(w.window_indel) COUNT,
       STD(w.window_indel) STD_indel
FROM gsw g, window w
WHERE g.gsw_type = ?
AND g.window_id = w.window_id
AND g.gsw_density BETWEEN ? AND ?
AND g.gsw_distance <= ?
GROUP BY CONCAT(g.gsw_type, g.gsw_distance)
ORDER BY g.gsw_distance

# processecd chromosome
SELECT c.chr_id, c.chr_name
FROM target t, sequence s, chromosome c
WHERE s.chr_id = c.chr_id
AND t.seq_id = s.seq_id
GROUP BY c.chr_id

#----------------------------------------------------------#
# coding & non_codig                                
#----------------------------------------------------------#
# indel number coding--non_coding
SELECT 'coding', coding.COUNT, 'non_coding', non_coding.COUNT
FROM   (SELECT COUNT(i.indel_id) COUNT
        FROM indel i
        WHERE 1 = 1
        AND (i.left_extand >= 100 OR i.right_extand >= 100)
        AND i.indel_coding BETWEEN 1 AND 1) coding,
       (SELECT COUNT(i.indel_id) COUNT
        FROM indel i
        WHERE 1 = 1
        AND (i.left_extand >= 100 OR i.right_extand >= 100)
        AND i.indel_coding BETWEEN 0 AND 0) non_coding

# isw length coding--non_coding
SELECT 'coding', coding.SUM, 
       'non_coding', non_coding.SUM,
       'isw', isw.SUM
FROM   (SELECT SUM(i.isw_length) SUM
        FROM isw i
        WHERE 1 = 1
        AND i.isw_distance != -1
        AND i.isw_coding BETWEEN 1 AND 1) coding,
       (SELECT SUM(i.isw_length) SUM
        FROM isw i
        WHERE 1 = 1
        AND i.isw_distance != -1
        AND 1.isw_coding BETWEEN 0 AND 0) non_coding,
       (SELECT SUM(i.isw_length) SUM
        FROM isw i
        WHERE i.isw_distance != -1) isw

# snp number coding--non_coding 
SELECT 'coding', coding.COUNT, 'non_coding', non_coding.COUNT
FROM   (SELECT COUNT(*) COUNT
        FROM isw i, snp s
        WHERE i.isw_id = s.isw_id
        AND i.isw_coding BETWEEN 1 AND 1) coding,
       (SELECT COUNT(*) COUNT
        FROM isw i, snp s
        WHERE i.isw_id = s.isw_id
        AND i.isw_coding BETWEEN 0 AND 0) non_coding

#----------------------------------------------------------#
# segment                                
#----------------------------------------------------------#
# extreme numbers in segments

# add column
ALTER TABLE `segment` ADD COLUMN `segment_feature4` DOUBLE;

# extreme numbers in segments
UPDATE segment,
       (SELECT s.id id, COUNT(*) count
        FROM 
               (SELECT s.segment_id id, 
                       w.window_start start, 
                       w.window_end end
                FROM   segment s, window w
                WHERE  s.window_id = w.window_id
                AND    w.align_id = 1) s,
               (SELECT e.extreme_id id, 
                       w.window_start start, 
                       w.window_end end
                FROM   extreme e, window w
                WHERE  e.window_id = w.window_id
                AND    w.align_id = 1) e
        WHERE e.start BETWEEN s.start AND s.end
        OR e.end BETWEEN s.start AND s.end
        GROUP BY s.id) seg
SET segment.segment_feature4 = seg.count
WHERE segment.segment_id = seg.id

#----------------------------------------------------------#
# JC_correction.sql                                        
#----------------------------------------------------------#
# real
create table isw_b (select * from isw)

update isw
set isw_pi =  -0.75 * log2(1- (4.0/3.0)* isw_pi)
where isw_id between 1 and 1000000

update align
set pi =  -0.75 * log2(1- (4.0/3.0)* pi)

# restore
update align
set pi =  differences / comparable_bases

update isw, isw_b
set isw.isw_pi =  isw_b.isw_pi
where isw.isw_id between 1 and 1000000
and isw.isw_id = isw_b.isw_id

# test
SELECT avg(isw_pi) FROM isw

SELECT -0.75 * log2(1- (4.0/3.0)* avg(isw_pi)) FROM isw

create table isw_b (select * from isw)

update isw_b
set isw_pi =  -0.75 * log2(1- (4.0/3.0)* isw_pi)

SELECT avg(isw_pi) FROM isw_b i

#----------------------------------------------------------#
# Fix errors                                         
#----------------------------------------------------------#
# isw_distance error
UPDATE isw
SET isw_distance = 1
WHERE isw_density >= 1
AND isw_distance = 0
AND isw_length >= 100


# indel_insert deletion
UPDATE indel
SET indel_insert = 'D'
WHERE indel_insert = 'S288C'
# insertion
UPDATE indel
SET indel_insert = 'I'
WHERE indel_insert = 'RM11'
# check
SELECT indel_insert, COUNT(*)
FROM indel i
GROUP BY indel_insert;

# pick up align
SELECT a.align_id, a.align_length, s.segment_gc_std
FROM align a, segment s, window w
where s.segment_type = '1'
AND s.window_id = w.window_id
AND w.align_id = a.align_id
AND a.align_length < 50000
order by s.segment_gc_std DESC

#----------------------------------------------------------#
# 3-sequences test                                      
#----------------------------------------------------------#
SELECT isw_distance, avg(isw_d_indel), avg(isw_d_noindel)
FROM isw i
where isw_density > 10
group by isw_distance

SELECT indel.indel_occured, isw_distance, (indel.indel_occured = snp.snp_occured), count(*)
FROM indel, isw, snp
where indel.indel_id = isw.indel_id
AND isw.isw_id = snp.isw_id
AND isw.isw_distance between 0 and 10
group by indel.indel_occured,
isw_distance,
(indel.indel_occured = snp.snp_occured)
order by indel.indel_occured desc,
isw_distance,
(indel.indel_occured = snp.snp_occured) desc

select snp_pos, ref_base, target_base, query_base, isw_tyep, isw_distance
from snp, isw
where isw.isw_id = snp.isw_id
and align_id = 1

#----------------------------------------------------------#
# regression                                         
#----------------------------------------------------------#
# dd_group
SELECT isw_distance,
       isw_density,
       AVG (isw_pi) avg_pi
FROM isw i
WHERE isw_distance >= 0 AND
      isw_distance <= 30 AND
      isw_density <= 60
GROUP BY isw_distance,
         isw_density

#----------------------------------------------------------#
# coding & CV                                         
#----------------------------------------------------------#
SELECT w.window_coding `coding`,
       AVG(s.segment_gc_cv) `avg_cv`,
       STD(s.segment_gc_cv) `std_cv`,
       count(*) `count`
FROM segment s,
     window w
WHERE s.window_id = w.window_id AND
      s.segment_type = '3' AND
      w.window_coding IN (0, 1)
GROUP BY coding;


