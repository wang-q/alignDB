# alignDB: Analyze the relationship between indels and substitutions in genomes

## STEPS starting from .axt files (two-way)

```
cd ~/Scripts/alignDB
```

1. `perl init/init_alignDB.pl -d S288cvsRM11`

    ```
    ==============================
    Create DB skeleton
    ==============================
    
    ==============================
    Init S288cvsRM11...
    Start at: Fri Jul 17 22:57:41 2015
    
    Use /Users/wangq/Scripts/alignDB/init/../data/taxon.csv to Init table taxon
    Found a row: taxon_id = 9606, genus = Homo, species = sapiens, common_name = Human
    
    Use /Users/wangq/Scripts/alignDB/init/../data/chr_length.csv to Init table chromosome
    Found a row: taxon_id = 3702, chr_name = chrUn, chr_length = 999999999
    
    End at: Fri Jul 17 22:57:41 2015
    Runtime 0 seconds.
    ==============================
    ```

2. `perl init/gen_alignDB.pl -d S288cvsRM11 -t "559292,S288c" -q "285006,RM11" -da data/S288CvsRM11 -lt 5000 --parallel 2`

    ```
    ----Total .axt Files:    0----
    
    ----Total .axt.gz Files:   17----
    
    ===Do task 1 out of 17===
    ===Do task 2 out of 17===
    
    ==============================
    Process data/S288CvsRM11/axtNet/chrI.net.axt.gz...
    ==============================
    
    ==============================
    Process data/S288CvsRM11/axtNet/chrII.net.axt.gz...
    ==============================
    
    Prosess align [1] at S288c.chrI(+):17221-24930
    Prosess align [2] at S288c.chrII(+):9425-29638
    Prosess align [3] at S288c.chrI(+):27070-160233
    
    ...
    
    Prosess align [204] at S288c.chrXVI(+):927529-936894
    
    ==============================
    data/S288CvsRM11/axtNet/chrXVI.net.axt.gz has been processed.
    Runtime 10 seconds.
    ==============================
    
    All files have been processed.
    End at: Fri Jul 17 23:05:29 2015
    Runtime 1 minute and 3 seconds.
    ==============================
    ```

3. `perl init/insert_isw.pl -d S288cvsRM11 --parallel 2`

    ```
    ==============================
    Update isw-indel relationship of S288cvsRM11...
    Start at: Fri Jul 17 23:07:58 2015
    
    Emptying tables...
    ===Do task 1 out of 5===
    ===Do task 2 out of 5===
    Process align [1] at chrI(+):17221-24930
    Process align [51] at chrIV(+):758404-871820
    Process align [2] at chrII(+):9425-29638
    
    ...
    
    Process align [200] at chrXVI(+):572158-720458
    
    End at: Fri Jul 17 23:09:36 2015
    Runtime 1 minute and 38 seconds.
    ==============================
    ```

4. `perl init/insert_gc.pl -d S288cvsRM11 --parallel 2`

    ```
    ==============================
    Update GC tables of S288cvsRM11...
    Start at: Fri Jul 17 23:16:49 2015
    
    Emptying tables...
    ===Do task 1 out of 5===
    ===Do task 2 out of 5===
    Process align [1] at chrI(+):17221-24930
    Process align [51] at chrIV(+):758404-871820
    Process align [2] at chrII(+):9425-29638
    
    ...
    
    Process align [200] at chrXVI(+):572158-720458
        Can't determine right vicinity.
    
    End at: Fri Jul 17 23:20:19 2015
    Runtime 3 minutes and 30 seconds.
    ==============================
    ```

5. `perl init/update_sw_cv.pl -d S288cvsRM11 --parallel 2`

    ```
    ==============================
    Update S288cvsRM11...
    Start at: Fri Jul 17 23:56:05 2015
    
    Table codingsw, ofgsw, isw and gsw altered
    ===Do task 1 out of 5===
    ===Do task 2 out of 5===
    Process align [1] at chrI(+):17221-24930
    Process align [51] at chrIV(+):758404-871820
    Process align [2] at chrII(+):9425-29638
    
    ...
    
    Process align [200] at chrXVI(+):572158-720458
    
    End at: Fri Jul 17 23:57:14 2015
    Runtime 1 minute and 9 seconds.
    ==============================
    ```

6. `perl init/update_feature.pl -d S288cvsRM11 -e yeast_65 --parallel 2`

    ```
    ==============================
    Update annotations of S288cvsRM11...
    Start at: Fri Jul 17 23:59:25 2015
    
    ===Do task 1 out of 5===
    ===Do task 2 out of 5===
    Process align [1] at chrI(+):17221-24930
    Process align [51] at chrIV(+):758404-871820
    Process align [2] at chrII(+):9425-29638
    
    ...
    
    Process align [200] at chrXVI(+):572158-720458
    
    End at: Sat Jul 18 00:03:44 2015
    Runtime 4 minutes and 19 seconds.
    ==============================
    ```

7. `perl init/update_indel_slippage.pl -d S288cvsRM11`

    ```
    ==============================
    Update indel-slippage of S288cvsRM11...
    Start at: Sat Jul 18 00:09:55 2015
    
    Processing align_id 1
    Processing align_id 2
    Processing align_id 3
    
    ...
    
    Processing align_id 204
    
    End at: Sat Jul 18 00:10:14 2015
    Runtime 19 seconds.
    ==============================
    ```

8. `perl stat/common_stat_factory.pl -d S288cvsRM11`

    ```
    ==============================
    Do stat for S288cvsRM11...
    Start at: Sat Jul 18 00:11:40 2015
    
    Sheet "basic" has been generated.
    Sheet "process" has been generated.
    Sheet "summary" has been generated.
    Sheet "d1_pi_gc_cv" has been generated.
    Sheet "d2_pi_gc_cv" has been generated.
    Sheet "d1_comb_pi_gc_cv" has been generated.
    Sheet "d2_comb_pi_gc_cv" has been generated.
    Sheet "group_distance" has been generated.
    Sheet "group_density" has been generated.
    Sheet "d1_comb_coding" has been generated.
    Sheet "d1_comb_non_coding" has been generated.
    Sheet "d2_comb_coding" has been generated.
    Sheet "d2_comb_non_coding" has been generated.
    Sheet "d1_comb_slippage" has been generated.
    Sheet "d1_comb_non_slippage" has been generated.
    Sheet "d2_comb_slippage" has been generated.
    Sheet "d2_comb_non_slippage" has been generated.
    Sheet "dd_group" has been generated.
    Sheet "dd_group_gc" has been generated.
    Sheet "indel_size_group" has been generated.
    Sheet "indel_size_asymmetry" has been generated.
    Sheet "indel_extand_group" has been generated.
    Sheet "indel_extand_asymmetry" has been generated.
    Sheet "indel_position_group" has been generated.
    Sheet "indel_coding_group" has been generated.
    Sheet "indel_repeat_group" has been generated.
    Sheet "indel_slip_group" has been generated.
    Sheet "indel_gc_group" has been generated.
    Sheet "snp_indel_ratio" has been generated.
    Sheet "indel_length" has been generated.
    Sheet "indel_length_100" has been generated.
    Sheet "snp_base_change" has been generated.
    Sheet "distance_snp" has been generated.
    Sheet "density_snp" has been generated.
    
    End at: Sat Jul 18 00:12:08 2015
    Runtime 28 seconds.
    ==============================
    ```
