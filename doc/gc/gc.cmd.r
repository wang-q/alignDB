# prepare database and annotations
if (FALSE) {
cd /home/wangq/Scripts/alignDB/ 

perl init/init_alignDB.pl -d S288Cvsself_gc 
perl init/gen_alignDB_genome.pl -d S288Cvsself_gc -t '4932,S288C' -dir ~/data/alignment/yeast65/S288C --length 2000000 --parallel 8

perl init/insert_gc.pl -d S288Cvsself_gc --one_level --parallel 8 --batch 1
perl gene/insert_gene.pl -d S288Cvsself_gc -e yeast_65 --parallel 8 --batch 1

perl init/update_sw_cv.pl -d S288Cvsself_gc --batch 1 --parallel 8 
perl init/update_feature.pl -d S288Cvsself_gc -e yeast_65 --batch 1 --parallel 8

perl gene/update_gene_yeast_ess.pl -d S288Cvsself_gc

perl util/query_sql.pl -s 114.212.202.159 -d S288Cvsself_gc -f doc/gc/query_gc_one.txt -o doc/gc/S288Cvsself_gc.csv
}

# in R
mydata <- read.csv("d:/wq/Scripts/alignDB/doc/gc/S288Cvsself_gc.1.csv")

head(mydata)
summary(mydata)

library(ggplot2)
ggplot(mydata, aes(x = segment_gc_mean)) + geom_histogram(colour = "black", fill = "white") + geom_vline(aes(xintercept = mean(segment_gc_mean, 
    na.rm = T)), color = "red", linetype = "dashed", size = 2) + labs(x = "segment_gc_mean")

ggplot(mydata, aes(x = segment_gc_cv)) + geom_histogram(colour = "black", fill = "white") + geom_vline(aes(xintercept = mean(segment_gc_cv, 
    na.rm = T)), color = "red", linetype = "dashed", size = 2) + labs(x = "segment_gc_cv")

library(changepoint)
y = c(rnorm(100, 0, 1), rnorm(100, 5, 1))
ansmean = cpt.mean(y)
plot(ansmean, cpt.col = "blue")
print(ansmean)


align_plot <- function(data, string) {
    library(changepoint)
    filename <- paste(c("d:/", string, ".png"), collapse = "")
    png(filename, width = 9, height = 3, units = "in", res = 300, pointsize = 9, )
    par(mfrow = c(2, 1), mar = c(2, 4, 2, 2) + 0.1)
    
    data.cv_pelt = cpt.mean(data$segment_gc_cv, method = "PELT", penalty = "Manual", pen.value = 0.1)
    plot(data.cv_pelt, ylab = "window_cv", xlab = "x 100 bp", cpt.width = 3)
    title(main = string)
    data.gc_pelt = cpt.mean(data$segment_gc_mean, method = "PELT", penalty = "Manual", pen.value = 0.1)
    plot(data.gc_pelt, ylab = "window_gc", xlab = "x 100 bp", cpt.width = 3)
    # plot(data$window_coding, ylab = 'window_coding', xlab = 'x 100 bp', )
    dev.off()
    
    list(string = string, nrow = nrow(data), n_cv = ncpts(data.cv_pelt), n_gc = ncpts(data.gc_pelt), cv_cpt = cpts(data.cv_pelt), cv_param = param.est(data.cv_pelt)$mean)
}

report_cv <- function(mydata, result) {
    filename <- paste(c("d:/", result$string, ".txt"), collapse = "")
    cat(file = filename, append = TRUE, paste(c(result$string, " has ", result$n_cv, " changepoints for CV and ", result$n_gc, " changepoints for GC", 
        "\n"), collapse = ""))
    
    cpt_with_ends <- result$cv_cpt
    if (cpt_with_ends[1] != 1) 
        cpt_with_ends = c(1, cpt_with_ends)
    if (tail(cpt_with_ends, n = 1) != result$nrow) 
        cpt_with_ends = c(cpt_with_ends, result$nrow)
    cat(file = filename, append = TRUE, paste(c("Points: ", result$cv_cpt, "\n\n"), collapse = " "))
    
    low_cv_idx <- which(result$cv_param < 0.1)
    for (i in low_cv_idx) {
        start = cpt_with_ends[i]
        end = cpt_with_ends[i + 1]
        cat(file = filename, append = TRUE, paste(c("Region: ", start, "-", end, "\n"), collapse = ""))
        
        seg_start = unlist(strsplit(as.character(mydata$window_runlist[start]), "-"))[1]
        seg_end = unlist(strsplit(as.character(mydata$window_runlist[end]), "-"))[2]
        cat(file = filename, append = TRUE, paste(c("Position: ", seg_start, "-", seg_end, "\n"), collapse = ""))
        
        cat(file = filename, append = TRUE, paste(c("CV: ", round(result$cv_param[i], digits = 4), "\n\n"), collapse = ""))
    }
}
 
myrun <- function(mydata, string, id) {    
    subdata = subset(mydata, mydata$align_id == id)
    result <- align_plot(subdata, string)
    report_cv(subdata, result)    
}


myrun(mydata, "chrI", 1) 
myrun(mydata, "chrII", 6)
myrun(mydata, "chrIII", 3)
myrun(mydata, "chrIV", 11)
myrun(mydata, "chrV", 5)
myrun(mydata, "chrVI", 2)
myrun(mydata, "chrVII", 8)
myrun(mydata, "chrVIII", 7)
myrun(mydata, "chrIX", 4)
myrun(mydata, "chrX", 9)
myrun(mydata, "chrXI", 10)
myrun(mydata, "chrXII", 12)
myrun(mydata, "chrXIII", 13)
myrun(mydata, "chrXIV", 14)
myrun(mydata, "chrXV", 15)
myrun(mydata, "chrXVI", 16)

    445901-450200
if(FALSE) {
    select g.*, w.window_runlist 
    from gene g, window w 
    where 1 = 1 
    AND g.window_id = w.window_id 
    AND w.align_id = 11 
    AND (w.window_start between 1436501 and 1449900 
    OR w.window_end between 1436501 and 1449900)  

    select concat(g.gene_stable_id, "(", g.gene_external_name, "; ", g.gene_feature4, ")")
    from gene g, window w 
    where 1 = 1 
    AND g.window_id = w.window_id 
    AND w.align_id = 16
    AND (w.window_start between 445901 and 450200 
    OR w.window_end between 445901 and 450200
    OR (w.window_start < 445901 and w.window_end > 450200) ) 
}
