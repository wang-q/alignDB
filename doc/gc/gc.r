
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

