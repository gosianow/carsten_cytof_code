
library(ggplot2)


# -----------------------------------------------------------------------------
# Plot PCA scores from 23_03 vs 29_03
# -----------------------------------------------------------------------------


rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/'
setwd(rwd)

s23 <- read.table("CK_2016-06-23_03/020_pcascores/23_03_princompscore_by_sample.xls", header = TRUE, as.is = TRUE)

s29 <- read.table("CK_2016-06-29_03/020_pcascores/29_03_princompscore_by_sample.xls", header = TRUE, as.is = TRUE)


s <- merge(s23[, c("marker", "avg_score")], s29[, c("marker", "avg_score")], by = "marker", suffixes = c("_23","_29"))


limmin <- min(s[, c("avg_score_23", "avg_score_29")], na.rm = TRUE)
limmax <- max(s[, c("avg_score_23", "avg_score_29")], na.rm = TRUE)

ggp <- ggplot(data = s, aes(x = avg_score_23, y = avg_score_29)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_abline(intercept = 0, slope = 1) +
  coord_cartesian(xlim = c(limmin, limmax), ylim = c(limmin, limmax)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
    axis.title = element_text(size=14, face="bold"))

pdf("avg_score_29vs23.pdf", w=7, h=7)
print(ggp)
dev.off()







