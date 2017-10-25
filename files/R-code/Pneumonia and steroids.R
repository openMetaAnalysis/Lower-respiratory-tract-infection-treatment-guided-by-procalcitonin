library(meta)
library(data.table) # For set names
library(grid)

setwd("../data")

#openMeta-Analysis
analytic.method = "Random effects model (Hartung-Knapp)"
setwd("../data")
filename = file.choose()
newdata<- read.table(filename, header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
meta1 <- metabin(experimental_events, experimental_total, control_events, control_total, data = newdata, sm="RR", hakn = TRUE, method = "Inverse", studlab = paste(author, ", ", year, sep = ""), byvar = Setting, comb.fixed = FALSE) # incr = "TACC", 
summary(meta1)
png(file = "Outcome-mortality.png", width = 800, height = 800, units = "px")
forest(meta1, print.I2.ci=TRUE, print.p=FALSE, print.tau2=FALSE, text.random = analytic.method, text.random.w = analytic.method)
grid.text("Outcome: mortality", 0.5, 0.97, gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()
funnel(meta1)

# Meta-regression
library(metafor)
dat <- escalc(measure='RR' , ai=experimental_events, bi=experimental_total-experimental_events, ci=control_events, di=control_total-control_events, data=newdata)
dat$cofactor <- dat$LRTI.rate/100
xlabel = "LRTI rate"
dat$cofactor <- dat$Adherence/100
xlabel = "Adherence rate to procalcitonin protocol"
res <- rma.uni(yi, vi, mods = ~ dat$cofactor, method = "DL", data=dat, knha = TRUE)
### calculate point sizes by rescaling the standard errors
wi    <- 1/sqrt(dat$vi)
size  <- 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi))
plot(dat$cofactor, exp(dat$yi), pch=19, cex=size, 
     xlab=xlabel, ylab="Relative Risk", main="Metaregression", ylim=c(0.5,2.5), xlim=c(0,1),
     las=1, bty="l", log="y")
abline(h=1, lty="dotted")
abline(res$beta[1], res$beta[2], lty="dotted")
text(par("usr")[2],  10^(par("usr")[4]-1.4*strheight("A"))                     ,cex=1,adj=c(1,0),paste("R2 = ",round(res$R2),"% (QM = ",sprintf(res$QM, fmt='%#.1f'),", p = ",sprintf(res$pval[2], fmt='%#.3f'), ")", sep=""), font=1)
#text(par("usr")[2], 10^(par("usr")[4]-1.2*strheight("A"))                     ,cex=1.2,adj=c(1,0),paste("p (correlation) = ",sprintf(res$pval[2], fmt='%#.3f'), sep=""), font=1)
text(par("usr")[2],  10^(par("usr")[4]-2.4*strheight("A")-0.6*strheight("A"))  ,cex=1,adj=c(1,0),paste("Residual I2 = ",sprintf(res$I2, fmt='%#.1f'),'%', sep=""), font=1)
summary(res)

plot(meta1$data$LRTI.rate, meta1$TE)

cor.test(meta1$TE, meta1$data$LRTI.rate, method= "pearson")

studyweights = meta1$data$.n.e + meta1$data$.n.c
metaregression <- lm(meta1$TE ~ meta1$data$LRTI.rate , data = meta1, weights = studyweights)
summary(metaregression)

