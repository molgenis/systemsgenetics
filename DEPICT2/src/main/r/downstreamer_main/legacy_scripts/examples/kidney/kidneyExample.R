setwd("C:\\Users\\patri\\Dropbox\\UMCG\\GrantApplications\\2020 Veni")

load("calcium_nephr_ROC_objects.RData")

library(pROC)
pdf(paste0(hpo_name,".pdf"), width = 10, height = 10, version = 1.2)
par(pty="s", bty = "n", cex = 2.2, cex.axis = 1, las = 1, cex.lab = 1.2, xpd = NA)
plot.roc(KN.roc, col = "green", main = hpo_name, mgp=c(2.6, 0.7, 0), lwd = 3)
lines.roc(GN.roc, col = "red", lwd = 3)
legend(x = "bottomright", legend = c("KN", "GN"), lwd=3, bty="n")
dev.off()


svg(paste0(hpo_name,".svg"), width = 10, height = 10)
par(pty="s", bty = "n", cex = 2.2, cex.axis = 1, las = 1, cex.lab = 1.2, xpd = NA)
plot.roc(KN.roc, col = "green", main = hpo_name, mgp=c(2.6, 0.7, 0), lwd = 3)
lines.roc(GN.roc, col = "red", lwd = 3)
legend(x = "bottomright", legend = c("KN", "GN"), lwd=3, bty="n")
dev.off()


pdf(paste0(hpo_name,"_5.pdf"), width = 10, height = 10, useDingbats = F, compress = F)
par(mar = c(10,10,10,10))
plot.new()
plot.window(xlim = c(1,0), ylim = c(0,1))
axis(1)
axis(2)

lines(x=KN.roc[["specificities"]], y = KN.roc[["sensitivities"]], cex = 0.5, pch = 20)
lines(x=GN.roc[["specificities"]], y = GN.roc[["sensitivities"]], cex = 0.5, pch = 20)
dev.off()

