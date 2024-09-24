suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(propagate))
suppressMessages(library(lsa))

srt.immune <- readRDS("E:/Storage/rna_seq_datasets/immune_seurat_norm_annotated.rds") %>%
  subset(subset = cell_label %in% c("lateMacs", "repMacs", "intMacs", "regMacs 1", "regMacs 2", "infMacs"))

xdata <- t(as.matrix(GetAssayData(srt.immune, "data")))

xdata <- as.data.frame(xdata)

ydata <- c(rep(0, length(srt.immune$Condition)))
for (i in 1:length(srt.immune$Condition)) {
  ydata[i] <- 0*(unname(srt.immune$Condition)[i] == "Young") + 1*(unname(srt.immune$Condition)[i] == "Old")
}

ylabel <- unname(srt.immune$cell_label)

data <- data.frame(xdata, y = ydata, label = ylabel)

rm(ylabel, srt.immune)

groups <- unique(data[,'label'])
group_list <- list()

for (i in 1:length(groups)) {
  group_list[[i]] <- data[which(data[,'label'] == groups[i]),!colnames(data) %in% c("y", "label")]
}

for (i in 1:length(groups)) {
  group_list[[i]] <- group_list[[i]][,sapply(seq(ncol(group_list[[i]])), function(ix) {sd(group_list[[i]][,ix])}) != 0]
}

col_common = colnames(xdata)
for (i in 1:length(group_list)){
  col_common = intersect(col_common, colnames(group_list[[i]]))
}

group_list <- lapply(group_list, function(x) x[, col_common])

xdata <- as.data.frame(scale(as.matrix(xdata[, col_common])))

group_list <- lapply(group_list, scale)

group_cor <- lapply(group_list, function (x) bigcor(x, y = NULL, fun = "cor", size = 2000, verbose=FALSE))
group_cor <- lapply(group_cor, function (x) as.data.frame(as.ffdf(x)))
group_cor <- lapply(group_cor, function (x) x %>% mutate(across(everything(), replace_na, 0)))
saveRDS(group_cor, "group_correlation_mat_macs.rds")

itrain <- sort(sample(nrow(xdata), nrow(xdata)*.9))

xdata_train <- makeX(xdata[itrain,], na.impute = TRUE)
xdata_test <- makeX(xdata[-itrain,], na.impute = TRUE)

fit_EN_cv_oldyg <- cv.glmnet(xdata_train, as.factor(ydata[itrain]), family="binomial", nfolds=10, alpha=0.9, penalty.factor=ang_weight, type.measure="mse")

var_selected_idx <- which(fit_EN_cv_oldyg$glmnet.fit$beta[,which(fit_EN_cv_oldyg$cvm == min(fit_EN_cv_oldyg$cvm))] != 0)

var_selected <- colnames(xdata)[var_selected_idx]

print("Evaluate classifier on training data")
assess.glmnet(fit_EN_cv_oldyg, newx = xdata_train, newy = as.factor(ydata[itrain]))

cnf <- confusion.glmnet(fit_EN_cv_oldyg, newx = xdata_train, newy = as.factor(ydata[itrain]))
print(cnf)

cfit_train <- cv.glmnet(xdata_train, as.factor(ydata[itrain]), family = "binomial", nfolds=10, alpha=0.9, penalty.factor=ang_weight, type.measure = "auc", keep = TRUE)
rocs <- roc.glmnet(cfit_train$fit.preval, newy = as.factor(ydata[itrain]))
best <- cfit_train$index["min",]
plot(rocs[[best]], type = "l")
invisible(sapply(rocs, lines, col="grey"))
lines(rocs[[best]], lwd = 2,col = "red")

print("Evaluate classifier on test data")
assess.glmnet(fit_EN_cv_oldyg, newx = xdata_test, newy = as.factor(ydata[-itrain]))

cnf <- confusion.glmnet(fit_EN_cv_oldyg, newx = xdata_test, newy = as.factor(ydata[-itrain]))
print(cnf)

pdf("disTwiner_macs_roc.pdf")
cfit_test <- cv.glmnet(xdata_test, as.factor(ydata[-itrain]), family = "binomial", nfolds=10, alpha=0.9, penalty.factor=ang_weight, type.measure = "auc", keep = TRUE)
rocs <- roc.glmnet(cfit_test$fit.preval, newy = as.factor(ydata[-itrain]))
best <- cfit_test$index["min",]
plot(rocs[[best]], type = "l")
invisible(sapply(rocs, lines, col="grey"))
lines(rocs[[best]], lwd = 2,col = "red")
dev.off()

macs_fit <- cv.glmnet(makeX(xdata, na.impute = TRUE), as.factor(ydata), family="binomial", nfolds=10, alpha=0.9, penalty.factor=ang_weight, type.measure="mse")

genes_selected_idx <- which(macs_fit$glmnet.fit$beta[,which(macs_fit$cvm == min(macs_fit$cvm))] != 0)

genes_selected <- colnames(xdata)[genes_selected_idx]

saveRDS(genes_selected, "macs_cvglmnet_genes_selected.rds")