suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(e1071)) 


srt.immune <- readRDS("E:/Storage/rna_seq_datasets/immune_seurat_norm_annotated.rds") %>%
  subset(subset = cell_label %in% c("lateMacs", "repMacs", "intMacs", "regMacs 1", "regMacs 2", "infMacs"))

xdata <- t(as.matrix(GetAssayData(srt.immune, slot = "data", assay = "RNA")))

xdata <- xdata[,sapply(seq(ncol(xdata)), function(ix) {sd(xdata[,ix])}) != 0]

ydata <- c(rep(0, length(srt.immune$Condition)))
for (i in 1:length(srt.immune$Condition)) {
  ydata[i] <- 0*(unname(srt.immune$Condition)[i] == "Young") + 1*(unname(srt.immune$Condition)[i] == "Old")
}

dataset <- data.frame(cbind(xdata, ydata))
dataset$ydata <- factor(dataset$ydata, levels = c(0, 1)) 

rm(srt.immune)

set.seed(123)
itrain <- sort(sample(nrow(xdata), nrow(xdata)*.9))

training_set = dataset[itrain,]
test_set = dataset[-itrain,]

training_set[,-16836] = scale(training_set[,-16836]) 
test_set[,-16836] = scale(test_set[,-16836]) 

classifier = svm(formula = ydata ~ ., 
                 data = training_set, 
                 type = 'C-classification', 
                 kernel = 'linear') 

y_pred = predict(classifier, newdata = test_set[-16836]) 

cm = table(test_set[, 16836], y_pred) 
saveRDS(cm, "SVM_confusionmatrix.rds")