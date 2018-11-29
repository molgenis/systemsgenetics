## ----Install-------------------------------------------------------------
library(devtools)
#install_github("raguirreg/DeconCell")

## ----Libraries-----------------------------------------------------------
library(DeconCell)
library(edgeR)
library(tidyverse)
library(ghibli)


data("count.table")
dCell.exp <- dCell.expProcessing(count.table, trim = TRUE)

## ----Data----------------------------------------------------------------
data("dCell.models")
prediction <- dCell.predict(dCell.exp, dCell.models, res.type = "median")
head(prediction$dCell.prediction)
head(prediction$Evaluation)

## ----evaluation, fig.align = "center", fig.width=7, fig.height=10--------
data("cell.proportions")
library(reshape2)
library(ggplot2)
data("dCell.names")
pData <- data.frame(PearsonCor= diag(cor(cell.proportions, prediction$dCell.prediction)), 
                    CTs = dCell.names[colnames(cell.proportions), "finalName"], 
                    Subpop = dCell.names[colnames(cell.proportions), "broadSubpopulations"])

ggplot(pData, aes(y=PearsonCor , x= CTs, fill=Subpop))+
  geom_bar(stat="identity", alpha=0.8)+
  geom_hline(yintercept = 0.5, alpha=0.5, color="red")+
  coord_flip()+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "bottom", text=element_text(size=11, family="Helvetica"))
  

## ----dCell.run-----------------------------------------------------------
library(DeconRNASeq)
data(multi_tissue)

## remove colums that are not needed.
datasets <- x.data[,2:11]
signatures <- x.signature.filtered.optimal[,2:6]
proportions <- fraction
exp <- datasets



## ----dCell.run models, message = FALSE, warning=FALSE--------------------
set.seed(1121)
sampled.train <- sample(colnames(exp), size = 6, replace = FALSE)
#use the rest of the samples for testing the models
sampled.test <- colnames(exp)[which(colnames(exp) %in% sampled.train == FALSE)]

new.dCell.models <- dCell.run(exp = exp[,sampled.train], 
                              proportions = proportions[sampled.train,], 
                              iterations = 5)



## ----dCell.run eval, fig.align = "center", fig.width=7, fig.height=5-----

test.prediction <- dCell.predict(exp[,sampled.test],
                                 dCell.models= new.dCell.models$deconCell.models.per.CT, 
                                 res.type = "median", custom = TRUE)
# we use custom=TRUE to keep the original names of the proportions.

# reshape the data for plotting 
pData <- reshape2::melt(as.matrix(proportions[sampled.test,]))
pData$Predicted <- reshape2::melt(as.matrix(test.prediction$dCell.prediction))$value

## Function to calculate the Root Mean Square Error
rmse.calculate <- function(x, x.pred){
  sqrt(mean((x - x.pred)^2))
}

tissues <- as.character(unique(pData$Var2))
rmse.per.tissue <- sapply(tissues, function(x){rmse.calculate(pData$value[which(pData$Var2 == x)], pData$Predicted[which(pData$Var2 == x)])})

pData$RMSE <- rmse.per.tissue[as.character(pData$Var2)]
pData$RMSE <- paste0("RMSE= ", format(pData$RMSE,digits= 3))

decon.cell.tissue.plot <- ggplot(pData, aes(x= value, y=Predicted))+
                          geom_point(alpha= 0.9, size=1.5, aes(color= Var2))+
                          facet_grid(facets = ~Var2+RMSE, scales = "free")+
                          geom_smooth(method='lm', lwd=0.5,aes(color= Var2, alpha= 0.5))+
                          ylab("Decon-cell predicted \n tissue proportions")+
                          xlab("Tissue proportions")+
                          scale_color_manual(values = ghibli_palette("KikiMedium")[1:5])+
                          theme_bw()+
                          theme(text = element_text(family = "Helvetica", size = 10), legend.position = "none")


plot(decon.cell.tissue.plot)

