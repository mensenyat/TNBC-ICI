load("~/Simplified Environment to Test your Samples.RData")

YourData = read.table("Your_Data.txt", sep = "\t", header = T, row.names = 1) # Datatable with genes in the columns and samples in the rows

# Train the RF algorithm using your available genes #
RFgenes = RFgenes[RFgenes %in% colnames(YourData)]
data_train = data_train[which(colnames(data_train) %in% RFgenes)]
data_train$Label = metaTrain$Response

RF <- randomForest(factor(Label) ~ ., data = data_train, importance = TRUE, type = "classification", ntree = 10000)

# Test your samples #
YourData = YourData[,which(colnames(YourData) %in% RFgenes)]
YourData = YourData[,order(match(colnames(YourData), colnames(data_train)))]

YourResult = predict(RF, YourData, type = "prob")

