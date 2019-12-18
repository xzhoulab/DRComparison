

## F-measure for rare cell type identification
FMeasure <- function(true_label, predict_label){
	cm = as.matrix(table(Actual = true_label, Predicted = predict_label)) # create the confusion matrix
	n = sum(cm) # number of instances
	nc = nrow(cm) # number of classes
	diag = diag(cm) # number of correctly classified instances per class 
	rowsums = apply(cm, 1, sum) # number of instances per class
	colsums = apply(cm, 2, sum) # number of predictions per class
	p = rowsums / n # distribution of instances over the actual classes
	q = colsums / n # distribution of instances over the predicted classes

	precision = diag / colsums 
	recall = diag / rowsums 
	f1 = 2 * precision * recall / (precision + recall)
	return(f1)
}# end 
