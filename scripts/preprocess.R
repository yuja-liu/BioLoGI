preprocess <- function(x, min.expressed.gene = 0, min.expressed.cell = 2, max.expressed.ratio = 1, normalize.by.size.effect = FALSE){

	if (class(x) == 'SummarizedExperiment')
		X <- assays(x)$count
	else if (class(x) == 'matrix')
		X <- x
	else if (is(x, 'sparseMatrix'))
		X <- x
	else
		stop(sprintf('unknown class(x): %s', class(x)))

	M <- ncol(X)
	N <- nrow(X)
	m <- Matrix::colSums(X > 1) >= min.expressed.gene	# cells that have at least min.expressed.gene expreseed genes
	n <- Matrix::rowSums(X > 1) <= max.expressed.ratio * M & Matrix::rowSums(X > 1) >= min.expressed.cell	# genes that are detected in at least min.expressed.cell or at most max.expressed.ratio cells
	if (normalize.by.size.effect){
	  sf <- apply((X[n, m] + 1) / exp(Matrix::rowMeans(log(X[n, m] + 1))), 2, median)
		X <- t(t(X[n, m]) / sf)
	}else
		X <- X[n, m]

	if (class(x) == 'SummarizedExperiment'){
		x <- x[n, m]
		assays(x)$count <- X
	}else if (class(x) == 'matrix'){
		x <- as.matrix(X)
	}else if (is(x, 'sparseMatrix')){
		x <- X
	}
	
	x
} # end of preprocess

