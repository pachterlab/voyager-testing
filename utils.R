
.checkpoint.counter <- 0
.sync.counter <- 0
.is.primary <- NULL
.root.dir <- NULL
.checkpoint.dir <- NULL
.sync.dir <- NULL
.quiet <- FALSE

init.dirs <- function() {
	args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
	is.primary <- args$primary
	root.dir <- args$dir

	if (is.null(args$q)){
		args$q <- FALSE
	}

	.quiet <<- as.logical(args$q)

	if (is.null(is.primary)){
		args$primary <- FALSE
	}
	if (is.null(root.dir)){
		stop("argument -dir must be specified")
	}
	.is.primary <<- as.logical(is.primary)
	.root.dir <<- root.dir
	checkpoint_name <- paste("checkpoints", as.integer(!.is.primary), sep="_")
	sync_name <- paste("sync", as.integer(!.is.primary), sep="_")
	.checkpoint.dir <<- file.path(root.dir, checkpoint_name)
	.sync.dir <<- file.path(root.dir, sync_name)
	
	if(!dir.exists(.checkpoint.dir)){
		dir.create(.checkpoint.dir, recursive=TRUE)
	}
	if(!dir.exists(.sync.dir)){
		dir.create(.sync.dir, recursive=TRUE)
	}
}

normalize.filename <- function(
	filename,
	ext=NULL,
	sync=FALSE
){
	if (sync){
		filename <- file.path(.sync.dir, filename)
	} else {
		filename <- file.path(.checkpoint.dir, filename)
	}

	if (is.null(ext)){
		ext <- tools::file_ext(fn)
	}
	filename <- tools::file_path_sans_ext(filename)
	
	if (!sync){
		filename <- paste(
			filename, .checkpoint.counter, sep="_"
		)
		.checkpoint.counter <<- .checkpoint.counter + 1
	} else {
		filename <- paste(
			filename, .sync.counter, sep="_"
		)
		.sync.counter <<- .sync.counter + 1
	}

	filename <- paste(filename, ext, sep=".")
	if(!.quiet){
		cat(paste(filename, "\n"))
	}
	filename
}

checkpoint.csv <- function(
	filename,
	dat,
	sync=FALSE,
	row.names=TRUE,
	na='',
	quote=FALSE,
	...
) {
	filename <- normalize.filename(filename, ext="csv", sync=sync)
	dat %>% as.data.frame %>%
		write.csv(filename, row.names=row.names, na=na, quote=quote, ...)
}

checkpoint.txt <- function(
	filename,
	dat,
	sync=FALSE,
	ncolumns=1,
	...
) {
	filename <- normalize.filename(filename, ext="txt", sync=sync)
	dat %>% write(filename, ncolumns=ncolumns, ...)
}

checkpoint.mtx <- function(
	filename,
	dat,
	sync=FALSE,
	...
) {
	filename <- normalize.filename(filename, ext="mtx", sync=sync)
	dat %>% Matrix::writeMM(filename, ...)
}

checkpoint.knn <- function(
	filename,
	listw,
	sync=FALSE,
	...
) {

	glist <- listw$weights
	foo_nb <- listw$neighbours
	n_nbrs <- length(glist[[1]])
	n_nodes <- length(glist)

	rows <- replicate(n_nbrs, seq(n_nodes)) %>% t() %>% as.vector()
	cols <- do.call(rbind, foo_nb) %>% t() %>% as.vector()
	mat_vals <- do.call(cbind, glist) %>% as.vector()
	sp <- Matrix::sparseMatrix(rows, cols, x=mat_vals)
	checkpoint.mtx(filename, sp, sync=sync, ...)
}