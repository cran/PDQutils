## ----'preamble', include=FALSE, warning=FALSE, message=FALSE, cache=FALSE----
library(knitr)

# set the knitr options ... for everyone!
# if you unset this, then vignette build bonks. oh, joy.
#opts_knit$set(progress=TRUE)
opts_knit$set(eval.after='fig.cap')
# for a package vignette, you do want to echo.
# opts_chunk$set(echo=FALSE,warning=FALSE,message=FALSE)
opts_chunk$set(warning=FALSE,message=FALSE)
#opts_chunk$set(results="asis")
opts_chunk$set(cache=TRUE,cache.path="cache/PDQutils")

#opts_chunk$set(fig.path="figure/",dev=c("pdf","cairo_ps"))
opts_chunk$set(fig.path="figure/PDQutils",dev=c("pdf"))
opts_chunk$set(fig.width=4.5,fig.height=4,dpi=300)

# doing this means that png files are made of figures;
# the savings is small, and it looks like shit:
#opts_chunk$set(fig.path="figure/",dev=c("png","pdf","cairo_ps"))
#opts_chunk$set(fig.width=4,fig.height=4)
# for figures? this is sweave-specific?
#opts_knit$set(eps=TRUE)

# this would be for figures:
#opts_chunk$set(out.width='.8\\textwidth')
# for text wrapping:
options(width=64,digits=2)
opts_chunk$set(size="small")
opts_chunk$set(tidy=TRUE,tidy.opts=list(width.cutoff=50,keep.blank.line=TRUE))

compile.time <- Sys.time()

# from the environment

# only recompute if FORCE_RECOMPUTE=True w/out case match.
FORCE_RECOMPUTE <- 
	(toupper(Sys.getenv('FORCE_RECOMPUTE',unset='False')) == "TRUE")

# compiler flags!

# not used yet
LONG.FORM <- FALSE

mc.resolution <- ifelse(LONG.FORM,1000,200)
mc.resolution <- max(mc.resolution,100)

library(PDQutils)

gen_norm <- rnorm
lseq <- function(from,to,length.out) { 
	exp(seq(log(from),log(to),length.out = length.out))
}

## ----'rsnak',eval=TRUE,echo=TRUE------------------------------
rsnak <- function(n,dfs) {
	samples <- Reduce('+',lapply(dfs,function(k) { sqrt(rchisq(n,df=k)/k) }))
}

## ----'testit',eval=TRUE,echo=TRUE,fig.width=4.5,fig.height=4,fig.scap=NA,fig.cap=paste0("A q-q plot of ",n.samp," draws from a sum of Nakagamis distribution is shown against normality.")----
n.samp <- 1e5
dfs <- c(8,15,4000,10000)
set.seed(18181)
# now draw from the distribution 
rvs <- rsnak(n.samp,dfs)
qqnorm(rvs)
qqline(rvs,col='red')

## ----'snakcu',eval=TRUE,echo=TRUE-----------------------------
# for the moment2cumulant function:
library(PDQutils)

# compute the first ord.max raw cumulants of the sum of Nakagami variates
snak_cumulants <- function(dfs,ord.max=10) {
	# first compute the raw moments
	moms <- lapply(dfs,function(nu) { 
		ords <- 1:ord.max
		moms <- 2 ^ (ords/2.0) * exp(lgamma((nu+ords)/2) - lgamma(nu/2))
		# we are dividing the chi by sqrt the d.f.
		moms <- moms / (nu ^ (ords/2.0))
		moms
	})
	# turn moments into cumulants
	cumuls <- lapply(moms,moment2cumulant)
	# sum the cumulants
	tot.cumul <- Reduce('+',cumuls)
	return(tot.cumul)
}

## ----'dpqsnak',eval=TRUE,echo=TRUE----------------------------

library(PDQutils)

dsnak <- function(x,dfs,ord.max=10,...) {
	raw.moment <- cumulant2moment(snak_cumulants(dfs,ord.max))
	retval <- dapx_gca(x,raw.moment,...)
	return(retval)
}
psnak <- function(q,dfs,ord.max=10,...) {
	raw.moment <- cumulant2moment(snak_cumulants(dfs,ord.max))
	retval <- papx_gca(q,raw.moment,...)
	return(retval)
}
qsnak <- function(p,dfs,ord.max=10,...) {
	raw.cumul <- snak_cumulants(dfs,ord.max)
	retval <- qapx_cf(p,raw.cumul,...)
	return(retval)
}

## ----'improvedqq',eval=TRUE,echo=TRUE,fig.width=4.5,fig.height=4,fig.scap=NA,fig.cap=paste0("A q-q plot of ",n.samp,"draws from a sum of Nakagamis distribution is shown against quantiles from the `qsnak' function.")----

qqplot(qsnak(ppoints(n.samp),dfs=dfs),rvs,
	main='Q-Q against qsnak (C-F appx.)')
qqline(rvs, distribution=function(p) qsnak(p,dfs),
	col='red')


