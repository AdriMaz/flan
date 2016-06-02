\name{number of mutant cells distribution}
\alias{dflan}
\alias{pflan}
\alias{qflan}
\alias{rflan}
\encoding{UTF-8}


\title{The distribution of mutant cell counts}

\description{
  Density, distribution function, quantile function and random generation for mutant cell counts.
}

\usage{
  dflan(m,mutations=1,fitness=1,death=0,model=c("LD","H"))
  pflan(m,mutations=1,fitness=1,death=0,model=c("LD","H"),lower.tail=TRUE)
  qflan(p,mutations=1,fitness=1,death=0,model=c("LD","H"),lower.tail=TRUE)
  rflan(n,mutations=1,mutprob=NULL,fitness=1,death=0,
        dist=list(name="lnorm",meanlog=-0.3795851,sdlog=0.3016223),
        mfn=1e9,cvfn=0) 
}

\arguments{
  \item{m}{ vector of quantiles. }
  \item{p}{ vector of probabilities. }
  \item{n}{ number of observations. If \code{length(n) > 1}, the length is taken to be the number required. }
  \item{mutations}{mean number of mutations: positive integer. Ignored if \code{mutprob} is non-empty. }
  \item{mutprob}{ mutation probability: numeric between 0 and 1. By default empty. See details. }
  \item{fitness}{ fitness parameter: numeric positive. }
  \item{death}{ death probability: numeric between 0 and 0.5. } 
  \item{dist}{ lifetime distribution for mutant cells. See Details. }
  \item{model}{ statistical lifetime model. Must be one of "LD" (default) for Luria-Delbrück model (exponential lifetimes), or "H" for Haldane model (constant lifetimes). }
  \item{mfn}{ mean final number of cells: numeric positive. }
  \item{cvfn}{ coefficient of variation of final numbers of cells: numeric, default 0. If non-zero and if \code{mutprob} is empty, compute the sample with \code{mutprob} = \code{mutations}/\code{mfn}}
  \item{lower.tail}{ logical; if TRUE (default), probabilities are P[X <= x] otherwise, P[X > m] }
}


\details{
  The argument \code{dist} is a list beginning with the distribution name followed by its parameters, and must be one the 4 following distributions: "dirac"(location), "exp"(rate), "lnorm"(meanlog, sdlog), "gamma"(shape, scale). 
    
  If \code{cvfn} is positive, the final numbers of cells are generated with the log-normal distribution with mean \code{mfn} and coefficient of variation \code{cvfn}.
  
  If \code{mutprob} is non-empty if \code{cvfn} is zero, the sample is computed with \code{mutations} as the product of \code{mutprob} by \code{mfn}. If \code{cvfn} is non-zero, the sample is computed with \code{mutations} as 
  
  \deqn{\alpha = \mathcal{L}\left[\pi\left(1-\frac{\delta}{1-\delta}\right)\right]}
  
  whith \eqn{\pi} is the ratio of \code{mutations} by \code{mfn}, \eqn{\delta} the death probability and \eqn{\mathcal{L}(x)} is the Laplace transform of the log-normal distribution with mean \code{mfn} and coefficient of variation \code{cvfn}.
  
}


\value{

    \code{dflan} gives the density, \code{pflan} gives the distribution
     function, \code{qflan} gives the quantile function, and \code{rflan}
     generates a random sample.

     \code{rflan} returns a list with two arguments, each with length \code{n}: a vector of integers 
     \code{$mc} (mutant counts), and a vector of numeric \code{$fn} (final numbers of cells). 
}


\seealso{
  \code{link{draw.clone}}
}

\examples{
#----------------------- distributions ----------------------------------

# Luria-Delbrück model, mean number of mutations 1, fitness parameter 1
dflan(0:10)
pflan(0:10)
qflan(c(0.95,0.99))

# Luria-Delbrück model, mean number of mutations 2, fitness parameter 0.5
qflan(c(0.95,0.99),mutations=2,fitness=0.5)
\donttest{
qflan(c(0.05,0.01),mutations=2,fitness=0.5,lower.tail=FALSE)

# Haldane model, mean number of mutations 2, fitness parameter 0.5
qflan(c(0.95,0.99),mutations=2,fitness=0.5,model="H")
}
#---------------------- random samples ----------------------------------
# lognormal lifetime distribution, lognormal final numbers
X <- rflan(100,cvfn=0.3)
X$mc              # mutant counts
X$fn              # final numbers

# mean number of mutations 2, fitness parameter 0.5 (realistic model, but slow)
\donttest{
rflan(1000,mutations=2,fitness=0.5)$mc
}
# exponential lifetimes (Luria-Delbrück model, fast)
rflan(1000,mutations=2,fitness=0.5,dist=list(name="exp",rate=1))$mc

# constant lifetimes (Haldane model, fast)
rflan(1000,mutations=2,fitness=0.5,dist=list(name="dirac",location=1))$mc

# specifying mutation probability and mean final number
\donttest{
rflan(1000,mutprob=2e-9,mfn=2e9,fitness=0.5,dist=list(name="dirac",location=1))$mc


# positive cell death probability
rflan(1000,mutprob=2e-9,mfn=2e9,death=0.1,fitness=0.5,dist=list(name="dirac",location=1))$mc
}
}

