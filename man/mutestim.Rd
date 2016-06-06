\name{mutestim}
\alias{mutestim}
\encoding{UTF-8}

\title{Fluctuation Analysis parametric estimation}

\description{
  Estimates mean number of mutations, mutation probability, and fitness parameter, with different methods, under different models. Returns the estimated means and standard deviations for each parameter.
}

\usage{
  mutestim(mc,fn=NULL,mfn=NULL,cvfn=NULL,              # user's data
                  fitness=NULL,death=0,                # user's parameters
                  method=c("ML","GF","P0"),winsor=512, # estimation method
                  model=c("LD","H"))                   # clone growth model
}

\arguments{
  \item{mc}{ a (non-empty) numeric vector of mutants counts. }
  \item{fn}{ an optional (non-empty) numeric vector with same length as \code{mc} of final numbers of cells.}
  \item{mfn}{ mean final number of cells. Ignored if \code{fn} is non-missing. }
  \item{cvfn}{ coefficient of variation of final number of cells. Ignored if \code{fn} is non-missing. }
  \item{fitness}{ fitness parameter: ratio of growth rates of normal and mutant cells. Default 1. If \code{fitness} is \code{NULL}, then the fitness will be estimated. Otherwise, the given value will 
be used to estimate the mean mutation number \code{mutations}. }
  \item{death}{ death probability. Must be smaller than 0.5. } 
  \item{method}{ estimation method as a character string: one of \code{ML} (default), \code{P0}, or \code{GF}. See details. }
  \item{winsor}{ winsorization parameter: positive integer. Only used when \code{method} is \code{ML} or when \code{method} is \code{P0} and \code{fitness} is NULL. See details. }
  \item{model}{ statistical lifetime model as a character string: one of \code{LD} (default) for Luria-Delbrück model (exponential lifetimes), or \code{H} for Haldane model with (constant lifetimes). }
}

\details{
  Method \code{ML} is the classic maximum likelihood estimation method. The maximum is computed with a BFGS (bounded) algorithm. 
  
  Method \code{P0} uses the number of null values in the sample, therefore it can be applied only if there is at least one zero in \code{mc}. The estimate of \code{fitness} is computed by maximum likelihood.

  Method \code{GF} uses the empirical generating function of \code{mc}. Since it is a very fast method, \code{"GF"} is used to initialize the values of the estimates for methods \code{"ML"} and \code{"P0"} (if the fitness is estimated).
  
  If \code{fn} is non-empty, then \code{with.prob} is TRUE. If \code{method} is \code{P0} or \code{GF}, then \code{mfn} and \code{cvfn} are computed from \code{fn}, and the estimate of \code{mutprob} is computed from the estimate of \code{mutations}. If \code{method} is \code{ML}, the estimate of \code{mutprob} is first computed and the estimate of \code{mutations} is deduced.
  
  The winsorization parameter \code{winsor} is used as a threshold for values in \code{mc} when maximum likelihood estimates are computed.
}


\value{ A list containing the following components:
\item{mutations}{ mean number of mutations}
\item{sd.mutations}{ estimated standard deviation on mean number of mutations}
\item{mutprob}{ mutation probability (if \code{with.prob=TRUE})}
\item{sd.mutprob}{ estimated standard deviation on mutation probability}
\item{fitness}{ estimated fitness (if argument \code{fitness} is NULL)}
\item{sd.fitness}{ estimated standard deviation on fitness}
}

\references{
  B. Ycart and N. Veziris: Unbiased estimates of mutation rates under fluctuating final counts. 
PLoS one 9(7) e101434 (2014) 
  
  B. Ycart: Fluctuation analysis with cell deaths.
 J. Applied Probab. Statist, 9(1):12-28 (2014)
 
  B. Ycart: Fluctuation analysis: can estimates be trusted? 
One PLoS one 8(12) e80958 (2013)
 
  A. Hamon and B. Ycart: Statistics for the Luria-Delbrück distribution. 
Elect. J. Statist., 6:1251-1272 (2012)
}

\seealso{
  \code{\link{rflan},\link{flan.test}}.
}

\examples{

# realistic random sample of size 100: mutation probability 1e-9,
# mean final number 1e9, coefficient of variation on final numbers 0.3,
# fitness 0.9, lognormal lifetimes, 5% mutant deaths
x <- rflan(100,mutprob=1e-9,mfn=1e9,cvfn=0.3,fitness=0.9,death=0.05)

# maximum likelihood estimates with mean final number
meanfn <- mutestim(x$mc,mfn=1e9)

# maximum likelihood estimates with final numbers
withfn <- mutestim(x$mc,x$fn)

# change model
Hmodel <- mutestim(x$mc,x$fn,model="H")

# faster methods
GFmethod <- mutestim(x$mc,x$fn,method="GF")
P0method <- mutestim(x$mc,x$fn,method="P0")

# take deaths into account
withdeaths <- mutestim(x$mc,x$fn,death=0.05,method="GF")

# compare results
rbind(meanfn,withfn,Hmodel,GFmethod,P0method,withdeaths)

# use known value of fitness
mutestim(x$mc,x$fn,fitness=0.9)

# extreme example
x <- rflan(10000,mutations=50,fitness=0.5,dist=list(name="exp",rate=1))$mc
summary(x)
mutestim(x,method="GF")
mutestim(x)
mutestim(x,winsor=2000)

\dontrun{
mutestim(x,method="P0")}
}

