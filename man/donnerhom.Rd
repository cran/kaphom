\name{donnerhom}
\alias{donnerhom}
\title{
    Donner GOF test for homogeneity of kappa statistics
}
\description{
    Donner GOF test for homogeneity of kappa statistics
}
\usage{
    donnerhom(pp, pm, mm)
}
\arguments{
    \item{pp}{
        a number vector carrying the number of positive maching rates for each study
    }
    \item{pm}{
        a number vector carrying the number of non-maching rates for each study
    }
    \item{mm}{
        a number vector carrying the number of negative maching rates for each study
    }
}
\details{
    This function can be used only for studies with binary output
}
\value{
    This function prints the Donner GOF test statistic, P-value for hypothesis test of the statistic and whether the difference is statistically
    significant or not.
}
\author{
    Muammer ALBAYRAK
}
\examples{

    library(kaphom)

    pp <- c(11,26,22)
    pm <- c(6,5,14)
    mm <- c(22,10,39)

    donnerhom(pp,pm,mm)
}
