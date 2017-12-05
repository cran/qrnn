\name{dummy.code}
\alias{dummy.code}
\title{
Convert a factor to a matrix of dummy codes
}
\description{
Converts a factor (categorical) variable to a matrix of dummy codes
using a 1 of C-1 binary coding scheme.
}
\usage{
dummy.code(x)
}
\arguments{
  \item{x}{
      a factor variable.
  }
}
\value{
      a matrix with the number of rows equal to the number of cases in \code{x}
      and the number of columns equal to one minus the number of factors in
      \code{x}. The last factor serves as the reference group.
}
\examples{
  print(dummy.code(iris$Species))
}
