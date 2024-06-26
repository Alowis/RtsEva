% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{tsEasyParseNamedArgs}
\alias{tsEasyParseNamedArgs}
\title{Parse named arguments and assign values to a predefined argument structure.}
\usage{
tsEasyParseNamedArgs(args, argStruct)
}
\arguments{
\item{args}{A list of named arguments.}

\item{argStruct}{A list representing the argument structure with named elements.}
}
\value{
A modified argument structure with values assigned from the list of named arguments.
}
\description{
This function takes a list of named arguments and assigns their values to a predefined argument structure.
The argument structure is a list with named elements representing the available arguments.
If an argument is present in the list of named arguments, its value is assigned to the corresponding element in the argument structure.
If an argument is not present, its value in the argument structure remains unchanged.
}
\examples{
args <- list(arg1 = 10, arg2 = "tanargue")
argStruct <- list(arg1 = 0, arg2 = "", arg3 = TRUE)
modifiedArgStruct <- tsEasyParseNamedArgs(args, argStruct)
modifiedArgStruct

}
