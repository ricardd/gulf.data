#' Length Weight Regression
#'
#' Performs length-weight regressions.
#'
#' A linear regression is performed on log-transformed variables. Weights of
#' one gram or less are removed from the analysis, due to the fact that small
#' fish with low weights are often rounded up to one gram, which can skew the
#' analysis.
#'
#' @param x A data frame with fields \sQuote{length} and \sQuote{weight}.
#' @param species A numeric vector of species codes of the data to be loaded if
#' \code{x} is not specified.
#' @param sex A numeric vector of specifying the sex of the data to be loaded
#' if \code{x} is not specified. Sex code \code{9} is generally interpreted as
#' being a pooled sex category, i.e. all sexes combined.
#' @param year A numeric vector of specifying the year of the data to be loaded
#' if \code{x} is not specified. Specifying this argument signifies that
#' length-weight coefficients are to be calculated from survey data rather than
#' looked up in a table.
#' @param survey A character string specifying the name of the survey data to
#' be loaded. This argument is passed onto the \code{\link[gulf]{read.gulf}}
#' function.
#' @param by A character vector specifying the names of grouping variables
#' (e.g.  \sQuote{sex}, \sQuote{year}, \sQuote{species}, ...)  by which to
#' perform the analyses.
#' @param log10 A logical values specifying whether the \sQuote{a}
#' length-weight coefficient should log-10 base rather then the natural
#' logarithm.
#' @param units Character value specifying the units of the weights. It may be
#' either in grams (\code{units = "g"}) or kilograms (\code{units = "kg"}).
#' @return A matrix of regression coefficients is produced. These correspond to
#' the scaling (\sQuote{a}) and the exponential (\sQuote{b}) parameters of the
#' allometric model.
#'
#' If \code{by} is specified, then its variables values are included in the
#' output.
#' @seealso \code{\link[gulf]{read.gulf}}
#' @examples
#'
#'    # Plaice coefficients:
#'    length.weight(species = 40, sex = 0:2, by = "sex") # From reference table.
#'    length.weight(species = 40, by = "sex", year = 2015) # From RV 2015, separated by sex.
#'    length.weight(species = 40, sex = 1:2, year = 2015) # From RV 2015,
#'      # grouped males and females.
#'    length.weight(species = 40, sex = 1:2, year = 2015, by = "sex") # From RV 2015,
#'      # separated males and females.
#'    length.weight(species = 40, year = 2015) # Pooled coefficients.
#'    length.weight(species = 40, year = 2015, sex = 9) # Same answer as previous example.
#'
#'    # Cod length-weight coefficients from 2013 RV data:
#'    length.weight(species = 10, year = 2013)
#'
#'    # Length-weight coefficients for multiple species:
#'    length.weight(species = c(10, 12, 40, 41, 42, 43), sex = 2)
#'
#'    # Length-weight coefficients for multiple species and years:
#'    length.weight(species = c(10, 12, 40, 41, 42, 43), year = 2010:2013,
#'      by = c("species", "year", "sex"))
#'
#'    # Length-weight coefficients for American plaice for years (2010-2013) combined:
#'    length.weight(species = 40, year = 2010:2013)
#'    length.weight(species = 40, year = 2010:2013, by = "sex")
#'
#'    # Load data separately and calculate coefficients:
#'    x <- read.gulf.bio(year = 2013, species = 40)
#'    length.weight(x, by = "sex")
#'
#' @export length.weight
length.weight <- function(x, species, sex, year, survey = "rv", by, log10 = TRUE, units = "g"){
   # LENGTH.WEIGHT - Extract length-weight regression coefficients.

   # Parse 'x' argument:
   if (!missing(x)){
      if (!is.data.frame(x)) stop("'x' must be a data frame.")
      names(x) <- tolower(names(x))
      if (!all(c("length", "weight") %in% names(x))) stop("'x' must contain 'length' and 'weight' fields.")
      if (missing(species) & ("species" %in% names(x))) species <- sort(unique(x$species))
   }

   # Parse 'species' argument:
   if (missing(species)) stop("'species' must be specified.")

   # Parse 'units' argument:
   units <- match.arg(tolower(units), c("grams", "kg", "kilograms"))
   if (units == "kg") units <- "kilograms"

   # Load data:
   if (missing(x) & !missing(species) & !missing(year))
   {
      x <- read.gulf.bio(species = species, year = year, survey = survey, password="R_GulfPKG_19")
#      if(year == 2018)
 #        x <- special.condition.2018(x)
   }
   if (missing(x)){
      # Look-up coefficients from table:
      data(length.weight.coefficients, envir = environment())
      x <- length.weight.coefficients
      if (missing(sex)) sex <- 9 else sex[sex == 4] <- 9
      t <- data.frame(species = species, sex = sex)
      index <- (x$species %in% t$species) & (x$sex %in% t$sex)
      if (length(which(index)) == 0){
         cat("No match found for length-weight coefficients found in reference table (no pooled estimate?).\n")
         return(NULL)
      }else{
         x <- x[index, ]
         x <- x[c("species", "sex", "a", "b", "number")]
         names(x) <- c("species", "sex", "a", "b", "n")
         if (!log10) x$a <- log(10^x$a)
         if ((units == "kilograms") & log10) x$a <- x$a - 3
         if ((units == "kilograms") & !log10) x$a <- x$a - log(1000)
         return(x)
      }
   }else{
      # Isolate subset of 'x':
      if (!missing(sex) & ("sex" %in% names(x))) if (!any(sex == 9)) x <- x[x$sex %in% sex, ]
      if (!missing(species) & ("species" %in% names(x))) x <- x[x$species %in% species, ]

#      if(year == 2018)
 #        x <- special.condition.2018(x)
      # Remove invalid data:
      index <- is.na(x$length) | is.na(x$weight) | (x$length <= 1) | (x$weight <= 1)
      x <- x[!index, ]

      # Convert to appropriate units:
      if (units == "kilograms") x$weight <- x$weight / 1000

      # Partition data set:
      if (!missing(by)){
         flag <- TRUE
         flag <- flag & ("sex" %in% by) & ("sex" %in% names(x))
         if (flag) flag <- flag & all(x$sex != 9 | is.na(x$sex ))
         if (flag & !missing(sex)) flag <- flag & any(sex == 9)
         if (flag){
            temp <- x
            temp$sex <- 9
            if (!missing(sex)) x <- x[x$sex %in% sex, ]
            x <- rbind(x, temp)
         }
         data <- by(x, x[by], function(x) x)
      }else{
         data <- list(x)
      }

      # Perform linear regression:
      v <- matrix(NA, nrow = length(data), ncol = 4)
      for (i in 1:length(data)){
         if (!is.null(dim(data[[i]])[1])){
            m <- stats::lm(log(weight) ~ log(length), data = data[[i]])
            v[i, 1] <- stats::coef(m)[1]
            v[i, 2] <- stats::coef(m)[2]
            v[i, 3] <- summary(m)$sigma
            v[i, 4] <- nrow(data[[i]])
         }
      }
      dimnames(v) <- list(1:dim(v)[1], c("a", "b", "sigma", "n"))

      # Convert to data frame:
      v <- as.data.frame(v)

      # Convert to base 10:
      if (log10) v$a <- base::log10(exp(v$a))

      # Add category labels:
      if (!missing(by)) v <- cbind(expand.grid(dimnames(data)), v)

      return(v)
   }
}
special.condition.2018 <- function(x){
   # Special condition put here to eliminate cods smaller than 15cm.
   # In the past these specimens were coded as 251 but this year they were not.
   index = x$year == 2018 & x$species == 10 & x$length < 15
   if(length(which(index)) > 0)
      x = x[!index,]
   return(x)
}
