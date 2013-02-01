##some utility functions

#######################################################################################
##characteristic function of a normally distributed random variable
#######################################################################################
cnorm <- function(t, mean, sd) exp(complex(imaginary = mean*t, real = - 0.5*(sd*t)^2))

#######################################################################################
##characteristic function of a t distributed random variable
#######################################################################################
##from QRMlib or prob package:
ct <- function(t, df){
  aux <- ifelse(t == 0, 0, log(besselK(abs(t)*sqrt(df), df/2, expon.scaled=FALSE) * ((abs(t)*sqrt(df))^(df/2))*(2^(1-df/2))) -lgamma(df/2))
  return(exp(aux))
}

#######################################################################################
##select and evaluate kernel
#######################################################################################
Kernel <- function(t, h, kernel=c("fan", "wand", "sinc"))
  {
    kernel <- match.arg(kernel)

    k <- switch(kernel,
                fan = ifelse(abs(h*t) > 1, 0, (1-(h*t)^2)^3),
                wand = ifelse(abs(h*t) < 1/2 | abs(h*t) > 1, 0, (2*(1-abs(h*t))^3)) +
                ifelse(abs(h*t) > 1/2, 0, (1-6*(h*t)^2+6*abs(h*t)^3)),
                sinc = ifelse(abs(h*t) > 1, 0, (1-abs(h*t))))
    k
  }

