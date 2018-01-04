## Various thermodynamic functions

## ---
## This function returns the saturated vapour pressure in Pascals
## for a temperature in T in K.

## Justin R. Peter CAWCR, BoM, 3 August 2010.

## Uses formula of Emanual (1994), 
## Atmospheric convection, p. 116, eqn(4.4.13).
## T is in K. In the range 0 C to 40 C this is accurate to 0.006%, In the range
## -30 C to 0 C it is accurate to within 0.3% and within 0.7% as low as -40 C.

es <- function(T){
          (exp(53.67957-(6743.769/T)-4.8451*log(T)))*100.
          }

## ---

## This function returns the saturated mixing ratio for a given pressure and
## saturated vapour pressure (both in pascals), in units of kg/kg.


qsat <- function(press,T){
            epsilon<-0.622
            (epsilon*es(T))/(press-es(T))
        }

## ---

## This function calculates the water vapour mixing ratio given
## the pressure and dewpoint temperature. This is a straightforward
## use of qsat.r (above) as the water vapour mixing is the saturation
## mixing ratio at the dewpoint temperature.

r_mix_td <- function(press,Td){
                qsat(press,Td)
            }


#This function calculates the water vapour pressure given the mixing ratio
#as calculated above
e_vap_press <- function(press,r){
                   epsilon = 0.622
                   press*r/(epsilon+r)
               }

cpd <- function(T){
           1005+(T-250)^2./3364
}

#Calculate the latent heat of evaporation
levap <- function(T){
             conv  = 4.18684
             gamma = 0.167+3.67e-4*T
             result<-(597.3*(273.15/T)^gamma)*conv*1000.
             return(result)
}

#Calculate the potential temperature
theta <- function(T,press){
#Convert pressure to hPa
             p  <- press/100.
# Gas constant of dry air, J/kg/K.
             Rd <- 287.04
             chi <- Rd/cpd(T)
             pot_temp<-T*(1000/p)^chi
             return(pot_temp)
}

#Calculate the equivalent potential temperature
theta_e_approx <- function(T,press,r){
                      pot_temp<-theta(T,press)*(1.+levap(T)*r/(cpd(T)*T))
                      return(pot_temp)
}

#Calculate the saturated equivalent potential temperature
theta_es_approx <- function(T,press){
                      pot_temp<-theta(T,press)*exp(levap(T)*qsat(press,T)/(cpd(T)*T))
                      return(pot_temp)
}
                      
#http://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
