##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap)
# Copyright 2017-2018 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################
### Transformation and inverse trnasformation for transforming the r.f's on the range [0,0.5] to [-Inf,Inf]
### Author: Timothy Bilton
### Date: 18/01/2017
### Edited: 23/03/2017

## Functions requred for transforming the recombination fraction parameters on the interval [0,0.5]
## to the interval [-inf,inf]
logit2 <- function(p) log(2*p/(1-2*p))

inv.logit2 <- function(logit_p) {
  p <- 1/(2*(1+1/exp(logit_p)))
  p.na <- is.na(p)
  if(sum(p.na)!=0)
    p[which(p.na)] <- 0.5
  return(p)
}

## Functions requred for transforming the recombination fraction parameters on the interval [0,1]
## to the interval [-inf,inf]
logit <- function(p) qlogis(p)
inv.logit <- function(logit_p) plogis(logit_p)

# Alternatives
#logit3 <- function(p) qlogis(2*p)
#inv.logit3 <- function(logit_p) plogis(logit_p)/2
