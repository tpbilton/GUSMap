##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap)
# Copyright 2017-2020 Timothy P. Bilton
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
#' Mapping functions
#' 
#' Transform recobination fractions into map distances using a specified mapping function.
#' 
#' Currrently, the mapping implemented are: 
#' \describe{
#' \item{Morgan mapping function (\code{fun="morgan"})}{\eqn{d = r}}
#' \item{Haldane mapping function (\code{fun="haldane"})}{\eqn{d = -(1/2)ln(1 - 2r)}}
#' \item{Kosambi mapping function (\code{fun="kosambi"})}{\eqn{d = (1/4)ln((1+2r)/(1-2r))}}
#' }
#' where \eqn{r} denotes the recombination fraction value and \eqn{d} denotes the mapping distance value.
#' 
#' @param rf Vector of recombination fraction values
#' @param fun Character vector giving the name of the mapping function (see details).
#' @param centiM Logical value. If \code{TRUE}, the mapping distances returned are in centiMorgans. Otherwise,
#' the mapping distances returned are in Morgans.
#' 
#' @keywords internal
#' @author Timothy P. Bilton
#' @seealso \code{\link{$orderLG}}, \code{\link{$plotLM}}

mfun <- function(rf, fun="morgan", centiM=FALSE){
  if(fun == "morgan")
    dis = rf
  else if (fun == "haldane")
    dis = -0.5 * log(1 - 2 * rf)
  else if (fun == "kosambi")
    dis = 0.25 * log((1 + 2 * rf)/(1 - 2 * rf))
  else
    stop("Map function not defined")
  if(centiM)
    return(100 * dis)
  else
    return(dis)
}
  
  