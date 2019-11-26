##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap)
# Copyright 2017-2019 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
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
#
#' Genotyping Uncertainty with Sequencing data for linkage Mapping (GUSMap)
#' 
#' Package for constructing genetic maps using outcrossed full-sib family
#' populations that have been sequenced using next generation multiplexing
#' sequencing methods \insertCite{@see @bilton2018genetics1}{GUSMap}.
#'
#' \tabular{ll}{ Package: \tab GUSMap\cr Type: \tab Package\cr Version: \tab
#' 2.0.0\cr Date: \tab 2019-11-27\cr License: \tab GPL 3\cr }
#'
#' @name GUSMap-package
#' @aliases GUSMap
#' @docType package
#' @author Timothy P. Bilton, Maintainer: Timothy P. Bilton
#' <tbilton@@maths.otago.ac.nz>
#' @references
#' \insertAllCited{}
#' @keywords package
#' @importFrom Rdpack reprompt
#' @importFrom R6 R6Class
#' @importFrom GUSbase logit
#' @importFrom GUSbase inv.logit
#' @importFrom GUSbase logit2
#' @importFrom GUSbase inv.logit2
#' @importFrom Rcpp sourceCpp
#' @import plotly foreach
NULL



