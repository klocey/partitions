rpartitions--Code for integer partitioning
=========================================

rpartitions is an R package is focused on the integer partitioning problem of randomly partitioning
some total q into n parts. These functions solve the main computational challenge of Locey and White (2013)
and extend the 'feasible set' approach based on integer partitioning to several ecological patterns and 
other distributions. 


Algorithms were derived by Ken Locey, and coded into R by Dan McGlinn. Ken Locey and Dan McGlinn both work to improve the readability
and organization of the code.

Installation
============
rpartitions, the R package, can be installed either from CRAN or directly from git hub. 

To install the most recent stable version from CRAN use the following R command:

    install.packages('rpartitions')
    
To install the most updated version directly from github use the following R commands:

    library(devtools)
    install_github('partitions', 'klocey', subdir='rpartitions')

Altneratively the pacakge can be downloaded from github and locally installed using:

    library(devtools)
    install()

GNU GENERAL PUBLIC LICENSE
==========================
Version 2, June 1991

Copyright (C) Kenneth J. Locey and Daniel J. McGlinn

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

