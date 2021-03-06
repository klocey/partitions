partitions--Code for integer partitioning
=========================================

This repository contains Python and R packages for the partitioning of integers with focus on generating random partitions of a total q into n parts.

The repository also contains Python scripts for generating figures in the following manuscript:
Locey KJ, McGlinn DJ. (2013) Efficient algorithms for sampling feasible sets of macroecological patterns. PeerJ PrePrints 1:e78v1

Algorithms were derived by Ken Locey and coded into R by Dan McGlinn.

Files & Folders
---------------
**./partitions/pypartitions/**

partitions.py --contains the primary partitioning functions.

test.py --contains functions that test the function of partitions.py

testfiles --folder containing .txt files used in test.py


**./partitions/Locey\_McGlinn\_ms/**

Locey\_McGlinn\_2013.py --code for regenerating figures and conducting analyses in Locey and McGlinn (2013)

time_files.py --code to generate files that write the amount of time the new algorithms take to generate random partitions into the time_files subdirectory.

time_files --folder containing .txt files generated by time_files.py


**./partitions/metrics/**

metrics.py --functions for computing evenness, diversity, inequality of an integer partition


Installation
============
pypartitions, the python package, can be installed by cloning the directory and running the following from the command line:

    python setup.py install

rpartitions, the R package, can be installed either from CRAN or directly from git hub. 

To install the most recent stable version from CRAN use the following R command:

    install.packages('rpartitions')
    
To install the most updated version directly from github use the following R commands:

    library(devtools)
    install_github('partitions', 'klocey', subdir='rpartitions')

Altneratively the pacakge can be downloaded from github and locally installed using:

    library(devtools)
    install('rpartitions')

License

partitions is licensed under the open source MIT License

Copyright (c) 2013 Kenneth J Locey

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

-------------------
Ken Locey's email: ken@weecology.org and locey@biology.usu.edu

*Go to [Ken's website](http://kenlocey.weecology.org)*
