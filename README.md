# make-large-DISSECT-xml
This R code assists in the assembly of a large XML comprising hundreds of sequence alignments for DISSECT analysis 

The code is simply a series of loops that generate blocks of text that need to be manually inserted into an XML file
in the appropriate places. A start-template can be produced by inserting ~5 sequences in Beauti and setting parameters as usual. This is a somewhat hackish solution to the problem, but is relatively efficient for inserting loggers, sequence alignments etc.
