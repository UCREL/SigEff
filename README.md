# SigEff
Significance and Effect size measures for corpus comparison

Developed by Paul Rayson at Lancaster University.

Part of Wmatrix and previously Xmatrix and Tmatrix.

See http://ucrel.lancs.ac.uk/llwizard.html for a description of the measures.

This repository contains both C code and the spreadsheet versions.

To compile the C version:

gcc -g -o sigeff sigeff.c -lm

To run the C version:

./sigeff -X semtag_summary.txt < spreadsheet > spreadsheet.sgf

Where X indicates the effect size measures to show from the list below:
1 = %DIFF
2 = Bayes
3 = ELL
4 = RRisk
5 = LogRatio
6 = OddRatio
0 = All six

