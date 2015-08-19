###

Supplemental Electronic Data files for Nature Manuscript 2015-04-04599C:

A Fully Resolved, Comprehensive Phylogeny of Birds (Aves) using Targeted Next Generation DNA Sequencing

Richard O. Prum 1,2,*, Jacob S. Berv 3,*, Alex Dornburg 1,2, Daniel J. Field 2,4, Jeffrey P. Townsend 1,5, Emily Moriarty Lemmon 6, and Alan R. Lemmon 7


1 Department of Ecology & Evolutionary Biology, Yale University, New Haven CT USA
2 Peabody Museum of Natural History, Yale University, New Haven CT USA
3 Department of Ecology & Evolutionary Biology, Cornell University, Ithaca NY USA
4 Department of Geology & Geophysics, Yale University, New Haven CT USA
5 Department of Biostatistics, and Program in Computational Biology and Bioinformatics, Yale University, New Haven, CT USA
6 Department of Biological Science, Florida State University, Tallahassee, FL USA
7 Department of Scientific Computing, Florida State University, Tallahassee, FL USA

 * These authors contributed equally to this work. 

Electronic Data Files prepared by:

Jacob Berv
jsb439@cornell.edu

and 

Alex Dornburg
alex.dornburg@naturalsciences.org

and 

Alan R. Lemmon
alemmon@fsu.edu

###

This directory contains assembled sequence data, newick formatted tree files, code and scripts for generating and analyzing phylogenetic informativeness, probe design, and notes on data assembly. Code used for generating figures is available on request.


~/Trees/best_scheme.txt is the output from PartitionFinder

~/Trees/Concatenated/ contains the final output from all analyses of concatenated sequence data, i.e. RAxML, ExaBayes

~/Trees/Concatenated/BEAST contains the time calibrated output from BEAST (“Avian-TimeTree-Vegavis.tre” includes Vegavis as a fossil calibration, “Avian-TimeTree.tre does not”). See supplemental text for more information

~/Trees/RAxML-GeneTrees/ contains individual gene trees from the 259 locus dataset as estimated and output from RAxML

~/Trees/Species_Trees/ contains the final output from ASTRAL, NJST, and STAR

~/Assembly/Data/ contains the individual locus alignments for the untrimmed and trimmed data in fasta format, as well as concatenated alignments in nexus and phylip formats for the trimmed data

~/Phylogenetic_Informativeness contains all code and files necessary to perform phylogenetic informativeness analyses, see annotations in PI_scripts for instructions

~/Assembly/Code contains the code used for the assembly and alignment generation

~/Assembly/References/ contains the reference files used in the assembly

~/ProbeDesign/ contains files specifying the alignments used in the probe design, the genomic coordinates of the probe design regions, and the probes sequences themselves

