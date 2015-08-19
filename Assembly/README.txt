Note that the scripts given will soon be published separately as an integrated assembly package. See www.anchoredphylogeny.com for updates.

*.java files need to be complied using e.g. javac Assembler.java

~16 GB of ram is required for assembly

Raw fastq files for each individual must be contained in a folder with the following format: e.g. I1001, with one folder per individual

The merger, which must be run before the assembler, will merge paired, overlapping reads and produce fastq files named I*/I*_M.fastq I*/I*_U1.fastq I*/I*_U2.fastq

The folder structure required is like the following:

Code
I1001
I1002
I1003
<more individual folders>
References

Run the programs from within the java folder

Contact alemmon@fsu.edu for more information