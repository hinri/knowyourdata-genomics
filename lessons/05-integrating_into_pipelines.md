---
layout: lesson
root: .
title: Lesson Title
minutes: 5
---


Getting to intergrate binary formats into your pilelines
===================

Learning Objectives:
-------------------
#### What's the goal for this lesson?
* You should get familiar with techniques to use binary data formats in custom pipelines and scripts
* You should be able to read and write binary data formats using your pipelines and scripts
* You should become aware that binary compression is optional
* You should get sense whether compression of binary output should be applied or not
* You should become aware about indexing of binary files


#### At the end of this lesson you should be able to:
* You should be able to make data pipelines that use of binary data format
* You should be able to read and write binary data formats using your scripts.
* You should be able to make the right decission in applying compression
* You should be abe to make effective use of binary file indexing

## Why using binary file formats?
Binary formats of a particular genomics/genetics data type have advantages over text format files:
These advantages include:
* File consisency checks (data type and truncation)
* Storage wise efficient
* Indexed for fast searches


##Reading from binary files
Lets assume we have simple sam file parser written in perl (maParser.pl) 

    while (<STDIN>){
        last if ($count >20);
        $count++;
        print $_;
    }
    
Then we can read bam just by using tool "samtools" and a shell pipe as explained in Software Carpentry tutorial - [The Unix shell](http://software-carpentry.org/v4/shell/index.html)
    
    samtools view -h in.bam |perl myParser.pl

Another option is intergrating the pipe construction in the perl code like:

    open (Bam, "samtools view -h $ARGV[0] |");
    my $count=0;
    while (<Bam>){
        last if ($count >20);
        $count++;
        print $_;
    }

Resulting perl script can be executed as:

    perl myParser.pl in.bam 

If multiple bam files belonging to the same experiment are involved we extend the pipeline like this:

    samtools cat in1.bam in2.bam|samtools view - |my Parser.pl

Or by integrating the pipe in the perl script:

    open Bam,"samtools cat  $ARGV[0] |samtools view -h - |";
    my $count=0;
    while (<Bam>){
        last if ($count >20);
        $count++;
        print $_;
    }

Which can be executed as:

    perl myParser.pl "in1.bam in2.bam"
    
Note that we use quotes here to tell the script that it should take the list of bam files as a single argument.





#Writing binary files
Now we want to store the parsed results, that are currently being printed to the screen, in a file using the binary format.

Using the shell pipes the command would look like this:
    samtools view in.bam |perl myParser.pl|samtools -Sb -o out.bam -

Whereas in the perl script we add an output filehandle which also decribes the output pipe:

    open (inBam, "samtools view -h $ARGV[0] |");
    open (outBam, "|samtools view -Sb -o $ARGV[0].parsed - ");

    my $count=0;
    while (<inBam>){
        last if ($count >20);
	    $count++; 
        print outBam;
    }


#The optional use of compression
Considder the following command in which paired reads RNA-seq data in read1.fastq.gz and read2.fastq.gz is aligned to the reference using gsnap[gsnap](http://research-pub.gene.com/gmap/)  with a minimal set of options.

    gsnap read1.fastq.gz read2.fastq.gz -d refrence -A sam -N --gunzip | samtools view -Sbu - |samtools sort -m 8000000000 - mapped.sorted.bam"
    
We make gsnap compression aware by using the optional --gunzip flag so there is no need to unzip the input data which saves a lot of (temporary) diskspace. Furthermore we tell gsnap to write sam (instead of gsnap native) output. Gsnap is (currently) not able to write bam directly so we pipe the out put to samtools view -Sb which makes bam out of sam. In the downstream analyses it might be useful to have coordinate sorted bam files which allows for merging and indexing. In order to sort the bamfile we tell samtools sort to read from standard input (-) and tel samtools view to pass the data UNCOMPRESSED (-u) to standard output. By explicitly telling samtools view to pass the data ucompressed we drastically reduce the number of required cpu cycles to complete the anlyses. Otherwise samtools view would have compressed the output (computation intensive) and samtools sort would have decompressed it anyway in order to perform sorting.

#Additional creative use of compression flags
In case you are having a bloody fast temp or scratch space on your computer you might tune algorithms that by default write copressed output to temp but allow for uncompressed writing of temporary files.  

#Advantages of using indexed binary files
Alignment files stored in bam and that come with an index (bam.bai) allow for random acces to the data. As an example we show here hw to query for the reads that map in region Chr5:500000-550000.
    samtools view -h mapped.sorted.bam Chr5:500000-550000


#more complex examples:
On the fly trimming of paired end reads prior to alignment. Major adavantage is that trimming and alignment occurs simultaniously and no "trimmed data" is written to disk.

finding pairs in fastq data and put them in the same order.


