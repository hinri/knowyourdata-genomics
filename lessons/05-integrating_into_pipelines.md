---
layout: lesson
root: .
title: Lesson Title
minutes: 5
---


Getting binary file formats integrated into your pilelines
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

## Where can I get a bamfile of a reasonable size to practice on?
[GSM847334]http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM847334

##Reading from binary files
Lets create a perl script that prints the first 41 lines of a file and we name it maParser.pl 

    while (<STDIN>){
        last if ($count >40);
        $count++;
        print $_;
    }

Then we can read the first 41 lines of a bam file by using our perl script, a tool named [samtools](http://www.htslib.org) and a shell pipe as explained in Software Carpentry tutorial - [The Unix shell](http://software-carpentry.org/v4/shell/index.html)
    
    samtools view -h in.bam |perl myParser.pl

Another option is to intergrate the pipe construction in the perl code like:

    open (Bam, "samtools view -h $ARGV[0] |");
    my $count=0;
    while (<Bam>){
        last if ($count >40);
        $count++;
        print $_;
    }

After modifying the perl script it can be executed like:

    perl myParser.pl in.bam
    
Resulting in the same 21 lines as in the previous example.

##Challenge
    
If an experiment involves multiple bam files we extend the pipeline like this:

    samtools cat in1.bam in2.bam|samtools view - |my Parser.pl

Or by integrating the pipe in the perl script:

    open Bam,"samtools cat  $ARGV[0] |samtools view -h - |";
    my $count=0;
    while (<Bam>){
        last if ($count >40);
        $count++;
        print $_;
    }

Which can be executed as:

    perl myParser.pl "in1.bam in2.bam"
    
Note that we use quotes here to tell the script that it should take the list of bam files as a single argument.


#Writing binary files
Now we want to store the 41 lines of output, that are currently being printed to the screen, in a file using the binary format.

Using the shell pipes the command would look like this:
    samtools view in.bam |perl myParser.pl|samtools -Sb -o out.bam -

Whereas in the perl script we add an output filehandle which decribes the output pipe:

    open (inBam, "samtools view -h $ARGV[0] |");
    open (outBam, "|samtools view -Sb -o out.bam - ");

    my $count=0;
    while (<inBam>){
        last if ($count >40);
	    $count++; 
        print outBam;
    }


#The optional use of compression
Considder the following command in which paired reads RNA-seq data in read1.fastq.gz and read2.fastq.gz is aligned to the reference using a popular multi pupose aligner [gsnap](http://research-pub.gene.com/gmap/)  with a minimal set of options.

    gsnap read1.fastq.gz read2.fastq.gz -d reference -A sam -N --gunzip
    
By adding the --gunzip flag we make gsnap compression aware so there is no need to unzip the input data which saves a lot of (temporary) diskspace.
Gsnap produces gsnap native output by default, however the -A sam tells gsnap to output sam formatted output (which is still not the binary compressed format that we prefer)

Gsnap is (currently) not able to write bam directly so we pipe the out put to samtools view -Sb which makes bam out of sam.
    
    gsnap read1.fastq.gz read2.fastq.gz -d reference -A sam -N --gunzip | samtools view -Sb -o mapped.bam -
    
In the downstream analyses it might be useful to have coordinate sorted bam files which allows for merging and indexing. In order to sort the bamfile we tell samtools sort to read from standard input (-)

    gsnap read1.fastq.gz read2.fastq.gz -d reference -A sam -N --gunzip | samtools view -Sb | samtools sort - mapped.sorted

The odd thing is now that samtools view compresses the output by default and samtools sort is decompressing its input in order to perform sorting. By explicitly telling samtools view to pass the data uncompressed we drastically reduce the number of cpu cycles required to complete the analyses.

    gsnap read1.fastq.gz read2.fastq.gz -d reference -A sam -N --gunzip | samtools view -Sbu | samtools sort - mapped.sorted

#Additional creative use of compression flags
In case you are having a bloody fast temp or scratch space on your computer you might tune algorithms that by default write copressed output to temp but allow for uncompressed writing of temporary files.  

#Advantages of using indexed binary files
Alignment files that are stored in bam and that come with an index (bam.bai) allow for random acces to the data. If you have a bam file but no index, the index file can be created by typing:

    samtools index in.bam
    
This will create an in.bam.bai file which is the index of in.bam

Now we have a bam file with index we can now query for the reads that map in a particular region e.g. Chr5:500000-550000.

    samtools view -h mapped.sorted.bam Chr5:500000-550000


