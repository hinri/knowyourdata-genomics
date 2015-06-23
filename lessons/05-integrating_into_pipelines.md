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

## Why using binary formats?
Binary formats of a particular genomics/genetics data type have advantages over text format files:
These advantages include:
* File consisency checks (data type and truncation)
* Storage wise efficient
* Indexed for fast searches


##

#reading from binary
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





#writing binary
Now we want to store the parsed results, that are currently printed to the screen, in a file using the binary format

Using the shell pipes the command would look like this:
    samtools view in.bam |perl myParser.pl|samtools -Sb -o out.bam -

Whereas in in the perl script we add an output filehandle which als decribeds the output pipe

open (inBam, "samtools view -h $ARGV[0] |");
open (outBam, "|samtools view -Sb -o $ARGV[0].parsed - ");

my $count=0;
while (<inBam>){
    last if ($count >20);
	$count++; 
    print outBam;
}





#writing binary


more complex examples:
On the fly trimming of paired end reads prior to alignment. Major adavantage is that trimming and alignment occurs simultaniously and no "trimmed data" is written to disk.

finding pairs in fastq data and put them in the same order.


Enter the command:

    git clone https://github.com/tracykteal/tutorials/
