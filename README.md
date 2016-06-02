# FluMatch
Automating Prokka annotations, BLAST runs and generation of BLAST reports

## Dependencies

The script will work and has been tested on Linux only. Windows implementation is being investigated, but it is not available yet. The installation of a Linux virtual machine (preferrably Ubuntu) on your Windows computer is strongly recommended for the meanwhile.

In order to run this script, you would need to have the following installed on your system:

1. __Python 2.7__

2. __NCBI BLAST+__
  
  This is the standalone version of BLAST that you can run locally on your machine. `flumatch.py` calls the `blastn` algorithm.

  Instructions for Linux:

  `sudo apt-get install ncbi-blast+`
  
3. __Prokka: microbial annotation software__
  
  Please follow the installation instructions in the official [Prokka repository](https://github.com/tseemann/prokka)

## Process

The script does the following:

1. Reads in your FASTA format file with contig sequences
2. 
  
## BLAST database set up 

`diagflu.py` will do a BLAST search of your annotated sequences vs. a local database. You need to build said database. The instructions to do this are:

1. Download the sequences from NCBI or [FluDB](fludb.org) that you want to compare to your contigs. In my case I downloaded all the Influenza Virus A sequences that existed in NCBI on May 26, 2016 (taxid:11320). I have a [Perl script](https://gist.github.com/ropolomx/1155bf740716d488f83b6f905fc2327d) that you can use to do that from the command line.
2. Once you have the database you want, you need to type in the following commands:

  `makeblastdb -in filename.fasta -dbtype nucl -title filename -out filename`

   For example, if you have a FASTA file named `avian.fasta` or `avian.fna` which contains the sequences that you want to use as your BLAST database, you could do the following:

  `makeblastdb -in avian.fasta -dbtype nucl -title avian -out avian`

  After you run that command you will see that three new files will be created in the directory with the extensions `*.nhr`, `*.nin`, and `*.nsq`. These are the files that BLAST will use to search your sequences against.

## Usage

`python flumatch.py --blast-db /path/to/blastdb -r name_of_report_file.txt -p name_of_prokka_folder contigs.fasta`

### Arguments and options:

The most important thing to keep in mind is that the name of the FASTA file that contains the contigs you want to analyze should be written at the very end of the command.

`--blast-db`: __This argument is required__ The local BLAST database that you want to search your annotated sequences against. See the __Set up__ section above.

`-p` or `--prokka-dir`: The name of the directory where the Prokka output will be stored. __Default =__ the program will create a sub-directory with the name of the file you are 

`-t` or `--top-hits` : The number of top BLAST hits for each annotated CDS that you want to see in the final report. __Default = 10__

`-r` or `--report-out`: The name of the report text file. This is a tab-separated file that you can open in R or Excel. __Default = `TopBLASThits.txt`__

## Sample data

I have included a folder with sample data to try the script. The folder contains the following files and directories 

* `La_Habana_test.fasta` which contains the 8 segments of the publicly available strain A/swine/La/Habana/130/2010/H1N1 (Genbank Accession numbers HE584753.1 to HE584760.1).

* `TopBLASThits.txt`: an example report of the BLAST report of the annotated strain. The version of `blastn` used for this was `2.2.28+`. The current version of online BLAST is `2.3.1+`.

* A directory called __La\_Habana\_test__ that has the output of the __Prokka annotation process__
