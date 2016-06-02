# diagflu
Automating Prokka annotations, BLAST runs and generation of BLAST reports

## Dependencies

The script will work for Linux. Windows usage has not been implemented yet. The installation of a Linux virtual machine (preferrably Ubuntu) on your Windows computer is strongly recommended for the meanwhile.

In order to run this script, you would need to have the following installed on your system:

1. __Python 2.7__

2. __NCBI BLAST+__
  
  This is the standalone version of BLAST that you can run locally on your machine. `diagflu.py` calls the `blastn` algorithm.

  Instructions for Linux:

  `sudo apt-get install ncbi-blast+`
  
3. __Prokka: microbial annotation software__
  
  Please follow the installation instructions in the official [Prokka repository](https://github.com/tseemann/prokka)
  
## Set up 

`diagflu.py` will do a BLAST search of your annotated sequences vs. a local database. You need to build said database. The instructions to do this are:

1. Download the sequences from NCBI or [FluDB](fludb.org) that you want to compare to your contigs. In my case I downloaded all the Influenza Virus A sequences that existed in NCBI on May 26, 2016 (taxid:11320). I have a [Perl script](https://gist.github.com/ropolomx/1155bf740716d488f83b6f905fc2327d) that you can use to do that from the command line.
2. Once you have the database you want, you need to type in the following commands:

  `makeblastdb -in filename.fasta -dbtype nucl -title filename -out filename`
3. You are good to go.

## Usage

`python diagflu.py --blast-db /path/to/blastdb -r name_of_report_file.txt -p name_of_prokka_folder contigs.fasta`

## Sample data

The folder contains the file `La_Habana_test.fasta` which contains the 8 segments of the A/swine/La/Habana/130/2010/H1N1 strain.
