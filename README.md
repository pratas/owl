<p align="center"><img src="imgs/owl.png"
alt="OWL" width="264" border="0" /></p>
<br>

<p align="justify">
Under development...
</b>

## 1. INSTALLATION ##

Downloading and installing OWL:
<pre>
git clone https://github.com/pratas/owl.git
cd owl/src/
cmake .
make
</pre>

Cmake is needed for the installation (http://www.cmake.org/). You can download it directly from http://www.cmake.org/cmake/resources/software.html or use an appropriate packet manager, such as:
<pre>
sudo apt-get install cmake
</pre>

## 2. USAGE ##

To see the possible options of OWL type
<pre>
./OWL
</pre>
or
<pre>
./OWL -h
</pre>
These will print the following options:
<pre>
<p>
Usage: OWL [OPTIONS]... [FILE] [FILE]                                    
A tool to sort FASTQ reads using cluster mapping.                        
                                                                         
Non-mandatory arguments:                                                 
                                                                         
  -h                         give this help,                             
  -V                         display version number,                     
  -v                         verbose mode (more information),            
  -N                         does NOT sort reads,                        
  -W                         writes the full header,                     
  -D                         does NOT delete the temporary file,         
  -k <k-mer>                 k-mer size [1;20],                          
  -m <minimum>               minimum block size,                         
                                                                         
Mandatory arguments:                                                     
                                                                         
  <FILE>                     reference file,                             
                                                                         
  <  <FILE>                  stdin input FASTQ file,                     
  >  <FILE>                  stdout output sorted FASTQ file,            
                                                                         
Example:                                                                 
                                                                         
  ./OWL -v -k 16 -m 40 reference.fa < ex1.fq > ex1-sort.fq               
                                                                         
Report bugs to <pratas@ua.pt>.                            
</pre>
All the parameters can be better explained trough the following table:

| Parameters          | Meaning                                                     |
|---------------------|:------------------------------------------------------------|
| -h                  | It will print the parameters menu (help menu)                                        |
| -V                  | It will print the OWL version number, license type and authors information.    |
| -v                  | It will print progress information.    |
| -N                  | It will NOT sort the reads (for analysis purposes). |
| -W                  | It will write the full header in the output FASTQ file. Usually a very part of the header is not needed.    |
| -D                  | It will not delete the temporary file for ordering the reads (for analysis purposes).    |
| -k &#60;k-mer&#62;   | word size of the slidding window. From 1 to 20. Usually, larger values need more memory.    |
| -m &#60;minimum&#62;      | minimum size of proximity. Used in the elastic clustering.              |
| [FILE]           | Reference filename (DNA sequence or FASTA file). |

## 3. CITATION ##

Under development!

## 4. ISSUES ##

For any issue let us know at [issues link](https://github.com/pratas/owl/issues).

## 5. LICENSE ##

GPL v3.

For more information:
<pre>http://www.gnu.org/licenses/gpl-3.0.html</pre>

