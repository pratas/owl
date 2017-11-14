<p align="center"><img src="imgs/logo.png"
alt="OWL" width="2666666rder="0" /></p>
<br>

<p align="center"><b>OWL: a tool to order FASTQ reads using elastic cluster mapping.</b></p></br>
<p align="justify">
OWL is a new tool to order FASTQ reads, neglecting the original order. It maps the reads according to a reference sequence using k-mer positional hashing and, then, it orders the reads using elastic clustering. Its usage is only needed during compression, enabling a very fast and low memory decompression. The tool can be used to substantially improve the compression of the FASTQ files (See the following Figure for a pipeline using the general purpose GZIP compressor). The time complexity of the tool is approximately linear.<br></p>
<p align="center"><img src="imgs/description.png"
alt="OWL" width="670" border="0" /></p>
<br>
A human reference genome can be downloaded using the script GetHuman.sh contained in the scripts folder.


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
An alternative to cmake, but limited to Linux, can be set using the following instructions:
<pre>
cp Makefile.linux Makefile
make
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
A tool to order FASTQ reads using elastic cluster mapping.                        
                                                                         
Non-mandatory arguments:                                                 
                                                                         
  -h                         give this help,                             
  -V                         display version number,                     
  -v                         verbose mode (more information),            
  -N                         does NOT sort reads,                        
  -W                         writes the full header,                     
  -D                         does NOT delete the temporary file,         
  -k &#60k-mer&#62                 k-mer size [1;20],                          
  -m &#60minimum&#62               minimum block size.
                                                                         
Mandatory arguments:                                                     
                                                                         
  &#60FILE&#62                     reference file,                             
                                                                         
  &#60  &#60FILE&#62                  stdin input FASTQ file,                     
  &#62  &#60FILE&#62                  stdout output sorted FASTQ file.
                                                                         
Example:                                                                 
                                                                         
  ./OWL -v -k 16 -m 40 reference.fa &#60 ex1.fq &#62 ex1-sort.fq               
                                                                         
Report bugs to &#60{pratas,ap}@ua.pt&#62.                            
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
| -k &#60;k-mer&#62;   | The word size of the slidding window. From 1 to 20. Usually, larger values need more memory.    |
| -m &#60;minimum&#62;      | The minimum size of proximity. Used in the elastic clustering.              |
| [FILE]           | Reference filename (DNA sequence or FASTA file). |
| &#60; [FILE]           | Input FASTQ file with the arbitrary read order (standard input). |
| &#62; [FILE]           | Output FASTQ file with the reads ordered (standard output). |

## 3. CITATION ##

On using this tool/method, please, cite:

D. Pratas, A. J. Pinho (2017). v1.1 pratas/owl: A tool to order FASTQ reads using elastic cluster mapping.<br>
DOI: 10.5281/zenodo.1048947 

## 4. ISSUES ##

For any issue let us know at [issues link](https://github.com/pratas/owl/issues).

## 5. LICENSE ##

GPL v3.

For more information:
<pre>http://www.gnu.org/licenses/gpl-3.0.html</pre>

