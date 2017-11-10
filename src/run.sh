#!/bin/bash
rm fqz* REP* OWL* MINIMIZED* OUTPUT* DECOMPRESSED* FINAL* -fr
# GET OWL =====================================================================
rm -fr owl/
git clone https://github.com/pratas/owl.git
cd owl/src/
make
cp OWL ../../
cd ../../
# GET FQZ_COMP ================================================================
rm -fr fqzcomp/ fqzcomp-4.6/ fqzcomp-4.6.tar.gz
SFURL="https://downloads.sourceforge.net/project/";
wget $SFURL/fqzcomp/fqzcomp-4.6.tar.gz
tar -xzf fqzcomp-4.6.tar.gz
mv fqzcomp-4.6/ fqzcomp/
cd fqzcomp/
make
cp fqz_comp ../
cd ..
# DOWNLOAD ====================================================================
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00097/sequence_read/SRR741385.filt.fastq.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00097/sequence_read/SRR741385_1.filt.fastq.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00097/sequence_read/SRR741385_2.filt.fastq.gz

