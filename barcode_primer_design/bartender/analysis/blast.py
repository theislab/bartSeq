__author__ = 'Steffen'
import subprocess


def run_blat(query, reads, output):
    """

    :param query: path to query FASTA file
    :param reads:  path to reads FASTA file
    """
    blat_config = ConfigurationHandler.read_blat_config()
    subprocess.call(blat_config.path + "-t=dna q=dna -out=blast8 -minScore=" + blat_config.minScore + " -minIdentity=" +
                    blat_config.minIdentity + " -noHead " + query + " " + reads + " " + output, shell=True)
#
# cmd = paste("/home/icb/steffen.sass/software/blatSrc/bin/amd64/blat -t=dna q=dna -out=blast8 -minScore=0 -minIdentity=95 -noHead ",paste(sequence.library.path,"inserts.fa",sep="/")," ",paste(sequence.library.path,"reads.fa",sep="/")," ", paste(sequence.library.path,"/",name,"_res_inserts.txt",sep=""), sep="")
# blast.out=try(system(cmd, intern = TRUE))