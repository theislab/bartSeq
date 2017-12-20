from __future__ import print_function

class PrimerPrinter:
    @staticmethod
    def print_set(gene_set, blast = True):
        for gene in gene_set:
            print(">"+gene.name+":\n")
            for i, amplicon in enumerate(gene.amplicons):
                print("sequence "+str(i + 1)+":\n\n")
                for primerpair in amplicon.primer_set.set:
                    print(primerpair.fwd_seq() + " | " + primerpair.rev_seq() + " ")
                    if blast:
                        print("(" + str(primerpair.fwd.blast_hits) + " | " + str(primerpair.rev.blast_hits) + ")\n")
                    else:
                        print("\n")
                print("\n")
            print("\n\n")