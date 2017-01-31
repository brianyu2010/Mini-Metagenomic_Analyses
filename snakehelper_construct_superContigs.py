import os, argparse

class snakehelper_construct_superContigs:
    """
    Class Description: Combines a sequence of fasta files representing
                       contigs. Node names will be changed for these contigs.
                       The prefix names will be passed in along with each file
                       Contig labels might not start with >NODE
    User Notes:        This file has been adapted for Snakemake use.
    Revision History:   2014.12.28 Brian Yu
                        2015.05.25 Brian Yu Adapted for Snakemake to operated 
                                            at biosample level renaming contigs
    """

    def __init__(self):
        """
        Create the necessary arrays
        :return:
        """
        self.fasta_label = []
        self.fasta_seq = []


    def read_in_file(self, filename):
        """
        Reads lines from file into a list.
        Add extra sequence identifier and \n if necessary
        :param filename:
        :return:
        """
        # get biosampleID 
        subsampleID = filename.split('.')[-2]

        with open(filename, 'r') as fid:

            self.fasta_label = []
            self.fasta_seq = []
            counter = 0

            for l in fid:

                # if line begins with '>NODE'
                if '>' in l:

                    counter = counter + 1

                    # Check if there's '\n'
                    if '\n' not in l:
                        l = ">SuperContig_" + str(counter) + "_" + l[1:] + "\n" 
                    else:
                        l = ">SuperContig_" + str(counter) + "_" + l[1:]
                    # Append do fasta_node
                    self.fasta_label.append(l)

                # else check 2 tings
                else:

                    if '\n' not in l:
                        l.append('\n')

                    # if length node > length seq then add to append to seq
                    # Basically this is the first line of sequence
                    if len(self.fasta_label) > len(self.fasta_seq):
                        self.fasta_seq.append(l)
                    # else append to existing seq (this is for multi line fasta.
                    else:
                        self.fasta_seq[-1] = self.fasta_seq[-1][0:-1] + l

        assert(len(self.fasta_seq) == len(self.fasta_label))


    def append_to_file(self, filename):
        """
        Append file to existing fasta file. If doesn't exist then create it.
        :param filename:
        :return:
        """
        fid = open(filename, 'a')
        for i in range(len(self.fasta_label)):
            fid.write('%s%s' %(self.fasta_label[i], self.fasta_seq[i]))
        fid.close()


    def write_name_list(self, namefile):
        """
        write name list to file
        :param filename:
        :return:
        """
        self.fasta_label.sort()
        fid = open(namefile, 'a')
        for i in range(len(self.fasta_label)):
            fid.write('%s' %self.fasta_label[i][1:])
        fid.close()


    def output_to_command(self):
        """
        Output joined fasta file to command line
        :return:
        """
        for i in range(len(self.fasta_label)):
            print '%s%s' % (self.fasta_label[i], self.fasta_seq[i])




# When running the script from command line, the following lines are executed
if __name__ == "__main__":
    usage = "USAGE: python snakehelper_construct_superContigs.py input files -o outputfile -l namelist"
    # Call with this line:
    # There should only be one input file

    # Making default argument list structures
    p = argparse.ArgumentParser(usage=usage)
    p.add_argument(dest='inputfile', action='store', type=str)
    p.add_argument('-o',dest='outputfile', action='store', type=str, default='0')
    p.add_argument('-l',dest='namelist', action='store', type=str, default='0')

    arguments = p.parse_args()

    try:
        if os.path.exists(arguments.outputfile):
            os.remove(arguments.outputfile)
        f = snakehelper_construct_superContigs()
        f.read_in_file(arguments.inputfile)
        if arguments.outputfile == '0':
            f.output_to_command()
        else:
            f.append_to_file(arguments.outputfile)
            f.write_name_list(arguments.namelist)

    except ValueError, e:
        print "ERROR: ValueError %s" % e
        print usage
    except TypeError, e:
        print "ERROR: TypeError %s" % e
        print usage
    except IOError, e:
        print "ERROR: IOError %s" % e
        print usage

