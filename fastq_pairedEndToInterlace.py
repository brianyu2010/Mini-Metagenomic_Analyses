import argparse, random

class pairedFastqToInterlace:
    """
    Class Description: Joins paired end fastq files saved in separate files into one
                       fastq file with interlaced format. 
    User Notes:
    Revision History:   2017.01.16 Brian Yu Created
                                   I do not check if input files are valid
    """

    def __init__(self, fastq1, fastq2, outputfile):
        """
        Reads in both fastq files, keeps index
        :param fastq1:
        :param fastq2:
        :return:
        """
        with open(fastq1, 'r') as fid:
            self.linenum = 0
            for l in fid:
                self.linenum += 1

        with open(fastq1, 'r') as f1, open(fastq2, 'r') as f2, open(outputfile, 'w') as fid:
            for i in range(self.linenum / 4):
                l1 = f1.readline()
                l2 = f1.readline()
                l3 = f1.readline()
                l4 = f1.readline()
                l5 = f2.readline()
                l6 = f2.readline()
                l7 = f2.readline()
                l8 = f2.readline()
                fid.write('%s%s%s%s%s%s%s%s' %(l1,l2,l3,l4,l5,l6,l7,l8))
                if i % 1000000 == 1:
                    print('.')
            print('\n')
            

# When running the script from command line, the following lines are executed
if __name__ == "__main__":
    usage = "USAGE: python fastq_pairedEndToInterlace.py fastq1 fastq2 outputname"

    # Making default argument list structures
    p = argparse.ArgumentParser(usage=usage)
    p.add_argument(dest='fastq1', action='store', type=str)
    p.add_argument(dest='fastq2', action='store', type=str)
    p.add_argument(dest='outputname', action='store', type=str)

    arguments = p.parse_args()

    try:
        pairedFastqToInterlace(arguments.fastq1, arguments.fastq2, arguments.outputname)

    except ValueError, e:
        print "ERROR: ValueError %s" % e
        print usage
    except TypeError, e:
        print "ERROR: TypeError %s" % e
        print usage
    except IOError, e:
        print "ERROR: IOError %s" % e
        print usage

