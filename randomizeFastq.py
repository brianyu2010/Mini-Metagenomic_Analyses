import argparse, random

class randomizeFastq:
    """
    Class Description:
    User Notes:
    Revision History:   2014.12.12 Brian Yu Created
                        2015.05.28 Changed file open syntax to with open() as f
    """

    def __init__(self, fastq1, fastq2):
        """
        Reads in both fastq files, keeps index
        :param fastq1:
        :param fastq2:
        :return:
        """
        with open(fastq1, 'r') as fid:
            self.file1 = []
            for l in fid:
                self.file1.append(l)
        fastq1_length = len(self.file1) / 4
    

        with open(fastq2, 'r') as fid:
            self.file2 = []
            for l in fid:
                self.file2.append(l)
        fastq2_length = len(self.file2) / 4

        assert(fastq1_length == fastq2_length)
        # print "Done Init"


    def randomizeIndices(self):
        """
        Creates a random set of indices assuming corresponding
        lines in self.file1 and self.file2 form paired end reads
        :return: a new set of indices
        """
        l1 = len(self.file1) / 4
        ind = range(l1)
        # print ind
        random.shuffle(ind)
        ind = [4*x for x in ind]
        return ind


    def generateFastq(self, ind, fastq1, fastq2):
        """
        Outputs two new fastq files with specified file names
        :param ind: randomized indices
        :return:
        """

        with open(fastq1, 'w') as fid1, open(fastq2, 'w') as fid2:
            for i in ind:
                fid1.write('%s%s%s%s' %(self.file1[i],self.file1[i+1],self.file1[i+2],self.file1[i+3]))
                fid2.write('%s%s%s%s' %(self.file2[i],self.file2[i+1],self.file2[i+2],self.file2[i+3]))


# When running the script from command line, the following lines are executed
if __name__ == "__main__":
    usage = "USAGE: python randomizeFastq.py fastq1 fastq2 outputname"

    # Making default argument list structures
    p = argparse.ArgumentParser(usage=usage)
    p.add_argument(dest='fastq1', action='store', type=str)
    p.add_argument(dest='fastq2', action='store', type=str)
    p.add_argument(dest='outputname', action='store', type=str)

    arguments = p.parse_args()

    try:
        f = randomizeFastq(arguments.fastq1, arguments.fastq2)
        ind = f.randomizeIndices()
        f.generateFastq(ind, arguments.outputname+'1.fastq', arguments.outputname+'2.fastq')

    except ValueError, e:
        print "ERROR: ValueError %s" % e
        print usage
    except TypeError, e:
        print "ERROR: TypeError %s" % e
        print usage
    except IOError, e:
        print "ERROR: IOError %s" % e
        print usage

