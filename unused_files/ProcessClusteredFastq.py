import argparse

class ProcessClusteredFastq:
    """
    Description:    This script can be used to re-organize paired end fastq files after
                    dnaclust so that the pairedness of the reads are maintained. To achieve
                    this, each cluster file is processed to generate a list of read tags.
                    Then, each fastq file is parsed to select out every 4 lines corresponding
                    to those selected tags and stored into new fastq files.
    Method:         Currently, only the first field of the dnaclust output is considered.
    Revision:       2014.12.11 Brian Yu
                    2015.05.21 Edited to allow user to input 2 different output names
                    2015.05.27 Changed all the open(file) to with open(file) as f
    """

    def __init__(self, r1clust, r2clust):
        """
        Read in all the read tags.
        For now I'm only reading in the first column of the tags (sorted)
        :return:
        """
        with open(r1clust, 'r') as fid:
            r1 = []
            for l in fid:
                # for now only take the first element
                r1.append(l.split()[0])

        with open(r2clust, 'r') as fid:
            r2 = []
            for l in fid:
                # for now only take the first element
                r2.append(l.split()[0])
        
        # join the 2 sets - unsorted
        self.selected_tags = set(r1) | set(r2)
        print 'Processing clustered fastq: tag organization completed'


    def rearrange_fastq(self, fastq1_in, fastq2_in, out1_fastq, out2_fastq):
        """
        Takes in 2 cluster files and finds 2 reduced list of fastq files
        :param fastq1: original fastq1 file
        :param fastq2: original fastq2 file
        :param out1_fastq: name of the output R1 file. 
        :param out2_fastq: name of the output R2 file.
        :return:
        """
        # process fastq1,
        self.grep_fastq(fastq1_in, out1_fastq)
        print 'Processing clustered fastq: %s completed' %(out1_fastq)
        # process fastq2
        self.grep_fastq(fastq2_in, out2_fastq)
        print 'Processing clustered fastq: %s completed' %(out2_fastq)


    def grep_fastq(self, fastq_in, fastq_out):
        """
        Takes a input fastq, pull out all the lines related to tags found in self.selected_tags
        Output the file fastq_out
        :param fastq_in:
        :param fastq_out:
        :return:
        """
        # print fastq_in
        # print fastq_out

        with open(fastq_in, 'r') as fid_in, open(fastq_out, 'w') as fid_out:
            l = fid_in.readline()
            linenum = 0
            while l:
                if linenum == 0:
                    # grab the tag (fastq file tags have @ at beginning
                    fastq_tag = l.split()[0].split('@')[1]
                    # if tag is in the selected set
                    if fastq_tag in self.selected_tags:
                        # write 4 lines to new fastq file
                        l2 = fid_in.readline()
                        l3 = fid_in.readline()
                        l4 = fid_in.readline()
                        linenum = (linenum + 3) % 4
                        fid_out.write('%s%s%s%s' %(l, l2, l3, l4))
                l = fid_in.readline()
                linenum = (linenum + 1) % 4


# When running the script from command line, the following lines are executed
if __name__ == "__main__":
    usage = "USAGE: python ProcessClusteredFastq.py clustered1 clustered2 fastq1 fastq2 out1 out2"

    # Making default argument list structures
    p = argparse.ArgumentParser(usage=usage)
    p.add_argument(dest='R1cluster', action='store', type=str)
    p.add_argument(dest='R2cluster', action='store', type=str)
    p.add_argument(dest='fastq1', action='store', type=str)
    p.add_argument(dest='fastq2', action='store', type=str)
    p.add_argument(dest='out1', action='store', type=str)
    p.add_argument(dest='out2', action='store', type=str)

    arguments = p.parse_args()

    try:
        seqtags = ProcessClusteredFastq(arguments.R1cluster, arguments.R2cluster)
        seqtags.rearrange_fastq(arguments.fastq1, arguments.fastq2, arguments.out1, arguments.out2)

    except ValueError, e:
        print "ERROR: ValueError %s" % e
        print usage
    except TypeError, e:
        print "ERROR: TypeError %s" % e
        print usage
    except IOError, e:
        print "ERROR: IOError %s" % e
        print usage


