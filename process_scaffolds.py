import argparse, random
from threshold_scaffolds import threshold_scaffolds

class process_scaffolds(threshold_scaffolds):
    """
    Class Description:  Trim scaffolds, inheriting threshold_scaffolds
    User Notes:
    Revision History:   2017.01.22 Brian Yu Created
                        
    """

    
    def __init__(self, inputscaffolds):
        """
        Passes the argument to threshold_scaffolds, uses init method
        of threshold_scaffolds.
        """
        threshold_scaffolds.__init__(self, inputscaffolds)
        print('Number of scaffolds is: %d' % len(self.label))
    

    def trim_beginning(self, num):
        """
        Given number of basepairs num, trim from beginning of scaffolds.
        If too short, discard scaffold. Output is stored and not printed.
        """
        if num > 0:
            numScaffolds = len(self.label)
            temp_label = []
            temp_sequence = []
            for i in range(numScaffolds):
                if len(self.sequence[i]) > num:
                    temp_label.append(self.label[i])
                    temp_sequence.append(self.sequence[i][num:])
            self.label = temp_label
            self.sequence = temp_sequence
        print('Number of scaffolds after trim_beginning is: %d' % len(self.label))


    def trim_end(self, num):
        """
        Given number of bps num, trim num bps from end of scaffolds.
        If too short, discard scaffold. Output is stored in self.
        """
        if num > 0:
            numScaffolds = len(self.label)
            temp_label = []
            temp_sequence = []
            for i in range(numScaffolds):
                if len(self.sequence[i]) > num:
                    temp_label.append(self.label[i])
                    temp_sequence.append(self.sequence[i][:-num])
            self.label = temp_label
            self.sequence = temp_sequence
        print('Number of scaffolds after trim_end is: %d' % len(self.label))


    def trim_to_length(self, num):
        """
        Trim scaffolds to a certain length. If the scaffold is shorter than
        the specified length, the scaffold is not touched. Output is stored
        in self and is not printed.

        """
        if num > 0:
            numScaffolds = len(self.label)
            temp_label = []
            temp_sequence = []
            for i in range(numScaffolds):
                if len(self.sequence[i]) >= num:
                    temp_label.append(self.label[i])
                    temp_sequence.append(self.sequence[i][:num])
            self.label = temp_label
            self.sequence = temp_sequence
        print('Number of scaffolds trim_to_length is: %d' % len(self.label))


    def subsample_by_threshold(self, num):
        """
        Subsample scaffolds by length threshold
        """
        if num > 0:
            numScaffolds = len(self.label)
            temp_label = []
            temp_sequence = []
            for i in range(numScaffolds):
                if len(self.sequence[i]) >= num:
                    temp_label.append(self.label[i])
                    temp_sequence.append(self.sequence[i])
            self.label = temp_label
            self.sequence = temp_sequence
        print('Number of scaffolds after subsample_by_threshold is: %d' % len(self.label))


    def subsample_randomly(self, num):
        """
        Subsample scaffolds randomly
        """
        if num > 0 and num < len(self.label):
            numScaffolds = len(self.label)
            ind = range(numScaffolds)
            random.shuffle(ind)
            subsample_ind = ind[:num]
            self.label = [self.label[x] for x in subsample_ind]
            self.sequence = [self.sequence[x] for x in subsample_ind]
        print('Number of scaffolds subsample_randomly is: %d' % len(self.label))


    def output_scaffolds(self, outputname):
        """
        Write scaffolds to file
        """
        with open(outputname, 'w') as fid:
            # if you add something to the beginning to the contig id then length is no longer #3
            for i in range(len(self.label)):
                # if int(self.label[i].split('_')[3]) >= threshold:
                fid.write('%s\n%s\n' %(self.label[i], self.sequence[i]))


# When running the script from command line, the following lines are executed
if __name__ == "__main__":
    usage = "USAGE: python process_scaffolds.py [Options] input_scaffolds output_scaffolds"

    # Making default argument list structures
    p = argparse.ArgumentParser(usage=usage)
    p.add_argument('-f', '--trimHead', action='store', dest='leadTrim', type=int, default=0) 
    p.add_argument('-n', '--trimTail', action='store', dest='endTrim', type=int, default=0)
    p.add_argument('-l', '--trimToLength', action='store', dest='lengthTrim', type=int, default=0)
    p.add_argument('-t', '--lengthThresh', action='store', dest='lengthThresh', type=int, default=0)
    p.add_argument('-r', '--subsampleDepth', action='store', dest='randomThresh', type=int, default=0)
    p.add_argument(dest='inputname', action='store', type=str)
    p.add_argument(dest='outputname', action='store', type=str)

    arguments = p.parse_args()

    try:
        f = process_scaffolds(arguments.inputname)
        f.trim_beginning(arguments.leadTrim)
        f.trim_end(arguments.endTrim)
        f.subsample_by_threshold(arguments.lengthThresh)
        f.output_scaffolds(arguments.outputname)

    except ValueError, e:
        print "ERROR: ValueError %s" % e
        print usage
    except TypeError, e:
        print "ERROR: TypeError %s" % e
        print usage
    except IOError, e:
        print "ERROR: IOError %s" % e
        print usage

