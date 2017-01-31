######################################
# Functions used by rules in aseembly
######################################

def split_paired_fastq_into_subfiles(input, output, num_files):
  """
  input: a list containing 2 files representing the two input paired fastq
  output: a list of 2 files representing the two output paired fastq
  Could be optimized better for time
  """
  assert(num_files > 0)
  # these files are fastq so it's 4 lines per read
  r1_lines = sum(1 for l in open(input[0], 'r'))
  r2_lines = sum(1 for l in open(input[1], 'r'))
  print('The number of lines in '+input[0]+' is: '+str(r1_lines))
  print('The number of lines in '+input[1]+' is: '+str(r2_lines))
  assert(check_lines(input[0],input[1])),"Input files have different number of lines."
  # the 4 here is so that each file contains complete reads (4 lines per read)
  readsperfile = int(np.ceil(r1_lines/4/float(num_files)))
  print('The number of reads per file is: '+str(readsperfile))
  assert(readsperfile*4*num_files >= r1_lines)
  # split files
  # first open the 2 inputs. Go throught the inputs once. figure out which output to copy to.
  with open(input[0],'r') as r1, open(input[1],'r') as r2:
    linecnt = 0
    # all the output files should have the same number of lines except the last file 
    for filenum in range(num_files-1):
      print(str(filenum))
      print(output[filenum]+' '+output[filenum+num_files])
      file_line_cnt = 0
      with open(output[filenum],'w') as f1, open(output[filenum+num_files],'w') as f2:
        while file_line_cnt < readsperfile*4:
          t = f1.write(r1.readline())
          t = f2.write(r2.readline())
          linecnt += 1
          file_line_cnt += 1
    assert(file_empty([output[filenum]])),"Output file read 1 is empty."
    assert(file_empty([output[filenum+num_files]])),"Output file read 2 is empty."
    assert(check_fastq_ids(output[filenum],output[filenum+num_files])),"Output files have different read id orders."  
    assert(linecnt < r1_lines),"Line count is larger than number of lines in fastq file."
    # Now put the rest of the lines in the last split file for output (note the num_files-1)
    with open(output[num_files-1],'w') as f1, open(output[num_files-1+num_files],'w') as f2:
      while linecnt < r1_lines:
        t = f1.write(r1.readline())
        t = f2.write(r2.readline())
        linecnt += 1
    assert(file_empty([output[filenum]])),"Output file read 1 is empty."
    assert(file_empty([output[filenum+num_files]])),"Output file read 2 is empty."
    assert(check_fastq_ids(output[filenum],output[filenum+num_files])),"Output files have different read id orders."


def extract_subsample_fastq_segment(input, output, num_files, segment):
  # This function extracts the nth segment from 
  # two fastq input files and puts them in 2 corresponding
  # fastq output files. Checks these output files.
  assert(num_files > 0),"Number of files to split into must be positive."
  assert(type(segment) is int),"Segment has to be an integer."
  assert(segment <= num_files),"Segment must be smaller or equal to total file numbers."
  # these files are fastq so it's 4 lines per read
  r1_lines = sum(1 for l in open(input[0], 'r'))
  r2_lines = sum(1 for l in open(input[1], 'r'))
  print('The number of lines in '+input[0]+' is: '+str(r1_lines))
  print('The number of lines in '+input[1]+' is: '+str(r2_lines))
  assert(check_lines(input[0],input[1])),"Input files have different number of lines."
  # the 4 here is so that each file contains complete reads (4 lines per read)
  readsperfile = int(np.ceil(r1_lines/4/float(num_files)))
  print('The number of reads per file is: '+str(readsperfile))
  assert(readsperfile*4*num_files >= r1_lines)
  # split files
  print('Segment is: '+str(segment))
  # Split R1 files
  with open(output[0],'w') as f1, open(input[0],'r') as r1:
    linecnt = 0
    file_line_cnt = 0
    for line in r1:
      # print(line)
      if linecnt >= segment*readsperfile*4 and linecnt < (segment+1)*readsperfile*4:
        t = f1.write(line)
        linecnt += 1
        file_line_cnt += 1
      else: 
        linecnt += 1
  print('The number of lines in '+output[0]+' is: '+str(file_line_cnt))
  # Split R2 files
  with open(output[1],'w') as f2, open(input[1],'r') as r2:
    linecnt = 0
    file_line_cnt = 0
    for line in r2:
      # print(line)
      if linecnt >= segment*readsperfile*4 and linecnt < (segment+1)*readsperfile*4:
        t = f2.write(line)
        linecnt += 1
        file_line_cnt += 1
      else: 
        linecnt += 1
  print('The number of lines in '+output[1]+' is: '+str(file_line_cnt))
  assert(check_lines(output[0],output[1])),"Output files have different number of lines."
  assert(file_empty(output)),"One of the output files is empty."
  assert(check_fastq_ids(output[0],output[1])),"Output files have different read id orders."
  


def check_fastq_ids(file1, file2):
  # This function checks sequence id of 2
  # fastq files and sees if they are the same
  with open(file1,'r') as f1:
    id1 = []
    linecnt = 0
    for l in f1:
      if linecnt % 4 == 0:
        id1.append(l.split()[0])
        linecnt += 1
      else:
        linecnt += 1
  with open(file2,'r') as f2:
    id2 = []
    linecnt = 0
    for l in f2:
      if linecnt % 4 == 0:
        id2.append(l.split()[0])
        linecnt += 1
      else:
        linecnt += 1
  if len(id1) == len(id2):
    for i in range(len(id1)):
      if id1[i] != id2[i]:
        print('Different fastq ids '+id1[i]+' and '+id2[i])
        return 0
  else:
    print('Different fastq file lengths.')
    return 0
  return 1


def file_empty(files):
  # check if files are empty
  # imput is a list of files
  for f in files:
    if sum(1 for l in open(f,'r')) == 0:
      return 0
  return 1


def check_lines(file1, file2):
  """
  check that the number of lines in 2 files are the same
  """
  l1 = sum(1 for l in open(file1,'r'))
  l2 = sum(1 for l in open(file2,'r'))
  if l1 == l2:
    return 1
  else:
    return 0


def merge_fastq_files(input, output):
  """
  both input and output are list of strings representing file names
  input is of the form ['r1_0.fastq','r1_1.fastq','r2_0.fastq','r2_1.fastq'...]
  """
  numfiles = len(input)/2
  print('Number of files for read1 and read2 detected is: '+str(numfiles))
  with open(output[0],'w') as f1, open(output[1],'w') as f2:
    for i in range(int(numfiles)):
      with open(input[i],'r') as in1:
        linecnt1 = 0
        for l in in1:
          f1.write(l)
          linecnt1 += 1
      print('File '+input[i]+' has '+str(linecnt1)+' lines.')
      with open(input[i+int(numfiles)],'r') as in2:
        linecnt2 = 0
        for l in in2:
          f2.write(l)
          linecnt2 += 1
      print('File '+input[i+int(numfiles)]+' has '+str(linecnt2)+' lines.')
      assert(check_lines(input[i],input[i+int(numfiles)])),"Input files have different number of lines."
  assert(check_lines(output[0],output[1])),"Output files have different number of lines."
  assert(file_empty(output)),"One of the output files is empty."
  assert(check_fastq_ids(output[0],output[1])),"Output files have different read id orders."


def get_split_filenames(wildcards):
  # This functions generates a list filesnames for 
  # the subsample split step. Different number of 
  # files for different subsamples
  numfiles = int(depth_table.ix[wildcards.subsample,'subsample_clust_file_numbers'])
  R1_names = []
  R2_names = []
  for i in range(numfiles):
    R1_names.append(wildcards.subsample+'/clusterR1_'+str(i)+'.fastq')
    R2_names.append(wildcards.subsample+'/clusterR2_'+str(i)+'.fastq')
  return R1_names+R2_names


def unzip_files(input):
  # This function takes in a list of input *.gz filenames
  # and unzips them. The input files must be a list, even
  # if it's just one file. The output files will not have
  # the .gz ending and will replace the .gz version in place
  for file in input:
    if os.path.isfile(file):
      subprocess.call(['gzip','-d',file])


def zip_back_files(input):
  # This function takes a list of modified input filenames
  # and zips files that are unzipped before back into 
  # the zipped versions. It can be used at the end of a rule
  # This function also checks if a file exists
  for file in input:
    if os.path.isfile(file):
      subprocess.call(['gzip',file])

def modify_zip_file_names(input):
  # This function removes the .gz ending from all input file names
  modified_filename = []
  for file in input:
    modified_filename.append(file[0:-3]) # gets rid of .gz ending
  return modified_filename
