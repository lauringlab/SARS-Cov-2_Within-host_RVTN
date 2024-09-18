import os
import subprocess
import argparse
import shutil
import glob



s=("data/aligned_output/consensus/")
f=("data/within_host_consensus/reference/")
print (f)
key= ("Within_host_ref.csv")
print = (key)
#make the output dir if it doesn't exist
if not os.path.exists(f):
    os.makedirs(f)

# add slash to final location if it was forgotten
if f[-1] != '/':
    f=f+'/'
#print f

if s[-1] != '/':
    s=s+'/'
#print(s)

    
# input argument is the sample sheet
junk_names = [] # The names given by the sequencing core sampleid_index
new_names = [] # the names we want, which are in the same row as the data taht is used to make the bad names. 

# so junk_names[i] should be the junk name for new_names[i] as they both come from the ith column.
# add the bad names from the sample sheet to a list. This uses columns of the sheet to construct the names according to the capricous nature of the sequencing core.Although they may have finally settled down a bit.
names =  open(key,"r")
next(names) # skip the header of the csv
for line in names:
    #print (line)
    line = line.strip()
    line = line.split(',')
    junk_names.append(line[0])    
    new=line[4]   
    new_names.append(new)
    #print(junk_names)
names.close()

outfile = open("data/within_host_consensus/Reference_processing.log",'w')

# Now to search through the fastq files or copy them to the new makes
for filename in glob.glob(s + "*.fa"):
    name=filename.split("/")
    bad_path = name[3] # the bad name is the first bit
    bad_file = bad_path.split(".")
    bad_name = bad_file [0]
    #print (bad_name)
    if bad_name in junk_names: # the new name and the bad name have the same index in their respective lists
        name_index = junk_names.index(bad_name)
        better_name= new_names[name_index]
        # Write file to new name
        outfile.write(bad_name + "\t COPIED to \t" + better_name + "\n") # log the move in a renaming_log.txt
        cp_cmd ="cp " + s + bad_name + ".fa "+ f+  better_name + ".fa" 
        subprocess.call(cp_cmd, shell=True)
		
outfile.close()

