#!/usr/bin/python3
import glob
import os
import re
import shutil
import sys
from pathlib import Path
import pandas as pd

#Importing necessary modules

## Intro to Script ##
os.system('clear')
print('This script is for sequence conservation analysis of user-defined protein sequences within a user-defined taxon.\nThe maximum allowed number of sequences is 1000.\nThe minumum allowed number of sequences is 3.\nInput your search paramaters below.')

## START OF ESEARCH ##
esearch_path1 = os.environ['HOME']
Path(esearch_path1+'/ICA2').mkdir(parents = True, exist_ok=True)
os.chdir(esearch_path1+'/ICA2')
esearch_path2 = os.getcwd()
esearch_path3 = esearch_path2+'/esearch_output'
#Defining the paths to be used in the script
if os.path.exists(esearch_path3):
    shutil.rmtree(esearch_path3)
os.makedirs(esearch_path3)
#Removes the output directory if it exits, then creates the output directory and any parent directories necessary

def esearch_input () :
#Defines the esearch_input function.

    protein_name_in= input('Enter protein name: ')
    taxonID= input('Enter taxon ID: ')
#Takes user input for protein name, then taxon ID

    print('Searching protein database, please wait...')
#Tells user the search is sarting

    protein_name_arg = protein_name_in.replace(" ", "_")
#Replaces the spaces with underscore so I can use it to name the output file

    esearch_cmd1 =  f"esearch -db protein -query '{protein_name_in}[PROTEIN]' 2>~/ICA2/error.txt | efilter -query txid'{taxonID}'[ORGANISM] 2>~/ICA2/error.txt | grep '<Count>' > ~/ICA2/esearch_temp.txt"
#First esearch command used to find the number of sequences. 2>~/ICA2/error.txt creates an error.txt file to pass the errors to if user enters an incorrect sequence

    esearch_cmd2 =  f"esearch -db protein -query '{protein_name_in}[PROTEIN]' 2>/dev/null | efilter -query txid'{taxonID}'[ORGANISM] 2>/dev/null | efetch -format fasta 2>/dev/null > ~/ICA2/esearch_output/{protein_name_arg}_{taxonID}.fasta"
#esearch command that takes the input from the user as the arguements for the search. 2>/dev/null mutes error messages as to not spam the screen if incorrect input has been entered. efilter allows for the search to be narrowed to a single taxon. Used f-string formatting to define the variable.

    os.system(esearch_cmd1)
#Calling the esearch_cmd variable defined above

    if os.stat(esearch_path1+'/ICA2/error.txt').st_size == 0:
        esearch_temp1 = open(esearch_path2+'/esearch_temp.txt')
        esearch_temp2 = esearch_temp1.read()
        esearch_temp3 = re.split(r"<|>", esearch_temp2) #<Count>234</Count>
        esearch_temp4 = int(esearch_temp3[2])
        os.remove(esearch_path2+'/error.txt')
#If the erorr.txt file has content (if the esearch command returned an error) then the sequence count is found and defined as the esearch_temp4 varaible. The intinsure it is an integer before being passed to the next if statment. Regex is used to split the <Count> line in the output file made by the esearch_cmd1 variable at '<' and '>'. The index of the number of sequences is known, which can then be defined to a varible to be used below.

        if esearch_temp4 > 1000:
            print('WARNING: Too many sequences detected, please narrow search paramters or try a different query.')
            os.remove(esearch_path2+'/esearch_temp.txt')
            return esearch_input()
#If the number of sequences is greater than 1000, then it prints there were too many seqeunces and returns out of the function after deleting the temp file.

        elif esearch_temp4 < 3:
            print('WARNING: Too few seuqences detected, please try a different query.')
            os.remove(esearch_path2+'/esearch_temp.txt')
            return esearch_input()
#If there were less than 3 sequence (Minimum for clusalto according to ebi website) then it prints a warning and returns the function after deleting the temp file.

        else:
            os.remove(esearch_path2+'/esearch_temp.txt')
            os.system(esearch_cmd2)
            return protein_name_arg, taxonID
#If there at more than 3 but less than 1000 sequences, then esearch_cmd2 is ran which fetches the file and outputs it in fasta format after deleting the temp file.
    else:
        os.remove(esearch_path2+'/error.txt')
        print('Invalid search input, please try again.')
        return esearch_input()
#Ending the function

protein_name_arg, taxonID= esearch_input()
#Calling the function

output_file_name = protein_name_arg+'_'+taxonID+'.fasta' 
#Assigning the output_file_name variable.

## REMOVE PARTIAL SEQUENCES ##
def remove_partial_question(question3= 'Do you want to remove partial sequences?'):
    reply = str(input(question3+' [y/n]: ')).lower().strip()
    if reply[0] == 'y':
        print('Removing partial sequences, please wait...')
        partial_lines = list()
        with open(esearch_path3+'/'+output_file_name) as partial_file:
            next_line = False
            for partial_line in partial_file.readlines():
                if next_line:
                    next_line = False
                    continue
                if "partial" in partial_line:
                    next_line = True
                    continue
                partial_lines.append(partial_line)
        with open(esearch_path3+'/'+output_file_name, "w") as partial_file:
            partial_file.writelines(partial_lines)
        return True
    if reply[0] == 'n':
        print('Continuing without removing partial seuqences...')
        return True
    else:
        return remove_partial_question("Invalid response, please try again.")
remove_partial_question()
#Asks the user if they would like to remove partial sequences. If they answer yes, then the code above find the lines with "partial" and removes them, along with the sequence by overwiting the original fasta file with a new one containing only non-partial sequences. 

## SKIP REDUNDANT ##
skip_redundant_arg1 = esearch_path3+'/'+output_file_name
skip_redundant_cmd = f"skipredundant -sequences {skip_redundant_arg1} -minthreshold 30.0 -maxthreshold 90.0 -gapopen 10.0 -gapextend 0.5 -redundantoutseq '' -mode 2 -outseq {skip_redundant_arg1}"
#Skip redudant removes redudnant sequences from fasta files based on defined thresholds.I chose all the default skipredundant values. I wrote the output file to the same file path as the fasta input file to replace the original file with the trimmed one to make subsequent analysis easier in the rest of the script. 

def skip_redundant_question(question3= 'Do you want to remove sequences with over 90% and below 30% similarity?'):
    reply = str(input(question3+' [y/n]: ')).lower().strip()
    if reply[0] == 'y':
        print('Removing sequences, this may take a while, please wait...')
        os.system(skip_redundant_cmd)
        return True
    if reply[0] == 'n':
        print('Continuing without removing seuqences...')
        return True
    else:
        return skip_redundant_question("Invalid response, please try again.")
skip_redundant_question()
#Asks the user if they would like to remove sequences with greater than 90% similarity and less than 30% similarity. These are the default values of skipredundant. If they choose yes the sequences are removed and the script continues. If they choose no then the sequences are left as they are and kept for subsequent analysis. 

print("Fasta file has been saved in esearch_output directory as "+output_file_name)
#If the file contains text, it tells the user the file has been saved in the esearch_output directory

## START OF FASTA -> DATAFRAME/CSV ##
fasta_path1 = os.environ['HOME']
os.chdir(fasta_path1+'/ICA2')
fasta_path2 = os.getcwd()
Path(fasta_path2+'/csv_file_dir').mkdir(parents=True, exist_ok=True)
shutil.rmtree(fasta_path2+'/csv_file_dir')
os.mkdir(fasta_path2+'/csv_file_dir')
fasta_path3 = fasta_path2+'/csv_file_dir'
#Since ~ wont work for homespace, I assign the path using os.environ['HOME'] + the name of the directroy I want to target to generalise the script

fasta_file = output_file_name
#Listing the directory gives the name of the file in list format
my_file = open(esearch_path3+'/'+fasta_file) 
#opening the fasta file for the script to work on

id = []
name = []
organism = []
seq = []
#Defining the lists I will be using in the for loop below

for eachline in my_file:
    if eachline.startswith('>'):
#This checks if the line starts which ">" indicating it is the title of the fasta sequence
        line = eachline.split()
        id.append(line[0])
#Splits the line into a list, then since the accession number is the first item in the list, I add it directly to the ID list defined earlier
        edit = re.search('(\[.+\])', eachline)
        line1 = eachline.replace(edit.group(1), '')
        line2 = line1.split()
        line3 = line2[1:]
        line4 = ' '.join(line3) 
        name.append(line4)
#Since the protein name can vary in terms of number of words, and so can the organism name, this is requried to generliase the script. It finds the organism name with regex. It then removes the organism name from the line, so that the only remaining are the ID and the protein name. Since protein name will always be 2nd place in the list, I can append from 2nd place until the end of the list to the name list defined earlier without adding the organsim name too. the join line simply make it look nicer in the title of the columns in the dataframe and csv made below
        organism.append(edit.group(1)) 
#Regex found the organism name above, this appends the name to the predefined list
    else:
        seq.append(eachline)
#If the line does not start with ">" then it is a sequence, so it simply appends the sequence to the sequence list
	
s_name = pd.Series(name)
s_id = pd.Series(id)
s_organism = pd.Series(organism)
s_seq = pd.Series(seq)
#Converts the lists to series

df = pd.DataFrame({'Accessoin Number' : s_id, 'Protein Name' : s_name, 'Organism Name' : s_organism, 'Sequence' : s_seq})
#Creates a dataframe with the series

Organism_count = df['Organism Name'].value_counts()
print(Organism_count)
number_of_organisms = len(df['Organism Name'].value_counts())
#Prints the organism name column to the screen, due panda being nice it ranks them in order of frequency automatically. Then prints the number of organisms in the file.

def user_input(question= 'There are '+str(number_of_organisms)+' organisms represented in the FASTA file, do you want to continue to clustalo alignment?'):
    reply = str(input(question+' [y/n]: ')).lower().strip()
    if reply[0] == 'y':
        csv_arg = fasta_path3+'/'+fasta_file+'.csv'
        df.to_csv(csv_arg, sep='\t')
#Outputs dataframe to csv file.

        return True
    if reply[0] == 'n':
        return sys.exit()
    else:
        return user_input("Invalid response, please try again.")
user_input()
#Asks user if they would like to continue to clustalo alignment based on the number of organsims represented in the fasta file.

## START OF CLUSTALO ##
print('Starting clustalo alignment, please wait...')
#Tells the user clustalo alignment is starting

clustalo_path1 = os.environ['HOME']
os.chdir(clustalo_path1+'/ICA2')
clustalo_path2 = os.getcwd()
clustalo_path3 = clustalo_path2+'/clusalto_output'
if os.path.exists(clustalo_path3) :
    shutil.rmtree(clustalo_path3)
Path(clustalo_path3).mkdir(parents=True, exist_ok=True)
fasta_file_name = fasta_file
fasta_file_name_only = fasta_file_name.replace('.fasta', '')
#Defining variables that will be used for clustalo command.

clustalo_arg1 = esearch_path3+'/'+fasta_file_name
clustalo_arg2 = clustalo_path3+'/'+fasta_file_name_only+'.msf'
#If I put these paths into the below clustalo_input variable without defining them as variables first it would not work. Defined these variables as a fix.

clustalo_input = f"clustalo -i {clustalo_arg1} -o {clustalo_arg2} --outfmt=msf --wrap=80 --force --threads=32"
os.system(clustalo_input)
#Defines the clustalo input as a variable and runs it with os.system. --outfmt=msf outputs in msf format for plotcon to work. --threads=32 is to speed it up by multithreading. --force overwrites the file if it already exists. --wrap=80 tells it to allow 80 residues before a line wrap in the output file.

print('clustalo alignment complete, saving msf file to clustalo_output directory...')
#Prints that clustalo alignment is complete and has been saved to the output directory.

## INFO ALIGN ##
if os.path.exists(clustalo_path2+'/info_align_results'):
    shutil.rmtree(clustalo_path2+'/info_align_results')
Path(clustalo_path2+'/info_align_results').mkdir(parents=True, exist_ok=True)
infoalign_out = clustalo_path2+'/info_align_results/'+fasta_file_name_only
infoalign_input1 = f"infoalign -sequence {clustalo_arg2} -outfile {infoalign_out}.infoalign -only -change -idcount -simcount -diffcount -heading"
infoalign_input2 = f"infoalign -sequence {clustalo_arg2} -outfile {infoalign_out}.infoalign -only -name -change -idcount -simcount -diffcount -heading"
#Defining the variables to be used for infoalign emboss programme. It shows basic information about sequence alignment in a text file. I thought this might be useful the user to get extra information about the alignment of the protein sequences. The first one is for creating the dataframe, the second command is used for saving the infoalign output file for the user to view if it is needed.

os.system(infoalign_input1)
#Running the first command variable 

df1 = pd.read_csv(infoalign_out+'.infoalign', sep='\t')
Ident = str(df1['# Ident'].agg(['max','min','mean','std']))
Similar = str(df1['Similar'].agg(['max','min','mean','std']))
Differ = str(df1['Differ'].agg(['max','min','mean','std']))
Change = str(df1['% Change'].agg(['max','min','mean','std']))
#Assigning the aggregates of the max, min, mean and standard deviation to their respective variables to be printed to the user below

print('Aggregate of number of identical residues =\n'+Ident+'\nAggregate of number of similar residues =\n'+Similar+'\nAggregate of number of different residues =\n '+Differ+'\nAggregate of percentage changed positions =\n'+Change) 
os.system(infoalign_input2)
#Running the 2nd command variable
print('An infoalign file has been created displaying information about the sequence alignment in the info_align_ esults directory.\nNavigate to this directory to view results.')
#Prints where to find the infoalign output file to the user.
#This is my "wildcard" option. It converts the info align output file into a dataframe, then uses that dataframe to get an aggregate of stats for the sequences. It then prints that info to the user, which gives extra, easily accessible info about the sequence alignments to the user.

## SHOW ALIGN ##
if os.path.exists(clustalo_path2+'show_align_results'):
    shutil.rmtree(clustalo_path2+'show_align_results')
Path(clustalo_path2+'/show_align_results').mkdir(parents=True, exist_ok=True)
show_align_out = clustalo_path2+'/show_align_results/'+fasta_file_name_only
show_align_input = f"showalign -sequence {clustalo_arg2} -outfile {show_align_out}.showalign"
#Defining the variables to be used for showalign. This is an emboss tool that visually represents the alignmnet of the protein sequences.

os.system(show_align_input)
print('A show align file has been created displaying the alignment of the protein sequences visually in the show_align_results directory.\nNavigate to this directory to view results')
#Running the command variable and telling the user how to find the output file for showalign.

def clustalo_input_question(question3= 'Do you want to continue to plotcon?'):
    reply = str(input(question3+' [y/n]: ')).lower().strip()
    if reply[0] == 'y':
        return True
    if reply[0] == 'n':
        print('Exiting script...')
        return sys.exit()
    else:
        return clustalo_input_question("Invalid response, please try again.")
clustalo_input_question()
#Asks user if they would like to continue to plotcon analysis, if not then the script exists.

## START OF PLOTCON ##
print('Starting plotcon analysis, please wait...')
#Tells user that plotcon is starting.

plotcon_path1 = os.environ['HOME'] 
os.chdir(plotcon_path1+'/ICA2')
plotcon_path2 = plotcon_path1+'/ICA2/plotcon_output'
if os.path.exists(plotcon_path2):
    shutil.rmtree(plotcon_path2)
Path(plotcon_path2).mkdir(parents=True, exist_ok=True)
msf_file_name = fasta_file_name_only+'.msf'
msf_file_name_only = msf_file_name.replace('.fasta.msf', '')
#Defining the variables I will use in the plotcon command.

plotcon_arg1 = clustalo_arg2
plotcon_arg2 = plotcon_path2+'/'+msf_file_name_only
#Similarly to the clustalo command, plotcon would not work without me defining these paths as variables before being used in the plotcon_input variables. 

plotcon_input1 = f"plotcon -sequences {plotcon_arg1} -graph png -winsize 4 -goutfile {plotcon_arg2} >/dev/null"
#First plotcon input that is run. It saves the plot as a png in the plotcon_output directory.

plotcon_input2 = f"plotcon -sequences {plotcon_arg1} -graph x11 -winsize 4 -goutfile {plotcon_arg2} >/dev/null"
#Second plotcon input that is run. Instead of saving the plot, it runs it as a pop up on the screen if the display supports it.

os.system(plotcon_input1) 
print('PNG of plot saved to plotcon_output directory')
#Calling the first input command for plotcon and telling user where the file was saved.

def plotcon_input_question(question3= 'Do you want to view the plot now?'):
    reply = str(input(question3+' [y/n]: ')).lower().strip()
    if reply[0] == 'y':
        os.system(plotcon_input2)
        return True
    if reply[0] == 'n':
        print('To view plot, download the png file from the plotcon_output directory')
        return False
    else:
        return plotcon_input_question("Invalid response, please try again.")
plotcon_input_question()
#Asks the user if they would like the view the plotcon plot now. If the answer is "y", then the 2nd plotcon input is run and a pop up of the plot is displayed if supported.

def prosite_continue_question(question3= 'Do you want to scan protein sequences against PROSITE motif database?'):
    reply = str(input(question3+' [y/n]: ')).lower().strip()
    if reply[0] == 'y':
        return True
    if reply[0] == 'n':
        print('Exiting script...')
        return sys.exit()
    else:
        return prosite_continue_question("Invalid response, please try again.")
prosite_continue_question()
#Asks the user if they want to continue to PROSITE motif database scan.

## START OF PROSITE DATABASE SCAN ##
os.chdir(esearch_path1+'/ICA2')
prosite_path2 = os.getcwd()
prosite_path3 = prosite_path2+'/split_fasta_dir'
prosite_path4 = esearch_path3+'/'+output_file_name
if os.path.exists(prosite_path3):
    shutil.rmtree(prosite_path3)
Path(prosite_path3).mkdir(parents=True, exist_ok=True)
#Defining the variables used for seqretsplit.

seqretsplit_cmd = f"seqretsplit -sequence {prosite_path4} -outseq *.fasta -osdirectory2 {prosite_path3}"
os.system(seqretsplit_cmd)
#Using seqretsplit to split the sinlge fasta file into individual fasta files with 1 sequence each as patmatmotif only accepts 1 sequence at a time.

if os.path.exists(prosite_path2+'/PROSITE_dir'):
        shutil.rmtree(prosite_path2+'/PROSITE_dir')
Path(prosite_path2+'/PROSITE_dir').mkdir(parents=True, exist_ok=True)
prosite_path5 = prosite_path2+'/PROSITE_dir'
pathlist = Path(prosite_path3).glob('*.fasta')
#Assigning variables to be used below. Used Path module to iterate throughthe the split fasta files.

print('Scanning protein sequences against PROSITE database of motifs...')
#Print statment to tell the user PROSITE database scan is begining.

for files in pathlist:
    patmatmotif_cmd = f"patmatmotifs -sequence {files} -outfile {files}.dbmotif"
    os.system(patmatmotif_cmd)
for f in glob.iglob(prosite_path3+'/*.dbmotif'):
    shutil.move(f, prosite_path5)
#Used a loop to iterate through the files in the split_fasta_dir directory and pass them through patmatmotif to scan for the PROSITE database motifs. Used glob module to select the files ending in .dbmotif, which are the output of patmatmotif, to move them from the split_fasta_dir to the PROSITE_out directory.I had issues with patmatmotif output so this was the best solution I found that didn't involve using os.system() with a bash command.

prosite_pathlist = Path(prosite_path5).glob('*')
move_path = prosite_path2+'/PROSITE_hits'
if os.path.exists(move_path):
        shutil.rmtree(move_path)
Path(move_path).mkdir(parents=True, exist_ok=True)
#Defining varibles that will be used below to move all files with hits to a seperate directory.
for i in prosite_pathlist:
    open_file = open(i)
    read_file = open_file.read()
    list_file = read_file.split('\n')
#Loops through files in the output directory from the patmatmotif search, opens the file and splits the file based on new lines.

    if any("HitCount" in s for s in list_file):
#Using the any function to identiy substring in the list of lines in the file to find the line that contains "HitCount".

        matching = [s for s in list_file if "HitCount" in s]
#Adds any index in the list that contains "HitCount" to the list "matching".

        hitcount_join = ' '.join(matching)
        hitcount_list = hitcount_join.split(':')
        hitcount_int = int(hitcount_list[1])
#join the matching list, then split again using  ':' as the split character as I know the number of counts will come after the ':',and only 1 ':' is in the line, then the count number will be the 2nd position in the list, aka the first index. Making sure the count number is an integer as well.

        if hitcount_int > 0:
            shutil.move(str(i), move_path)
#If the count number is more than 0, then the file gets moved to the PROSITE_hits directory for the user to view. str() had to be used to convert the path back to a string as shutil will only recognise it as a string.

print('Sequences with PROSITE motif database hits have been saved to PROSITE_hits directory.\nThe accession number of the sequence is the file name of the .dbmotif file.')
#Telling the user where to find the sequences with the PROSITE_hits and how the files are named.

## PEP STATS ##
if os.path.exists(prosite_path2+'/pep_stats_results'):
    shutil.rmtree(prosite_path2+'/pep_stats_results')
Path(prosite_path2+'/pep_stats_results').mkdir(parents=True, exist_ok=True)
pepstats_input = f"pepstats -sequence {clustalo_arg2} -outfile {prosite_path2}/pep_stats_results/{fasta_file_name_only}.pepstats"
#Defining the vairbales to be used in pepstats. This is another wildcard option showing basic info about the proteins as I thought the user might be interested in it. 

os.system(pepstats_input)
print('A pepstats file has been created in the pep_stats_results directory showing protein statistis such as molecular weight, number of residues, etc.\nNavigate to this directory to view results.')
#Running pepstats and telling the user where to find the results
