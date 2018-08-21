import sys
import argparse
from collections import defaultdict
from Bio import Phylo
import os


#Finds where Reconcliation starts and ends
#This is the section underneath the "Reconcilation:" and
#right above the whitespace after the last event
def findRec(lines) :
    start = 0
    end = 0
    for lineNum, current in enumerate(lines) :
        if (current.find("Reconciliation:") != -1) :
           start = lineNum+1
        if (not current.strip() and (start != 0)) :
            end = lineNum
            return (start,end)

#Takes the newick form of the trees from the file and returns them
def findSpTree(lines) :
    for lineNum, current in enumerate(lines) :
        if(current.find("Species Tree:") != -1) :
            return lines[lineNum+1]

def findGeneTree(lines) :
    for lineNum, current in enumerate(lines) :
        if(current.find("Gene Tree:") != -1) :
            return lines[lineNum+1]

#All XML generators found here
#Extracts data from input lines and generates the appropriate XML
def transferXML(line,genetree) :
    
    #Data Extraction
    subNode = line[0:line.find(' =')]
    mapper = line[line.find('Mapping') + 12:line.find(', Recipient')]
    recip = line[line.find(', Recipient') + 16:len(line)]
    lca = line[line.find('[')+1:line.find(']:')].split(', ')
    subString = "<name>" + subNode + "</name>"
    for num, line in enumerate(genetree):
        if line.find(subString) != -1:
            leadingTabs = len(line) - len(line.lstrip())
            tabs = '\t'*(int(leadingTabs/2)+1)
            newLine = tabs + '<eventsRec>\n' + tabs + '\t<branchingOut speciesLocation=' + '\"' + mapper.rstrip() + '\"' + '></branchingOut>\n' + tabs + '</eventsRec>\n'
            print(leadingTabs)
            genetree.insert(num+1, newLine)
        
    return genetree

def duplicationXML(line,genetree) :

    #Data Extraction
    subNode = line[0:line.find(' =')]
    mapper = line[line.find('Mapping') + 12:len(line)]
    lca = line[line.find('[')+1:line.find(']:')].split(', ')
    subString = "<name>" + subNode + "</name>"
    for num, line in enumerate(genetree):
        if line.find(subString) != -1:
            leadingTabs = len(line) - len(line.lstrip())
            tabs = '\t'*(int(leadingTabs/2)+1)
            newLine = tabs + '<eventsRec>\n' + tabs + '\t<duplication speciesLocation=' + '\"' + mapper.rstrip() + '\"' + '></duplication>\n' + tabs + '</eventsRec>\n'
            print(leadingTabs)
            genetree.insert(num+1, newLine)
        
    return genetree

def speciationXML(line, genetree) :
    
    #Data Extraction
    subNode = line[0:line.find(' =')]
    mapper = line[line.find('Mapping') + 12:len(line)]
    lca = line[line.find('[')+1:line.find(']:')].split(', ')
    #insert into the genetree
    subString = "<name>" + subNode + "</name>"
    for num, line in enumerate(genetree):
        if line.find(subString) != -1:
            leadingTabs = len(line) - len(line.lstrip())
            tabs = '\t'*(int(leadingTabs/2)+1)
            newLine = tabs + '<eventsRec>\n' + tabs + '\t<speciation speciesLocation=' + '\"' + mapper.rstrip() + '\"' + '></speciation>\n' + tabs + '</eventsRec>\n'
            print(leadingTabs)
            genetree.insert(num+1, newLine)
        
    return genetree

def leafXML(line,genetree) :
    leafName = line[0:line.find(': ')]
    subString = "<name>" + leafName + "</name>"
    for num, line in enumerate(genetree):
        if line.find(subString) != -1:
            leadingTabs = len(line) - len(line.lstrip())
            tabs = '\t'*(int(leadingTabs/2)+1)
            newLine = tabs + '<eventsRec>\n' + tabs + '\t<leaf speciesLocation=' + '\"' + leafName.rstrip() + '\"' + '></leaf>\n' + tabs + '</eventsRec>\n'
            print(leadingTabs)
            genetree.insert(num+1, newLine)
        
    return genetree

#Creates XML tree using biopython based on inputted newick tree
def buildTree(tree, qualifier) :
    with open("temp","w+") as temper :
        temper.writelines(tree)
    Phylo.convert('temp', 'newick', 'temp2', 'phyloxml')
    
    specFile = open("temp2", "r")
    lines = specFile.readlines()
    lines.pop(0) #remove first line
    lines.pop() #remove last line
    specFile.close()
    os.remove("temp")
    os.remove("temp2")

    #Inserting the recPhylo information to the gene/species tree
    for line in lines:
        line = "\t" + line
    if(qualifier == "s") :
        lines.insert(0,"<recPhylo>\n\t<spTree>\n")
        lines.append("\t</spTree>\n")
    elif(qualifier == "g") :
        lines.insert(0,"\t<recGeneTree>\n")
        lines.append("\t</recGeneTree>\n</recPhylo>\n")
    
    return lines

#This takes the locations of each event and creates the appropriate XML
def buildXML(recLines,geneTree) :
    events = ("Transfer", "Duplication", "Speciation")
    for line in recLines :
        if events[0] in line :
            #Transfer XML
            geneTree = transferXML(line,geneTree)

        elif events[1] in line :
            #Duplication XML
            geneTree = duplicationXML(line,geneTree)

        elif events[2] in line :
            #Speciation XML
            geneTree = speciationXML(line,geneTree)

        else :
            geneTree = leafXML(line,geneTree)


#Argument parser, checks input from command line, arguments are optional
parser = argparse.ArgumentParser(description='Commands')
parser.add_argument('-i', '--input', help="Input file path", default = sys.argv[1])
parser.add_argument('-o', '--output', help="Output file name", default = 'Output')
args = parser.parse_args()
outputFile = args.outputFile

#Hard coded input file for when biopython doesn't want to work
#inputFile = r'C:\Win\CorePrograms\test'


#Start of program, takes input and runs!
with open(args.input,'r') as file :
    xmlLines = []
    lines = file.readlines()
    s, e = findRec(lines)[0], findRec(lines)[1]
    recLines = lines[s:e]
    spTreeNewick = findSpTree(lines)
    geneTreeNewick = findGeneTree(lines)
    spTree = buildTree(spTreeNewick, "s")
    geneTree = buildTree(geneTreeNewick, "g")
    buildXML(recLines,geneTree)


    #If args.output doesn't work, change to any string file for the output
    with open(args.output, "w+") as outFile:
        outFile.writelines(spTree)
        outFile.writelines(geneTree)
