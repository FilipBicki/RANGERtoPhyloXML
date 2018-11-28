#########################################
##  Author:         Filip Bicki
##  Created:        20-Jun-2018
##  Last modified:  13-Nov-2018
##
##  What is does:
##  Takes a RANGER-DTL file output and converts it
##  into the recPhyloXML format.
##
##  requires : biopython (https://biopython.org/)
##
##
##  developped for python3.7
##
#########################################

import sys
import argparse
from collections import defaultdict
from Bio import Phylo
import os
import newick

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

#Checks to see if the tree was rooted (trivial?)
def findRooted(lines) :
    for current in lines :
        if(current.find("(rooted)") != -1) :
            return True
    return False

#Checks if a gene on the gene tree has had an event already added
def eventsRec(name,genetree) :
    for num, line in enumerate(genetree) :
        if(line.find(name) != -1) :
            if(genetree[num+1].strip() != "</clade>" and genetree[num +1].strip() != "<clade>") :
                return True
    return False
    

#All XML generators found here
#Extracts data from input lines and generates the appropriate XML
def transferXML(line,genetree,reclines) :
    
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
            #If the gene already has an <eventsRec> tag, don't add another one
            if eventsRec(subNode,genetree) :
                newLine = tabs + '\t<branchingOut speciesLocation=' + '\"' + mapper.rstrip() + '\"' + '></branchingOut>\n'
                genetree.insert(num+3, newLine)
            else :
                genetree.insert(num+1, tabs + '<eventsRec>\n')
                genetree.insert(num+2, tabs + '\t<branchingOut speciesLocation=' + '\"' + mapper.rstrip() + '\"' + '></branchingOut>\n')
                genetree.insert(num+3, tabs + '</eventsRec>\n')
    #Transfer events require a second part, the transferback, done here
    
    #find children of current gene node
    #find if genes are
    #leaves -> string matching
    #direct mapping -> pick that one
    #no mapping -> look for all children

    for line in reclines:
        print(lca[1])
        print(recip)
        if (lca[0] + "= " or lca[1] + "= " in line) and "Mapping --> " + recip in line:
            genetree = transferBackXML(line,genetree)
        elif lca[0].find(recip) != -1 :
            genetree = transferBackXML(line,genetree)
        elif lca[1].find(recip) != -1 :
            genetree = transferBackXML(line,genetree)

        #find other ways to get recipient species
    return genetree

def transferBackXML(line,genetree) :

    subNode = line[0:line.find(' =')]
    if line.find("Transfer") :
        mapper = line[line.find('Mapping') + 12:line.find(', Recipient')]
    else :
        mapper = line[line.find('Mapping') + 12:len(line)]
    lca = line[line.find('[')+1:line.find(']:')].split(', ')
    subString = "<name>" + subNode + "</name>"
    for num, line in enumerate(genetree):
        if line.find(subString) != -1:
            leadingTabs = len(line) - len(line.lstrip())
            tabs = '\t'*(int(leadingTabs/2)+1)
            #If the gene already has an <eventsRec> tag, don't add another one
            if eventsRec(subNode,genetree) :
                newLine = tabs + '\t<transferBack destinationSpecies=' + '\"' + mapper.rstrip() + '\"' + '></transferBack>\n'
                genetree.insert(num+2, newLine)
            else :
                genetree.insert(num+1, tabs + '<eventsRec>\n')
                genetree.insert(num+2, tabs + '\t<transferBack destinationSpecies=' + '\"' + mapper.rstrip() + '\"' + '></transferBack>\n')
                genetree.insert(num+3, tabs + '</eventsRec>\n')
        
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
            #If the gene already has an <eventsRec> tag, don't add another one
            if eventsRec(subNode,genetree) :
                newLine = tabs + '\t<duplication speciesLocation=' + '\"' + mapper.rstrip() + '\"' + '></duplication>\n'
                genetree.insert(num+3, newLine)
            else :
                genetree.insert(num+1, tabs + '<eventsRec>\n')
                genetree.insert(num+2, tabs + '\t<duplication speciesLocation=' + '\"' + mapper.rstrip() + '\"' + '></duplication>\n')
                genetree.insert(num+3, tabs + '</eventsRec>\n')
        
    return genetree

def speciationXML(line, genetree) :
    
    #Data Extraction
    subNode = line[0:line.find(' =')]
    mapper = line[line.find('Mapping') + 12:len(line)]
    lca = line[line.find('[')+1:line.find(']:')].split(', ')
    subString = "<name>" + subNode + "</name>"
    for num, line in enumerate(genetree):
        if line.find(subString) != -1:
            leadingTabs = len(line) - len(line.lstrip())
            tabs = '\t'*(int(leadingTabs/2)+1)
            #If the gene already has an <eventsRec> tag, don't add another one
            if eventsRec(subString,genetree) :
                newLine = tabs + '\t<speciation speciesLocation=' + '\"' + mapper.rstrip() + '\"' + '></speciation>\n'
                genetree.insert(num+3, newLine)
            else :
                genetree.insert(num+1, tabs + '<eventsRec>\n')
                genetree.insert(num+2, tabs + '\t<speciation speciesLocation=' + '\"' + mapper.rstrip() + '\"' + '></speciation>\n')
                genetree.insert(num+3, tabs + '</eventsRec>\n')
        
    return genetree

def leafXML(line,genetree) :
    leafName = line[0:line.find(': ')]
    subString = "<name>" + leafName + "</name>"
    for num, line in enumerate(genetree):
        if line.find(subString) != -1:
            leadingTabs = len(line) - len(line.lstrip())
            tabs = '\t'*(int(leadingTabs/2)+1)
            #If the gene already has an <eventsRec> tag, don't add another one
            if eventsRec(subString,genetree) :
                newLine = tabs + '\t<leaf speciesLocation=' + '\"' + leafName.rstrip() + '\"' + '></leaf>\n'
                genetree.insert(num+3, newLine)
            else :
                genetree.insert(num+1, tabs + '<eventsRec>\n')
                genetree.insert(num+2, tabs + '\t<leaf speciesLocation=' + '\"' + leafName.rstrip() + '\"' + '></leaf>\n')
                genetree.insert(num+3, tabs + '</eventsRec>\n')
        
    return genetree

#Creates XML tree using biopython based on inputted newick tree
def buildTree(tree, qualifier, rooted) :
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
    if rooted :     
        lines.pop(0)
        lines.insert(0,'\t<phylogeny rooted="true">\n')

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
            print("Transfer")
            geneTree = transferXML(line,geneTree,recLines)

        elif events[1] in line :
            #Duplication XML
            #print("Duplication")
            geneTree = duplicationXML(line,geneTree)

        elif events[2] in line :
            #Speciation XML
            #print("Speciation")
            geneTree = speciationXML(line,geneTree)

        else :
            #print("Leaf")
            geneTree = leafXML(line,geneTree)

    return geneTree


#Argument parser, checks input from command line, arguments are optional
parser = argparse.ArgumentParser(description='Commands')
parser.add_argument('-i', '--input', help="Input file path", default = sys.argv[1])
parser.add_argument('-o', '--output', help="Output file name", default = 'output.xml')
args = parser.parse_args()
#outputFile = args.outputFile


#Start of program, takes input and runs!
with open(args.input,'r') as file :
    xmlLines = []
    lines = file.readlines()
    s, e = findRec(lines)[0], findRec(lines)[1]
    recLines = lines[s:e]
    rooted = findRooted(lines)
    spTree = buildTree(findSpTree(lines), "s", rooted)
    geneTree = buildTree(findGeneTree(lines), "g", rooted)
    geneTree = buildXML(recLines,geneTree)

    #Combining into one list
    spTree.extend(geneTree)

    #Fixing the format
    for line in spTree :
        line.rstrip()
        line.strip()

    #Writing the output file based on the arguments/default
    with open(args.output, "w+") as outFile:
        outFile.writelines(spTree)
        
