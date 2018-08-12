import sys, xml.etree.cElementTree as ET
import argparse
from collections import defaultdict

def is_valid_file(parser,arg):
    return open(arg,'r')

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

#All XML generators found here
#Extracts data from input lines and generates the appropriate XML
def transferXML(line, root) :
    
    #Data Extraction
    subNode = line[0:line.find(' =')]
    mapper = line[line.find('Mapping') + 12:line.find(', Recipient')]
    recip = line[line.find(', Recipient') + 16:len(line)]
    lca = line[line.find('[')+1:line.find(']:')].split(', ')

    #Creating XML
    #next = ET.SubElement(root, "clade")
    #next = ET.SubElement(next, "transfer" + lca, mapper + " " + recip)
    return(next)

def duplicationXML(line, prev) :

    #Data Extraction
    subNode = line[0:line.find(' =')]
    mapper = line[line.find('Mapping') + 12:len(line)]
    lca = line[line.find('[')+1:line.find(']:')].split(', ')

    #Creating XML
    next = prev
    return(next)


def speciationXML(line, prev) :
    
    #Data Extraction
    subNode = line[0:line.find(' =')]
    mapper = line[line.find('Mapping') + 12:len(line)]
    lca = line[line.find('[')+1:line.find(']:')].split(', ')
    
    #Creating XML
    next = prev
    return(next)

def leafXML(line, prev) :
    subNode = line[0:line.find(': ')]
    next = prev
    return(next)

#This takes the locations of each event and creates the appropriate XML
def buildXML(recLines) :
    root = ET.Element("recGeneTree")
    rooted = ET.SubElement(root, "phylogeny", rooted="true")
    events = ("Transfer", "Duplication", "Speciation")
    prev = root
    for line in recLines :
        #if any(x in line for x in events) :
        if events[0] in line :
            #Transfer XML
            next = transferXML(line,prev)
            #print(events[0])
        elif events[1] in line :
            #Duplication XML
            next = duplicationXML(line,prev)
            #print(events[1])
        elif events[2] in line :
            #Speciation XML
            next = speciationXML(line,prev)
            #print(events[2])
        else :
            next = leafXML(line,prev)
        prev = next
    #hardcoded example
    clade = ET.SubElement(rooted, "clade")
    name = ET.SubElement(clade, "name").text = "m3"

    tree = ET.ElementTree(root)
    tree.write(outputFile)

#Argument parser, checks input from command line       
parser = argparse.ArgumentParser(description='Commands')
parser.add_argument('-i', '--input', help="Input file path", default = sys.argv[1])
parser.add_argument('-o', '--output', help="Output file name", default = 'Output')
args = parser.parse_args()
outputFile = args.output

#Start of program, takes input and runs!
with open(args.input,'r') as file :
    xmlLines = []
    lines = file.readlines()
    s, e = findRec(lines)[0], findRec(lines)[1]
    recLines = lines[s:e]
    buildXML(recLines)