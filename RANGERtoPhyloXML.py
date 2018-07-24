import sys, xml.etree.cElementTree as ET

inputFile = sys.argv[1]

def findRec(lines) :
    start = 0
    end = 0
    for lineNum, current in enumerate(lines) :
        if (current.find("Reconciliation:") != -1) :
           start = lineNum+1
        if (not current.strip() and (start != 0)) :
            end = lineNum
            return (start,end)

def buildXML() :
    root = ET.Element("recGeneTree")
    rooted = ET.SubElement(root, "phylogeny", rooted="true")

    #hardcoded example
    clade = ET.SubElement(rooted, "clade")
    name = ET.SubElement(clade, "name").text = "m3"


    tree = ET.ElementTree(root)
    tree.write("file.xml")
        
with open(inputFile, 'r') as file :
    xmlLines = []
    lines = file.readlines()
    s, e = findRec(lines)[0], findRec(lines)[1]
    recLines = lines[s:e]
    buildXML()

    
