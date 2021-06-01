import pubmed_lookup as pl
import os
import json

def pidFromName(name):
    startIdx = 10000000
    maxLength = 0
    maxIdx = 0
    for i,c in enumerate(name):
        if c < '0' or c > '9':
            startIdx = 10000000
            continue
        if startIdx > i:
            startIdx = i
        if i-startIdx+1 > maxLength:
            maxIdx = startIdx
            maxLength = i-startIdx+1
    
    return name[maxIdx:(maxIdx+maxLength)]

def buildDict(folders, outputPath):
    outputDict = {}
    for folder in folders:
        for fname in os.listdir(folder):
            if not fname[-4:] == ".txt":
                continue
            pid = pidFromName(fname)

            try:
                lookup = pl.PubMedLookup('http://www.ncbi.nlm.nih.gov/pubmed/'+pid, '')
            except:
                print("could not convert ", fname, " because a valid pid was not found in title (pid assumed to be the longest sequence of digits in the name)")
                continue

            publication = pl.Publication(lookup)
            title = publication.title
            outputDict[pid] = title
            print("converted ", fname)
    
    outputFile = open(outputPath, "w")
    outputFile.write(json.dumps(outputDict))
    outputFile.close()


buildDict(["update_rules_cell_collective", "update_rules_models_in_literature_we_randomly_come_across", "update_rules_not_to_load"], "pidToTitle.txt")



