#!/usr/bin/env python3

def rawNetMHCToTable(rawData):
    import re
    data = re.sub("\n\-{15,}\n{1,2}.+\n{1,2}\-{15,}", "", rawData)
    data = data.split("\n")
    validLines = []
    for line in data:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            continue
        validLines.append(line)
    return validLines

def getPredictionTable(rawPredictionFile):
    import variantDataHandler
    file = open(rawPredictionFile)
    rawData = file.read()
    file.close()
    rawTable = rawNetMHCToTable(rawData)
    predictionTable = {}
    for line in rawTable:
        try:
            prediction = variantDataHandler.NetMHCPrediction(line)
        except RuntimeError:
            print(rawPredictionFile)
            quit()
        if not prediction.peptide in predictionTable:
            predictionTable[prediction.peptide] = {}
        predictionTable[prediction.peptide][prediction.hla] = prediction
    return predictionTable

if __name__ == '__main__':
    data = getPredictionTable("predictions.txt")
    print("something")