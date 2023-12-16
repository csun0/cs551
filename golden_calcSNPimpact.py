import argparse
import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def readCoupMat(inputCouplingMatrix):
    matrix = np.loadtxt(open(inputCouplingMatrix, "rb"), delimiter=',')
    return matrix

def readMSA(inputMSAFile):
    emptyPosPos = {}
    psuedoCount = 1/441
    jointFreq = {}
    aminoAcids = "ACDEFGHIKLMNPQRSTVWY-"
    for aa1 in aminoAcids:
        for aa2 in aminoAcids:
            emptyPosPos[(aa1, aa2)] = 0

    def initializeJointFreq(sequenceLength):
        for i in range(1,sequenceLength+1):
            for j in range(i+1,sequenceLength+1):
                jointFreq[(i,j)] = copy.deepcopy(emptyPosPos)


    def updateJointFreq(sequence):
        for i, aa1 in enumerate(sequence):
            if(aa1 not in aminoAcids):
                aa1 = '-'
            for j, aa2 in enumerate(sequence):
                if(aa2 not in aminoAcids):
                    aa2 = '-'
                if(j <= i):
                    continue
                jointFreq[(i+1,j+1)][(aa1,aa2)] += 1


    seq = ""
    entryCount = 0
    with open(inputMSAFile,'r') as msa:
        for line in msa:
            if(line[0] == '>'):
                entryCount += 1
                # First instace of >
                if(seq == ""):
                    continue
                # Second instance of >
                if(not bool(jointFreq)):
                    initializeJointFreq(len(seq))
                # Future instances of >
                updateJointFreq(seq)
                seq = ""
            else:
                seq += line.strip()

    for PosPosCount in jointFreq.values():
        for AApair in PosPosCount.keys():
            PosPosCount[AApair] = (PosPosCount[AApair]/entryCount) + psuedoCount

    return jointFreq

def normAAWeights(aaFreqAtPos):
    norm_weights = np.array(aaFreqAtPos)
    norm_weights[norm_weights < 0] = 0
    norm_weights = norm_weights/sum(norm_weights)
    return norm_weights

def calcMutantCoVchange(weightArr, jointFreq, mutation, position, nativeSeq):
    scoreSum = 0
    initRes = nativeSeq[position-1]
    ppm_joint_native = []
    ppm_joint_mutant = []
    aaList = "ACDEFGHIKLMNPQRSTVYW"
    for index, weight in enumerate(weightArr[position-1]):
        j = index + 1
        if(position == j):
            continue
        if(position < j):
            posPosMap = jointFreq[(position, j)]
            print(position)
            print(j)
            probNative = posPosMap[initRes, nativeSeq[j-1]]
            ppm_joint_native.append(normAAWeights([posPosMap[initRes, a] for a in aaList]))

            probMutant = posPosMap[mutation, nativeSeq[j-1]]
            ppm_joint_mutant.append(normAAWeights([posPosMap[mutation, a] for a in aaList]))
        else:
            posPosMap = jointFreq[(j, position)]
            probNative = posPosMap[nativeSeq[j-1], initRes]
            ppm_joint_native.append(normAAWeights([posPosMap[a, initRes] for a in aaList]))

            probMutant = posPosMap[nativeSeq[j-1], mutation]
            ppm_joint_mutant.append(normAAWeights([posPosMap[a, mutation] for a in aaList]))

        scoreSum += np.abs(weight)*np.log2(probMutant/probNative)
    ppm_native = np.transpose(np.asarray(ppm_joint_native))
    ppm_mutant = np.transpose(np.asarray(ppm_joint_mutant))
    return scoreSum

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(prog='Analyze MSA freq', description='Convert an MSA to a ij freq table and analyze impact of mutation')
    argparser.add_argument('-i', '--inputMSA', type=str, help='MSA input file', required=True)
    argparser.add_argument('-w', '--weightMatrix', type=str, help='Input paramters coupling scores matrix from plmc', required=True)
    argparser.add_argument('-p', '--position', type=str,
                           help='Relevant position to examine indexed from MSA', required=True)
    argparser.add_argument('-m', '--mutation', type=str, help='Mutation in one letter code of mutant residue')
    args = argparser.parse_args()
    nativeSeq = '----------------------------------------------------------------------DLMSELQKDSIQLDEDS---ERKVVKMLLRLLEDKNGEVQNLAVKCLGPLVVKVKEYQVETIVDTL'

    weightsMat = readCoupMat(args.weightMatrix)
    jointFreq = readMSA(args.inputMSA)
    aaList = "ACDEFGHIKLMNPQRSTVYW"
    matplotlib.rcParams.update({'font.size': 40})

    scoresPerAA = []
    for a in aaList:
        score = calcMutantCoVchange(weightsMat, jointFreq, a, int(args.position), nativeSeq)
        if(a == args.mutation):
            print(score)
        scoresPerAA.append(score)
    fig, ax = plt.subplots()
    bar_chart = ax.bar(list(aaList), scoresPerAA)
    plt.show()
    #ppm_native_obj = seqlogo.Ppm(ppm_native, alphabet_type='AA')
    #ppm_mutant_obj = seqlogo.Ppm(ppm_mutant, alphabet_type='AA')
    import os
    #seqlogo.seqlogo(ppm_native_obj, ic_scale=True, format='png', size='xlarge', filename=os.path.join('params_out','native_seqlogo.png'))
    #seqlogo.seqlogo(ppm_mutant_obj, ic_scale=True, format='png', size='xlarge', filename=os.path.join('params_out','mutant_seqlogo.png'))

