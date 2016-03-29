import os
import sys
import optparse
import re        # regular expressions
import string

from tt_log import logger

import numpy as np

import Annotations as anno
import Best        as best
import Cluster     as cl
import CigarString as cs
import pandas as pd
from sklearn.cluster import KMeans

VERSION = '20150529.01'

DEF_OUTPUT = 'exons.png'        # default plot filename
DEF_YSCALE = 1.0                # default  Y-axis scale factor

FIG_WIDTH = 14
FIG_HEIGHT_PER_TRANS = 0.2      # figure height depends on the number of rows
MAX_LABELS = 20                 # how many labels fit across the X axis
EXTRA_SPACE = 1.0               # added vertical space in figure for labels, etc

FASTA_WRAP = 60                 # bases per fasta line
Q_THRESHOLD = 20.0              # passing grade for Q score
MIN_REGION_SIZE = 50

REGEX_NAME = re.compile ('(c\d+)')      # cluster ID in cluster name
REGEX_LEN  = re.compile ('\/(\d+)$')     # cluster length in cluster name
COMPLTAB   = string.maketrans ('ACGTacgt', 'TGCAtgca')     # for reverse-complementing reads

def getGeneFromAnnotation (opt, tranList, exonList):
    '''Add to lists of transcripts and exons: annotations for gene of interest.'''

    if opt.gtf == None:
        return tranList, exonList

    omits = [] if opt.omit is None else opt.omit.split(',')            # transcripts which must not be included

    if opt.format == 'pickle':
        annotList   = anno.AnnotationList.fromPickle (opt.gtf)
    elif opt.format == 'alt':
        annotList   = anno.AnnotationList (opt.gtf, altFormat=True)
    else:     # standard format
        annotList   = anno.AnnotationList (opt.gtf)

    allGenes = annotList.getGeneDict()
    if opt.gene not in allGenes:
        raise RuntimeError ('gene %s is not in the annotation file' % opt.gene)
    geneList = allGenes[opt.gene]       # a list of Annotation objects
    if len(geneList) > 1:
        logger.warning('gene %s appears %d times in annotations, first occurrence plotted' \
                           % (opt.gene, len(geneList)))
    myGene = geneList[0]

    for tran in myGene.getChildren():                       # tran is an Annotation object

        if tran.name not in omits:                          # if not in ignore list

            myTran = Transcript(tran.name, annot=True)

            if hasattr(tran, 'startcodon'):
                myTran.startcodon = tran.startcodon
            if hasattr(tran, 'stopcodon'):
                myTran.stopcodon = tran.stopcodon

            for exon in tran.getChildren():                 # exon is an Annotation object
                myExon = Exon(myTran, exon.name, exon.start, exon.end, exon.strand)     # no Q score
                if hasattr (exon, 'polyAs'):
                    print exon.name
                    myExon.polyAs = exon.polyAs
                exonList.append (myExon)
                myTran.exons.append(myExon)

            tranList.append (myTran)

    return tranList, exonList

def getGeneFromMatches (opt, tranList, exonList):
    '''Add to lists of transcripts and exons: clusters which matched gene of interest.'''

    if opt.matches == None:
        return tranList, exonList

    omits = [] if opt.omit is None else opt.omit.split(',')            # clusters which must not be included
    shows = [] if opt.show is None else opt.show.split(',')            # clusters which must be included

    fullThreshold = -1
    partialThreshold = -1
    if opt.highsupport:                                                          # trim clusters with low support
        if opt.full is not None:
            fullThreshold = opt.full               # clusters should has higher full support than input
        if opt.partial is not None:
            partialThreshold =  opt.partial      # clusters should has higher partial support than input


    localList = list()                                                 # temporary list of clusters
    totClusters = 0

    for matchFile in opt.matches:                                      # --matches may have been specified more thn once

        clusterDict = cl.ClusterDict.fromPickle (matchFile)            # pickle file produced by matchAnnot.py

        for cluster in clusterDict.getClustersForGene(opt.gene):       # cluster is Cluster object

            totClusters += 1

            match = re.search (REGEX_NAME, cluster.name)

            if match is not None and match.group(1) in shows:          # if this is a force-include cluster
                localList.append ( [cluster, 'f999999p999999'] )       # fake sort key to push it to the front

            elif match is None or match.group(1) not in omits:         # shows and omits trump length filter

                full, partial = cluster.getFP()
                sortKey = 'f%06dp%06d' % (full, partial)               # single key includes full and partial counts

                matchLen = re.search(REGEX_LEN, cluster.name)          # filter by cluster length, if requested


                if full >= fullThreshold and partial >= partialThreshold:
                    if matchLen is None:
                        raise RuntimeError ('no length in name: %s' % cluster.name)
                        localList.append ( [cluster, sortKey] )            # shouldn't happen -- but let it slide
                    else:
                        cLen = int(matchLen.group(1))
                        if opt.minlen is None or cLen >= opt.minlen:
                            if opt.maxlen is None or cLen <= opt.maxlen:
                                localList.append ( [cluster, sortKey] )

    localList.sort(key=lambda x: x[1], reverse=True)                   # sort by full/partial counts

    if opt.nodups:                                                     # eliminate exact dups?

        tempList = list()
        uniqueClusters = set()
        totDups = 0

        for ent in localList:
            cluster = ent[0]
            key = '%9d %s | %s' % (cluster.start, cluster.cigar.prettyPrint(), cluster.cigar.MD)
            if key in uniqueClusters:
                totDups += 1
            else:
                tempList.append(ent)                                   # keep this
                uniqueClusters.add(key)                                # remember it

        localList = tempList
        logger.debug('discarded %d clusters as exact duplicates of each other' % totDups)

    if opt.howmany is not None:
        localList = localList[:opt.howmany]                            # keep the top N entries (which will include the forces)

    totFull = 0
    totPartial = 0

    for ent in localList:

        cluster = ent[0]
        myTran = Transcript(cluster.name, score=cluster.bestScore)

        full, partial = cluster.getFP()
        totFull += full
        totPartial += partial

        leading, trailing = cluster.cigar.softclips()

        for exonNum, exon in enumerate(cluster.cigar.exons()):         # exon is a cs.Exon object

            exonName = '%s/%d' % (myTran.name, exonNum)                # exons don't have names: make one up

            if cluster.cigar.MD is not None:                           # if MD string was supplied
                myExon = Exon(myTran, exonName, exon.start, exon.end, cluster.strand, QScore=exon.QScore())
            else:
                myExon = Exon(myTran, exonName, exon.start, exon.end, cluster.strand)

            if exonNum == 0:
                myExon.leading = leading              # add leading softclips to first exon

            exonList.append (myExon)
            myTran.exons.append(myExon)

        myExon.trailing = trailing                    # add trailing softclips to last exon

        tranList.append (myTran)

        if opt.fasta is not None:
            writeFasta (opt, cluster)

    logger.debug('kept %d of %d clusters for gene %s' % (len(localList), totClusters, opt.gene))
    logger.debug('kept clusters include %d full + %d partial reads' % (totFull, totPartial))

    return tranList, exonList


def assignBlocks (opt, exonList):
    '''
    Assign exons to blocks, separated by sequence which is intronic in
    all transcripts. exonList is assumed to be sorted by ascending
    start position.
    '''

    adjust  = 0
    blockNo = 0
    exonIx  = 0
    blocks = list()

    while exonIx < len(exonList):               # can't use enumerate here, it's a double loop

        blockStart = exonList[exonIx].start     # block start = start of first exon in block
        blockEnd   = exonList[exonIx].end       # initial value, updated in the loop below
        blockStartIx = exonIx

        while exonIx < len(exonList) and exonList[exonIx].start <= blockEnd:
            myExon = exonList[exonIx]
            if myExon.end > blockEnd:
                blockEnd = myExon.end
            myExon.block = blockNo
            myExon.tran.blocks.add(blockNo)    # transcript has an exon in this block
            myExon.adjStart = myExon.start - blockStart + adjust
            exonIx += 1

        adjust += blockEnd - blockStart + 1
        blocks.append(Block(blockStart, blockEnd, adjust))
        blockNo += 1

    return blocks

def assignBlocksReverse (opt, exonList):
    '''
    Like assignblocks, but for the reverse strand, ordering blocks
    from the 5' end of the transcript. exonList is assumed to be
    sorted by decreasing exon end position.
    '''

    # I did this as a separate mirror image of assignBlocks, rather
    # than clutter the scenery with lots of forward/reverse checks.

    adjust  = 0
    blockNo = 0
    exonIx  = 0
    blocks = list()

    while exonIx < len(exonList):               # can't use enumerate here, it's a double loop

        blockStart = exonList[exonIx].end       # block start = end of last exon in block
        blockEnd   = exonList[exonIx].start     # initial value, updated in the loop below
        blockStartIx = exonIx

        while exonIx < len(exonList) and exonList[exonIx].end >= blockEnd:

            myExon = exonList[exonIx]
            if myExon.start < blockEnd:
                blockEnd = myExon.start
            myExon.block = blockNo
            myExon.tran.blocks.add(blockNo)    # transcript has an exon in this block
            myExon.adjStart = blockStart - myExon.end + adjust
            exonIx += 1

        adjust += blockStart - blockEnd + 1
        blocks.append(Block(blockStart, blockEnd, adjust))
        blockNo += 1

    return blocks

def findRegions (tranList):
    '''Find breakpoints where coverage by exons changes.'''

    # Why are we doing this? See the note in the Transcript class
    # definition below.

    breaks = list()

    for tranIx, tran in enumerate(tranList):

        for exon in tran.exons:
            breaks.append ([exon.start, 0, tranIx, tran.name, exon.name])
            breaks.append ([exon.end,   1, tranIx, tran.name, exon.name])

    breaks.sort (key=lambda x: x[0])
    curPos = breaks[0][0]
    curTranSet = set()
    region = 0

    for ix in xrange(len(breaks)):

        posit, flag, tranIx, tranName, exonName = breaks[ix]

        if posit > curPos + MIN_REGION_SIZE:             # this is a new region
            if len(curTranSet) > 0:
                for ix in curTranSet:
                    tranList[ix].regions.add(region)     # update set of regions hit by this transcript
                region += 1
            curPos = posit

        if flag == 0:                                    # exon start
####            print '%9d  start  %s' % (posit, exonName)
            curTranSet.add (tranIx)
        else:                                            # exon end
####            print '%9d  end    %s' % (posit, exonName)
            curTranSet.remove (tranIx)

    logger.debug('found %d regions' % region)

    return

def orderTranscripts (tranList):
    '''
    Order the transcripts (i,e., assign each a Y coordinate) so similar
    transcripts are close to each other.
    '''

    # The measure of similarity used here is block occupancy: The
    # distance between two transcripts is the number of blocks where
    # one transcript has exons, and the other doesn't. How many exons
    # there are, or how similar they are in length, is not looked at.

    # The ordering is done using a greedy nearest-neighbor
    # heuristic. To do it optimally turns it into a Traveling Salesman
    # problem.

    tranNames = list()
    curTran = tranList[0]            # arbitrarily start with the first transcript
    tranIx = 0

    while True:                                     # loop until break below

        tranNames.append(curTran.name)              # needed for yticks call
        curTran.tranIx = tranIx

        bestTran = best.Best(reverse=True)
        for myTran in tranList:                     # find the next closest transcript
            if myTran.tranIx is None:               # if transcript hasn't been indexed yet
                diff = len(curTran.regions.symmetric_difference(myTran.regions))
                bestTran.update(diff, myTran)

        if bestTran.which is None:                  # every transcript has its index: we're done
            break
####        else:
####            logger.debug('%2d  %s' % (bestTran.value, bestTran.which.name))

        curTran = bestTran.which
        tranIx += 1

    return tranNames

def groupTran(tranList, exonList, cluster_num):
    df = pd.DataFrame()
    Name = list()
    # Get the names of all distinct transcript
    for tran in tranList:
        if '-' not in tran.name:
            Name.append(tran.name)
    # Append distinct name to the dataframe
    df['name'] = Name
    length = len(df)

    # minVal is the minimum starting point of all the transcripts, maxVal stands for maximum
    AllExon = list()
    minVal = float("inf")
    maxVal = 0

    # Need improvement!!! # for each distint transcript, append [start, end] of all its exons
    for name in Name:
        tranExons = list()
        for exon in exonList:
            if name in exon.name:
                start = exon.start
                end = exon.end
                if start > maxVal or end > maxVal:
                    maxVal = max(start, end)
                if start < minVal or end < minVal:
                    minVal = min(start, end)
                tranExons.append([start, end])
        AllExon.append(tranExons)

    # Set the minimum starting point to 0
    for i in AllExon:
        for j in i:
            j[0] = j[0] - minVal
            j[1] = j[1] - minVal
    df['exon'] = AllExon

    # Build a matrix contains only true and false
    #
    #   Transcript1:    -----    ----  -- -------
    #   Transcript2:  ----  ------  ------- ---
    #   Superimpose:  ---------------------------
    #   booleanTran:  FFFFFFFFFFFFFFFFFFFFFFFFFFF
    #   Overlap1:     FFTTTTTFFFFTTTTFFTTFTTTTTTT

    #     booleanTran has the same length with superimpose, which is a list of False. Change the overlap region
    #   between the booleanTran and each transcript to True.

    booleanMatrix = list()
    for tranExon in df.exon:
        booleanTran = [False for x in range(maxVal-minVal)]
        for exon in tranExon:
            booleanTran[exon[0]:exon[1]+1] = [True for x in range(exon[1]+1-exon[0])]
        booleanMatrix.append(booleanTran)
    df['boolean'] = booleanMatrix

    # Create a distance table that can be used in K-Means.
    #
    #                        c225/f26p50/6117  c483/f8p23/6083  c20615/f3p27/6185
    # c225/f26p50/6117           0.000000         0.029911           0.012941
    # c483/f8p23/6083            0.029911         0.000000           0.019231
    # c20615/f3p27/6185          0.012941         0.019231           0.000000

    #    The number can be interpreted as the similarity between each two transcript. 0 means they are exactly same
    #  while 1 means they have no overlap region.

    distanceTable = pd.DataFrame([[-1 for x in range(length)] for x in range (length)])
    index = list(df['name'])
    distanceTable.columns = index
    distanceTable.index = index
    for tran1 in distanceTable.index:
        for tran2 in distanceTable.columns:
            boolean1 = list(df[df['name'] == tran1].boolean)[0]
            boolean2 = list(df[df['name'] == tran2].boolean)[0]
            length1 = float(sum(boolean1))
            length2 = float(sum(boolean2))
            overlapLength = sum([a and b for a, b in zip(boolean1, boolean2)])
            # How the number in the distance table is calculated:
            distance = (length1 + length2 - 2 * overlapLength) / (length1 + length2-overlapLength)
            distanceTable.loc[tran1,tran2] = distance

    # Group transcripts, n_clusters set how mant groups should be assigned
    for i in range(cluster_num):
        group = KMeans(n_clusters=i+1).fit_predict(distanceTable)
        global groupName
        groupName = 'group%s' %str(i+1)
        colorName = 'color%s' %str(i+1)
        df[groupName] = group
        df[colorName] = df.apply(assignColor, axis=1)
    return df

def assignColor(row):
    if row[groupName] == 0:
        return 'blue'
    elif row[groupName] == 1:
        return 'green'
    elif row[groupName] == 2:
        return 'red'
    elif row[groupName] == 3:
        return 'orange'
    elif row[groupName] == 4:
        return 'yellow'

def FP(tranname):
    if '-' not in tranname:
        FullPar = tranname.split('/')[1][1:]
        FPList = FullPar.split('p')
        fullScore = int(FPList[0])
        partialScore = int(FPList[1])
        return fullScore, partialScore
    else:
        return 0, 0

class Transcript (object):
    '''Just a struct actually, containing data about a transcript.'''

    def __init__ (self, name, score=None, annot=False):

        self.name    = name
        self.score   = score
        self.annot   = annot            # transcript comes from annotations?
        self.tranIx  = None             # y-axis coordinate of transcript
        self.exons   = list()           # Exon objects for this transcript
        self.blocks  = set()            # blocks where this transcript has exon(s)
        self.regions = set()            # regions where this transcript has exon(s)

        # What's the difference between a block and a region? Every
        # exon boundary defines a new region. A new block occurs only
        # when exon coverage transitions from 0 to >0. The example
        # below comprises 5 regions, but only one block.

        #    ==============
        #           ==============
        #                      ==================
        #    |      |     |    | |              |

        # I originally used block occupancy to group similar
        # transcripts in orderTranscripts. But that turned out not to
        # work very well, so I invented regions as a finer-grain
        # version of that idea.

class Exon (object):
    '''Struct containing data about an exon.'''

    def __init__ (self, tran, name, start, end, strand, QScore=None):

        self.tran     = tran            # Transcript object containing this exon
        self.name     = name
        self.start    = start
        self.end      = end
        self.strand   = strand
        self.QScore   = QScore
        self.block    = None            # block number where this exon resides
        self.adjStart = None            # start of exon in phony x-axis coordinates
        self.leading  = 0               # number of leading softclipped bases
        self.trailing = 0               # number of trailing softclipped bases

class Block (object):
    '''Struct for plot block.'''

    # A plot block is a vertical span representing a range contiguous
    # bases. Plot blocks are separated by vertical lines representing
    # regions of the reference, of unspecified length, which contain
    # no exons.

    # one implication of that scheme is that the x axis of the plot is
    # meaningless: it represents neither genomic nor RNA sequence range.

    def __init__ (self, start, end, boundary):

        self.start    = start          # actual genomic start coord
        self.end      = end            # actual genomic end coord
        self.boundary = boundary       # right-hand boundary x-coord in phony space
        self.annot    = False          # block contains annotation exons?
