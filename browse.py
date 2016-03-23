path = "/home/jasper/Projects/MatchAnnot-master"
import sys
sys.path.append(path)
import getGene
from getGene import *
from bokeh.plotting import Figure, output_file
from bokeh.models import *
from bokeh.io import curdoc
from bokeh.models.widgets import Slider, Select, TextInput
from sklearn.cluster import KMeans
import pandas as pd

df = pd.DataFrame()
tickDF = pd.DataFrame()

tranNames = [None]

def selectGene(opt):
    tranList = list()                                      # list of Transcript objects
    exonList = list()                                      # list of Exon objects
    if opt.gtf is not None:
        getGeneFromAnnotation (opt, tranList, exonList)    # lists will be changed
    if opt.matches is not None:
        getGeneFromMatches (opt, tranList, exonList)       # lists will be changed
    if len(exonList) == 0:
        raise RuntimeError ('no exons found for gene %s in annotation or match files' % opt.gene)
    return tranList, exonList

def updateGene(attrname, old, new):
    opt.gene = Gene.value.strip()
    tranList, exonList = selectGene(opt)
    forwardStrand = '-' if opt.flip else '+'
    if exonList[0].strand == forwardStrand:
        exonList.sort(key=lambda x: x.start)               # sort the list by start position
        blocks = assignBlocks (opt, exonList)              # assign each exon to a block
    else:
        exonList.sort(key=lambda x: x.end, reverse=True)   # sort the list by decreasing end position
        blocks = assignBlocksReverse (opt, exonList)       # assign each exon to a block -- backwards

    findRegions (tranList)                       # determine regions occupied by each transcript
    global tranNames
    tranNames = orderTranscripts (tranList)

    xs = list()
    ys = list()
    colors = list()
    full = list()
    partial = list()
    maxX = 0

    groupdf = groupTran(tranList, exonList, Cluster.value)
    for myExon in exonList:
        exonSize = myExon.end - myExon.start + 1
        adjStart = myExon.adjStart
        exonName = myExon.name
        fullScore, partialScore = FP(myExon.name)
        full.append(fullScore)
        partial.append(partialScore)
        for index, row in groupdf.iterrows():
            name = row['name']
            groupColor = 'purple'
            if name in myExon.name:
                groupColor = row['color']
                break
        if maxX < adjStart+exonSize:
            maxX = adjStart+exonSize
        xs.append([adjStart, adjStart+exonSize])
        ys.append([myExon.tran.tranIx+1, myExon.tran.tranIx+1])
        colors.append(groupColor)

    exonNum = len(xs)
    df['xs'] = xs
    df['ys'] = ys
    df['alpha'] = [1 for x in range(exonNum)]
    df['colors'] = colors
    df['full'] = full
    df['partial'] = partial

    numTran = len(tranNames)

    tickSource.data=dict(
        x = [maxX+50 for x in range(numTran)],
        y = [x+1 for x in range(numTran)],
        text = tranNames,
    )

    source.data = dict(
        xs=xs,
        ys=ys,
        color=colors,
        line_alpha=list(df['alpha']),
    )

def updateFP(attrname, old, new):
    df['alpha'] = df.apply(greaterFP, axis=1)
    source.data = dict(
        xs=list(df['xs']),
        ys=list(df['ys']),
        color=list(df['colors']),
        line_alpha=list(df['alpha']),
    )

def greaterFP(row):
    alphaVal = Alpha.value
    if row['full'] < Full.value or row['partial'] < Partial.value:
        return alphaVal
    else:
        return 1.0

class getParams(object):
    format='standard'
    output = 'exons.png'
    yscale = 1.0
    def __init__(self, gtf, matches, gene, format="standard", omit=None,
                 show=None, howmany=None, nodups=None, minlen=None, maxlen=None, output="exon.png",
                 flip=None, yscale=1.0, details=None, fasta=None, title=None, notes=None, full= None,
                 partial=None, highsupport=None):
        self.gtf = gtf
        self.matches = matches
        self.gene = gene
        self.omit = omit
        self.show = show
        self.howmany = howmany
        self.nodups = nodups
        self.minlen = minlen
        self.maxlen = maxlen
        self.output = output
        self.flip = flip
        self.yscale = yscale
        self.details = details
        self.fasta = fasta
        self.title = title
        self.notes = notes
        self.highsupport = highsupport
        self.full = full
        self.partial = partial

opt = getParams('gencode.vM8.annotation.gtf', ['matches.pickle'], 'Trak2')

Gene = TextInput(title="Select gene to visualize")
Alpha = Slider(title="Alpha value of exons", value=1.0, start=0, end=1.0, step=0.1)
Full = Slider(title="Full support threshold", value=0, start=0, end=10, step=1.0)
Partial = Slider(title="Partial support threshold", value=0, start=0, end=10, step=1.0)
Cluster = Slider(title="The number of clusters", value=3, start=1, end=10, step=1.0)

source = ColumnDataSource(data=dict(xs=[], ys=[], color=[], line_alpha=[]))
tickSource = ColumnDataSource(data=dict(x=[], y=[], text=[]))

p = Figure(plot_height=600, plot_width=800, title="")
p.multi_line(xs="xs", ys="ys", source=source, color="color", line_width=20, line_alpha='line_alpha')
p.text(x="x", y="y", text="text", source=tickSource)

controls = [Gene, Alpha, Full, Partial, Cluster]
Gene.on_change('value', updateGene)
Full.on_change('value', updateFP)
Partial.on_change('value', updateFP)

inputs = HBox(VBoxForm(*controls), width=300)
curdoc().add_root(HBox(inputs, p, width=1100))
