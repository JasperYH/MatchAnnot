import os, sys
sys.path.append(os.getcwd())

import getGene
from getGene import *
from bokeh.plotting import Figure
from bokeh.models import ColumnDataSource, HoverTool, HBox, VBoxForm, PanTool, WheelZoomTool, BoxZoomTool, ResetTool, ResizeTool, PreviewSaveTool
from bokeh.io import curdoc
from bokeh.models.widgets import Slider, Select, TextInput, Button
from sklearn.cluster import KMeans
import pandas as pd

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
    opt = getParams(GTF.value.strip(), [Matches.value.strip()], Gene.value.strip(), annotations=Annotations)
    opt.gene = Gene.value.strip()
    global tranList
    tranList, exonList = selectGene(opt)
    forwardStrand = '-' if opt.flip else '+'
    if exonList[0].strand == forwardStrand:
        exonList.sort(key=lambda x: x.start)               # sort the list by start position
        blocks = assignBlocks (opt, exonList)              # assign each exon to a block
    else:
        exonList.sort(key=lambda x: x.end, reverse=True)   # sort the list by decreasing end position
        blocks = assignBlocksReverse (opt, exonList)       # assign each exon to a block -- backwards

    findRegions (tranList)                       # determine regions occupied by each transcript
    tranNames = orderTranscripts (tranList)

    global length, df, colorDF,boundaryDF, height
    length = len(tranNames)
    height = 900
    if round(900/float(2*length)) < 6:
        height = 12*length
    p.plot_height = height
    print p.plot_height
    p.title = "Transcript of %s" % opt.gene
    p.y_range.factors = tranNames[::-1]

    if Group.value == "on":
        colorDF = groupTran(tranList, exonList, 5)
    else:
        colorDF = None

    getExonData(exonList, colorDF)

    source.data = dict(
        xs=df['xs'],
        ys=df['ys'],
        color=df['colors'],
        line_alpha=df['alpha'],
        x=df['circlex'],
        y=df['circley'],
        QScore=df['QScore'],
        start=df['start'],
        end=df['end'],
        line_width=df['line_width'],
    )

    getBoundaryData(blocks)
    boundarySource.data = dict(
        xs=boundaryDF['x'],
        ys=boundaryDF['y'],
    )

def updateFP(attrname, old, new):
    global df
    df['alpha'] = df.apply(greaterFP, axis=1)
    source.data = dict(
        xs=df['xs'],
        ys=df['ys'],
        color=df['colors'],
        line_alpha=df['alpha'],
        x=df['circlex'],
        y=df['circley'],
        QScore=df['QScore'],
        start=df['start'],
        end=df['end'],
        line_width=df['line_width'],
    )

def updateGroup(attrname, old, new):
    colors = list()
    global colorDF, df
    for index, row in df.iterrows():
        color = getColor(row['tran'], colorDF)
        colors.append(color)
    df['colors'] = colors
    source.data = dict(
        xs=df['xs'],
        ys=df['ys'],
        color=df['colors'],
        line_alpha=df['alpha'],
        x=df['circlex'],
        y=df['circley'],
        QScore=df['QScore'],
        start=df['start'],
        end=df['end'],
        line_width=df['line_width'],
    )

def getExonData(exonList, colorDF):
    global df
    columns = ['name', 'xs', 'ys', 'colors', 'circlex', 'circley', 'QScore',
                'start', 'end', 'tran', 'full', 'partial']
    df = pd.DataFrame(columns=columns)
    for myExon in exonList:
        exonSize = myExon.end - myExon.start + 1
        adjStart = myExon.adjStart
        if colorDF is not None:
            color = getColor(myExon.tran.name, colorDF)
        else:
            if myExon.tran.annot:
                color = 'purple'
            else:
                color = 'blue'
        xs = [adjStart, adjStart+exonSize]
        ys = [length-(myExon.tran.tranIx), length-(myExon.tran.tranIx)]
        circlex = (adjStart + adjStart+exonSize)/2
        circley = length-(myExon.tran.tranIx)
        data = pd.Series([myExon.name, xs, ys, color, circlex, circley,
                        myExon.QScore, myExon.start, myExon.end,
                        myExon.tran.name, myExon.full, myExon.partial], index=[columns])

        df = df.append(data,ignore_index=True)
    df['line_width'] = round(height/float(2*length))
    df['alpha'] = 1

def getColor(exonName, colorDF):
    if exonName not in list(colorDF.name):
        color = 'purple'
    else:
        row = colorDF.loc[colorDF['name'] == exonName]
        colorName = 'color%s' %str(Cluster.value)
        color = row[colorName]
    return color

def getBoundaryData(blocks):
    boundaryX = list()
    for bound in blocks:
        boundaryX.append([bound.boundary, bound.boundary])
    boundaryY = [[0, length+1] for x in range(len(boundaryX))]
    global boundaryDF
    newBoundaryDF = pd.DataFrame()
    newBoundaryDF['x'] = boundaryX
    newBoundaryDF['y'] = boundaryY
    boundaryDF = newBoundaryDF
    return

def greaterFP(row):
    alphaVal = Alpha.value
    if row['full'] < Full.value or row['partial'] < Partial.value:
        return alphaVal
    else:
        return 1.0

def saveDF(attrname, old, new):
    outDF = pd.DataFrame()
    fileName = Save.value.strip()
    name = list()
    full = list()
    partial = list()
    start = list()
    end = list()
    score = list()
    for tran in tranList:
        if tran.annot == False:
            name.append(tran.name)
            full.append(tran.full)
            partial.append(tran.partial)
            start.append(tran.start)
            end.append(tran.end)
            score.append(tran.score)
    outDF['name'] = name
    outDF['full'] = full
    outDF['partial'] = partial
    outDF['start'] = start
    outDF['end'] = end
    outDF['score'] = score
    outDF.to_csv(fileName)

class getParams(object):
    def __init__(self, gtf, matches, gene, formatt="standard", omit=None,
                 show=None, howmany=None, nodups=None, minlen=None, maxlen=None, output="exon.png",
                 flip=None, yscale=1.0, details=None, fasta=None, title=None, notes=None, full=None,
                 partial=None, highsupport=None, annotations=None):
        self.gtf = gtf
        self.matches = matches
        self.gene = gene
        self.format = formatt
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
        self.annotations = annotations

Annotations = getAnnotations(getParams("gencode.vM8.annotation.gtf", None, None))
GTF = TextInput(title="Enter the name of GTF file", value="gencode.vM8.annotation.gtf")
Matches = TextInput(title="Enter the name of pickle file from MatchAnnot", value="matches.pickle")
Gene = TextInput(title="Select gene to visualize")
Alpha = Slider(title="Alpha value of exons", value=1.0, start=0, end=1.0, step=0.1)
Full = Slider(title="Full support threshold", value=0, start=0, end=30, step=1.0)
Partial = Slider(title="Partial support threshold", value=0, start=0, end=50, step=1.0)
Cluster = Slider(title="The number of clusters", value=3, start=1, end=5, step=1.0)
Group = Select(title="Group isoform or not", value="on", options=["on", "off"])
Save = TextInput(title="Enter the file name to save data (e.g. gene.csv)")

boundarySource = ColumnDataSource(data=dict(xs=[], ys=[]))
source = ColumnDataSource(data=dict(xs=[], ys=[], color=[], line_alpha=[],
                    x=[], y=[], QScore=[], start=[], end=[], line_width=[]))

df = pd.DataFrame()
boundaryDF = pd.DataFrame()
colorDF = pd.DataFrame()
outDF = pd.DataFrame()

hover = HoverTool(tooltips=[
    ("QScore", "@QScore"),
    ("start", "@start"),
    ("end", "@end")
])
tools = [PanTool(), WheelZoomTool(), BoxZoomTool(), ResetTool(), ResizeTool(), PreviewSaveTool(), hover]
p = Figure(plot_height=900, plot_width=1200, title="", y_range=['tran'], tools=tools)
p.xgrid.grid_line_color = None
p.ygrid.grid_line_color = None
p.circle(x="x", y="y", source=source, size=5, color="color", line_color=None)
p.multi_line(xs="xs", ys="ys", source=boundarySource, color="black",
            line_width=2, line_alpha=0.4, line_dash="dotted")
p.multi_line(xs="xs", ys="ys", source=source, color="color",
            line_width="line_width", line_alpha='line_alpha')

controls = [GTF, Matches, Gene, Alpha, Full, Partial, Cluster, Group, Save]
Gene.on_change('value', updateGene)
Full.on_change('value', updateFP)
Partial.on_change('value', updateFP)
Cluster.on_change('value', updateGroup)
Save.on_change('value', saveDF)

inputs = HBox(VBoxForm(*controls), width=300)
curdoc().add_root(HBox(inputs, p, width=1300))
