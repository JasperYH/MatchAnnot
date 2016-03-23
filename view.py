import getGene
from getGene import *
from bokeh.plotting import Figure, output_file, show
from bokeh.models import *
from bokeh.io import curdoc
from bokeh.models.widgets import Slider, Select, TextInput
from sklearn.cluster import KMeans
import pandas as pd

def main ():

    logger.debug('version %s starting' % VERSION)

    opt, args = getParms()

    # Find all the exons in all the transcripts for the gene, put them
    # in a list.

    tranList = list()                                      # list of Transcript objects
    exonList = list()                                      # list of Exon objects

    if opt.gtf is not None:
        getGeneFromAnnotation (opt, tranList, exonList)    # lists will be changed
    if opt.matches is not None:
        getGeneFromMatches (opt, tranList, exonList)       # lists will be changed
    if len(exonList) == 0:
        raise RuntimeError ('no exons found for gene %s in annotation or match files' % opt.gene)

    forwardStrand = '-' if opt.flip else '+'
    if exonList[0].strand == forwardStrand:
        exonList.sort(key=lambda x: x.start)               # sort the list by start position
        blocks = assignBlocks (opt, exonList)              # assign each exon to a block
    else:
        exonList.sort(key=lambda x: x.end, reverse=True)   # sort the list by decreasing end position
        blocks = assignBlocksReverse (opt, exonList)       # assign each exon to a block -- backwards

    findRegions (tranList)                       # determine regions occupied by each transcript

    tranNames = orderTranscripts (tranList)

    output_file("transcript.html")
    p = Figure(plot_width=1000, plot_height=750)
    df = groupTran(tranList, exonList, opt.group)
    length = len(tranNames)
    for myExon in exonList:
        exonSize = myExon.end - myExon.start + 1
        adjStart = myExon.adjStart
        for index, row in df.iterrows():
            name = row['name']
            groupColor = 'purple'
            if name in myExon.name:
                groupColor = row['color']
                break
        p.line([adjStart, adjStart+exonSize], [length-(myExon.tran.tranIx+1), length-(myExon.tran.tranIx+1)], line_width=20, line_color=groupColor)

    f_range = FactorRange(factors=tranNames[::-1])
    p.extra_y_ranges = {"Tran": f_range}
    new_axis = CategoricalAxis(y_range_name="Tran")
    p.add_layout(new_axis, 'left')
    show(p)


def getParms ():                       # use default input sys.argv[1:]

    # Note that --matches can be specified more than once, and opt.matches is a list.

    parser = optparse.OptionParser(usage='%prog [options]', version=VERSION)

    parser.add_option ('--gtf',     help='annotations file, in format specified by --format (optional)')
    parser.add_option ('--format',  help='format of annotation file: standard, alt, pickle (def: %default)', \
                           type='choice', choices=['standard', 'alt', 'pickle'])
    parser.add_option ('--matches', help='pickle file from matchAnnot.py (optional)', action='append')
    parser.add_option ('--gene',    help='gene to plot (required)')
    parser.add_option ('--omit',    help='clusters to ignore, e.g., c1234,c2345,c3456')
    parser.add_option ('--show',    help='clusters to force shown, even if underpopulated')
    parser.add_option ('--howmany', help='how many clusters to plot (def: all)', type='int')
    parser.add_option ('--nodups',  help='discard exact duplicate clusters (def: keep all)', action='store_true')
    parser.add_option ('--minlen',  help='minimum length of plotted cluster (def: all)', type='int')
    parser.add_option ('--maxlen',  help='maximum length of plotted cluster (def: all)', type='int')
    parser.add_option ('--output',  help='output plot file name (def: %default)')
    parser.add_option ('--flip',    help='reverse plot orientation (def: mRNA 5\' on left)', action='store_true')
    parser.add_option ('--yscale',  help='amount by which to scale Y axis (def: %default)', type='float')
    parser.add_option ('--details', help='output file name for details in text format (def: no output)')
    parser.add_option ('--fasta',   help='output directory for fasta files for chosen clusters (def: no output)')
    parser.add_option ('--title',   help='title for top of figure')
    parser.add_option ('--notes',   help='file of notations to add to plot (experts only!)')
    parser.add_option ('--full',    help='full support threshold', type='int')
    parser.add_option ('--partial', help='partial support threshold', type='int')
    parser.add_option ('--highsupport',   help='only clusters have higher full and partial support than threshold will be plotted', action='store_true')
    parser.add_option ('--group', help='classify transcripts into how many groups', type='int')
    parser.set_defaults (format='standard',
                         output=DEF_OUTPUT,
                         yscale=DEF_YSCALE,
                         group=2,
                         )

    opt, args = parser.parse_args()

    return opt, args

if __name__ == "__main__":
    main()
