import argparse
import pandas as pd

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadatacsv",
                        default=None,
                        type=str,
                        help="full path to metadata file produced from generate_fia_samples.py")
    parser.add_argument("--data_dir",
                        default=None,
                        type=str,
                        help="full path to directory containing all TREE csvs")
    parser.add_argument("--multipolyfile",
                        default=False,
                        type=str,
                        help="full path to a json containing a polygon or multipolygon geometry with which to subset points")
    parser.add_argument("--type",
                        default=False,
                        type=str,
                        help="one of 'hist' 'scatter' or 'bar' ")
    parser.add_argument("--cnlist",
                        default=False,
                        nargs='*',
                        help="Used to pull only a specific list of CNs")
    return(parser)

import fia_stats.poly_plt_hist as pph
import fia_stats.only_trees_bars as otb
from fia_sample_generation.shapetools import polys_from_jsons as pfj
from fia_sample_generation.shapetools import pull_plot_cns as ppc
import fia_sample_generation.sub_reals as sr

def main():
    args = create_parser().parse_args()
    fiametadata=pd.read_csv(args.metadatacsv)
    if args.multipolyfile:
        polylist= pfj(args.multipolyfile)
        print(polylist)
        args.cnlist= ppc(fiametadata,polylist,flatten=False)
    # while True:
    #     if args.cnlist is False:
    #         print('Please specify either a multipolyfile or a cnlist (single list or list of lists to compare)')
    #         break
    if len(args.cnlist[0]) > 1:
        dflist=[]
        treedflist=[]
        for list in args.cnlist:
            print(list)
            subdf=fiametadata[fiametadata['PLT_CN'].isin(list)]
            dflist.append(subdf)
    else:
        dflist=[fiametadata]
    pph.plot_hists(column='Total_Live_Bio', dflist=dflist, nonzeros=True, savehists=False, histtitles=['TAHOEMAMMOTH', 'SOUTHERNMOST', 'EASTNV'], comparepreds=False, log=False, randsample=False)
    # plot_scatters('qmd', 'counts', mergeddflist, nonzeros=True, colors='total',
    #               scattitles=['TAHOEMAMMOTH', 'SOUTHERNMOST', 'EASTNV'], bin=[1, 5, 10, 20, 30], sumbins=True)
    # plot_bars('qmd', 'counts', mergeddflist, bin=[1, 5, 10, 20, 30], nonzeros=False,
    #           bartitles=['Tahoe Area FIA Plots', 'SOUTHERNMOST', 'EASTNV'], xlabels=[2, 7, 15, 25, 40])


if __name__ == "__main__":
    main()
