'''
Given a directory full of downloaded FIA TREE, PLOT, and COND csvs,
generates samples for the Chimera CNN architecture for use with the EEcollection tool

'''
import argparse
import pandas as pd
import numpy as np
import sys

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir",
                        default=None,
                        type=str,
                        help="full path to directory containing all COND,TREE,PLOT csvs")
    parser.add_argument("--zerosonly",
                        default=False,
                        type=str,
                        help="returns csv of only corrsponding zeros")
    parser.add_argument("--minyear",
                        default=2004,
                        type=int,
                        help="Trim the fia samples to those inventoried after this year")
    parser.add_argument("--multipolyfile",
                        default=False,
                        type=str,
                        help="full path to a json containing a polygon or multipolygon geometry with which to subset points")
    parser.add_argument("--cnlist",
                        default=False,
                        type=str,
                        help="Used to pull only a specific list of CNs")
    parser.add_argument("--savemerged",
                        default=[False,False,False],
                        nargs='*',
                        help="A list of outputnames in the order: 'merged_plots' 'merged_trees' 'merged_conds']")
    parser.add_argument("--cpus",
                        default=1,
                        type=int,
                        help="Number of processes to split grouped processing into.")
    parser.add_argument("--bins",
                        default=[1,5,10,20,30],
                        nargs='*',
                        help="list of bin values including the less than and greater than of lowest and highest values, respectively")
    parser.add_argument("--mixthresh",
                        default=.8,
                        type=int,
                        help="mixed tree types ratio to classify as a mixed stand")
    parser.add_argument("--deadthresh",
                        default=.8,
                        type=int,
                        help="dead tree types ratio to classify as a dead stand")
    parser.add_argument("--genids",
                        default=0,
                        type=int,
                        help="set to 1 to create uniqid column")
    parser.add_argument("--heights",
                        default=False,
                        type=bool,
                        help="set to True to calculate average height")
    return(parser)

import fia_sample_generation.sub_reals as sr
import fia_sample_generation.generate_metadata as gm
def main():
    args = create_parser().parse_args()
    outnames=args.savemerged
    print(outnames)
    plotlist= sr.generate_plotlist(args.data_dir)
    masterplots = sr.merge_plots(plotlist, args.data_dir,args.minyear)
    statecdlist = sr.generate_statecodes(masterplots)
    realcoorddf= sr.merge_actuals(statecdlist)
    if args.cnlist:
        args.cnlist=pd.read_csv(args.cnlist)['CN'].tolist()
    if args.multipolyfile:
        from fia_sample_generation.shapetools import polys_from_jsons as pfj
        from fia_sample_generation.shapetools import pull_plot_cns as ppc
        from fia_sample_generation.shapetools import polys_from_anything as pfa
        polylist= pfa(args.multipolyfile)
        args.cnlist= ppc(realcoorddf,polylist,flatten=True)
        realcoorddf=realcoorddf[realcoorddf['CN'].isin(args.cnlist)]
    if args.genids == 1:
        args.genids=True
    finalplotdf = sr.sub_coords(realcoorddf, masterplots, outnames[0], genids=args.genids)
    print(finalplotdf.columns)
    if args.zerosonly:
        finalplotdf=finalplotdf[['CN','UNIQID','INVYR', 'STATECD', 'PLOT_STATUS_CD', 'LAT','LON']]
        print(finalplotdf.columns)
        extracolumns=['<1','1','5','10','20','>30','Class', 'QMD', 'Total_Live_Bio','CANOPY_CVR', 'BASAL_AREA', 'LAND_COVER_CLASS_CD', 'CONDPROP_UNADJ','DWM_FUELBED_TYPCD']
        for col in extracolumns:
            finalplotdf[col]=pd.Series(np.zeros((len(finalplotdf))))
        finalplotdf=finalplotdf[(finalplotdf['PLOT_STATUS_CD'] == 2) |(finalplotdf['PLOT_STATUS_CD'] == 3)]
        finalplotdf.to_csv('zeroscsv.csv')
        print('Zeros DF generated!')
        sys.exit()
    finaltreedf=sr.merge_trees_csvs(args.data_dir,args.minyear,args.cnlist,outname=outnames[1])
    finalconddf=sr.merge_conds_csvs(args.data_dir,args.minyear,args.cnlist,outname=outnames[2])
    all_metadata, binnames=gm.form_all_metadata(trees=finaltreedf,cpus=args.cpus,conds=finalconddf,bins=args.bins,deadthresh=args.deadthresh,mixthresh=args.deadthresh)
    final_metadata=pd.merge(all_metadata,finalplotdf,how='left',left_on='PLT_CN',right_on='CN')
    final_metadata = final_metadata.rename(index=str,columns={'INVYR_x':'INVYR'})
    columns=np.append(binnames,['Class', 'QMD', 'Total_Live_Bio','AVG_Height','CANOPY_CVR', 'BASAL_AREA', 'PLT_CN', 'INVYR', 'STATECD', 'PLOT_STATUS_CD', 'LAT','LON', 'LAND_COVER_CLASS_CD', 'CONDPROP_UNADJ','DWM_FUELBED_TYPCD'])
    final_metadata=final_metadata[columns].dropna()
    final_metadata.to_csv('11states_2019update.csv', columns=columns, index= False)
    print('Processing complete! You processed ' + str(len(finaltreedf)) + ' trees across ' + str(len(finalplotdf)) + ' locations to produce ' + str(len(final_metadata)) + ' complete samples.')


if __name__ == "__main__":
    main()
