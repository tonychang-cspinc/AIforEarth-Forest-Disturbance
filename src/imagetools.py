import os 
import osgeo
import ogr
import gdal
import osr
import numpy as np
import geopandas as gpd
import subprocess
import skimage.transform as st
from skimage.measure import compare_ssim as ssim
from skimage import exposure
import src.histogram_matching as hm
import src.eetools as eet
import matplotlib.pyplot as plt
from matplotlib import ticker

def rmse(x, y):
    return np.sqrt(np.linalg.norm(x - y))

def showim(image,size=(6,3),title=False):
    plt.figure(figsize=size)
    plt.imshow(image)
    if title:
        plt.title(title)
    plt.show()
    plt.close()

def read_and_transform_image(impath, scale=1):
    return np.swapaxes(gdal.Open(f'{impath}').ReadAsArray(),0,-1)*scale

def match_naip_and_sentinel(naipfile, sentinelfile, export_file=None, channels_first=False):
    SENTINEL_SCALE_CONST = 10000
    sentinel = np.swapaxes(read_and_transform_image(f'{sentinelfile}'),0,1)
    naip = np.swapaxes(read_and_transform_image(f'{naipfile}',scale=(1/255.)),0,1)
    #sentinel = np.nan_to_num(sentinel)/np.max(sentinel)
    sentinel = np.nan_to_num(sentinel)/SENTINEL_SCALE_CONST
    #scale sentinel between 0 and 1
    naip = np.nan_to_num(naip)
    sentinel_hist_matched = hm.match_histograms(image=sentinel,reference=naip, multichannel=True)
    if export_file is not None:
        export_tif(sentinel_hist_matched, sentinelfile, export_file)
    if channels_first:
        naip = np.moveaxis(naip,-1,0)
        sentinel = np.moveaxis(sentinel,-1,0)
        sentinel_hist_matched = np.moveaxis(sentinel_hist_matched,-1,0)
    return naip, sentinel, sentinel_hist_matched

#def match_histograms(src, tgt):
#    matchedsent = hm.match_histograms(sent, naip, multichannel=True)
    
def match_and_plot_sent(localpath=None):
    if localpath is None:
        localpath = './data/eeout'
    stretchedsent=[]
    pairs = [('sent_present.tif','naipmatch_present.tif'),\
            ('sent_postsent.tif','naipmatch_postsent.tif')]
    sent_pairs = pairs[0][0], pairs[1][0]
    naip_pairs = pairs[0][1], pairs[1][1]
    if all([os.path.exists(f'{localpath}/{n}') for n in naip_pairs]):
        for impair in pairs:
            sent = np.swapaxes(read_and_transform_image(f'{localpath}/{impair[0]}'),0,1)
            naip = np.swapaxes(read_and_transform_image(f'{localpath}/{impair[1]}',scale=(1/255.)),0,1)
            sent = np.nan_to_num(sent)/np.max(sent.flatten())
            naip = np.nan_to_num(naip)
            matchedsent = hm.match_histograms(sent, naip, multichannel=True)
            print('plotting showing',impair)
            print(sent.shape,naip.shape)
            ssim_val = ssim(matchedsent,naip,multichannel=True)
            rmselist = [ssim_val]
            for channel in range(3):
                rmselist.append(rmse(matchedsent[:,:,channel],naip[:,:,channel]))
            titlezip = list(zip(rmselist,['SSIM','R_RMSE','G_RMSE','B_RMSE']))
            counter = 0
            for im in [sent, naip, matchedsent]:
                if counter < 2:
                    showim(im)
                else:
                    showim(im,title=titlezip)
                counter+=1
            stretchedsent.append(matchedsent)
    elif os.path.exists(f'{localpath}/{naip_pairs[0]}'):
        for impair in sent_pairs:
            sent = np.swapaxes(read_and_transform_image(f'{localpath}/{impair}'),0,1)
            naip = np.swapaxes(read_and_transform_image(f'{localpath}/{naip_pairs[0]}',(1/255)),0,1)
            sent=np.nan_to_num(sent)/np.max(sent.flatten())
            naip=np.nan_to_num(naip)
            matchedsent=hm.match_histograms(sent, naip, multichannel=True)
            print('plotting showing',impair)
            print(sent.shape,naip.shape)
            ssim_val=ssim(matchedsent,naip,multichannel=True)
            rmselist=[ssim_val]
            for channel in range(3):
                rmselist.append(rmse(matchedsent[:,:,channel],naip[:,:,channel]))
            titlezip=list(zip(rmselist,['SSIM','R_RMSE','G_RMSE','B_RMSE']))
            counter=0
            for im in [sent, naip, matchedsent]:
                if counter < 2:
                    showim(im)
                else:
                    showim(im,title=titlezip)
                counter+=1
            stretchedsent.append(matchedsent)
    else:
        for impair in sent_pairs:
            sent = np.swapaxes(read_and_transform_image(f'{localpath}/{impair}'),0,1)
            sent = np.nan_to_num(sent)/np.max(sent.flatten())
            stretchedsent.append(sent)
    return(stretchedsent)

def write_stretchedsent(matchedsent, localpath=None):
    counter=0
    if localpath is None:
        localpath = './data/eeout'
    for sentarr in matchedsent:
        if counter==0:
            ref ='sent_present.tif'
        else:
            ref ='sent_postsent.tif'
        outname = f'{localpath}/{ref[:-4]}naipstretched.tif'
        reffile = gdal.Open(f'{localpath}/{ref}')
        xsize = sentarr.shape[1]
        ysize = sentarr.shape[0]
        driver = gdal.GetDriverByName('GTiff')
        srs = osr.SpatialReference()
        wkt = reffile.GetProjection()
        srs.ImportFromWkt(wkt)
        out = driver.Create(outname, xsize, ysize, 3, gdal.GDT_Float64)
        out.SetGeoTransform(reffile.GetGeoTransform())
        out.SetProjection(srs.ExportToWkt())
        print("Saving ",outname)
        for band in range(1, 4):
            exportband = sentarr[:, :, band - 1].astype(np.float64)
            out.GetRasterBand(band).WriteArray(exportband)
        del out
        counter+=1

def quick_qaqc(predlist,reshape,filelist):
    predims = []
    for tileind in range(len(predlist)):
        print('Showing outputs from',filelist[tileind])
        predvalues = np.reshape(predlist[tileind][1],(reshape[0],reshape[1],4))
        predims.append(predvalues)
        titlelist = ['Biomass','QMD','Basal Area', 'Canopy Cover']
        for channel in range(4):
            plt.imshow(predvalues[:,:,channel])
            plt.title(titlelist[channel])
            plt.colorbar(shrink=.5)
            plt.show()
            plt.close()
    return(predims)

def clip_preds(arr):
    clipvals = [(0,800),(0,70),(0,600),(0,100),(0,50),(0,30),(0,400),(0,1.5)]
    for band in range(arr.shape[-1]):
        arr[:,:,band]=np.nan_to_num(np.clip(arr[:,:,band],clipvals[band][0],clipvals[band][1]))
    return(arr)

def classify_lcpband(array,lowval,highval,sublow,subhigh):
    outshape=array.shape
    array=array.flatten()
    zerobin=np.where((array > -9999) & (array<=lowval))[0].copy()
    print(len(zerobin))
    lowbin=np.where((array>lowval) & (array <= highval))[0].copy()
    print(len(lowbin))
    highbin=np.where(array>highval)[0].copy()
    print(len(highbin))
    outarray=array.copy().flatten()
    for ind in [[zerobin,0],[lowbin,sublow],[highbin,subhigh]]:
        print(ind)
        outarray[ind[0]]=ind[1]
    print(np.max(outarray),np.min(outarray))
    return(outarray.reshape(outshape))

def write_preds(predims, reffile, prefix, output_dirc='output', names_list=None):
    reftif = gdal.Open(reffile)
    if names_list is None:
        names_list = ['biomass','QMD','basalarea','canopycover',\
                         'standheight','cbh','tpha','cbd']
    if not os.path.exists(output_dirc):
        os.mkdir(output_dirc)
    for im in range(len(predims)):
        outname = f'{output_dirc}/{prefix}_{names_list[i]}.tif'
        export_tif(image, ref_tif, outname, bands=1) 

def write_preds_deprecated(predims,output_dirc='output',split=True,localpath=None):
    if localpath is None:
        localpath = 'data/eeout'
    if not os.path.exists(output_dirc):
        os.mkdir(output_dirc)
    counter=0
    for tileind in range(len(predims)):
        predarr=clip_preds(predims[tileind])
        if counter == 0:
            reffile=gdal.Open(f'{localpath}/ancls_present.tif')
            outname=f'{output_dirc}/sent_pre_predictions.tif'
            if split:
                outvars=['biomass','QMD','basalarea','canopycover']
                if predarr.shape[-1] > 4:
                    outvars.extend(['standheight','cbh','tpha','cbd'])
                outnames=[outname[:-4] +'_sep_'+ x + '.tif' for x in outvars]
        else:
            reffile=gdal.Open(f'{localpath}/ancls_postsent.tif')
            outname=f'{output_dirc}/sent_post_predictions.tif'
            if split:
                outvars=['biomass','QMD','basalarea','canopycover']
                if predarr.shape[-1] > 4:
                    outvars.extend(['standheight','cbh','tpha','cbd'])
                outnames=[outname[:-4] +'_sep_'+ x + '.tif' for x in outvars]
        xsize = predarr.shape[0]
        ysize = predarr.shape[1]
        print(xsize,ysize)
        driver = gdal.GetDriverByName('GTiff')
        srs = osr.SpatialReference()
        wkt = reffile.GetProjection()
        srs.ImportFromWkt(wkt)
        if split:
            print(outnames)
            bandc=0
            hcbhbin=8.2
            lcbhbin=1.1
            lcbdbin=.015
            hcbdbin=.055
            cbhsublow=1.52
            cbhsubhigh=7.62
            cbdsublow=.035
            cbdsubhigh=.1
            for outname in outnames:
                if 'cbd' in outname:
                    predarr[:, :,bandc]=classify_lcpband(predarr[:, :,bandc],lcbdbin,hcbdbin,cbdsublow,cbdsubhigh)
                if 'cbh' in outname:
                    predarr[:, :,bandc]=classify_lcpband(predarr[:, :,bandc],lcbhbin,hcbhbin,cbhsublow,cbhsubhigh)
                out = driver.Create(outname, ysize, xsize, 1, gdal.GDT_Float32)
                out.SetGeoTransform(reffile.GetGeoTransform())
                out.SetProjection(srs.ExportToWkt())
                exportband = predarr[:, :,bandc].astype(np.float32)
                out.GetRasterBand(1).WriteArray(exportband)
                del out
                bandc+=1
        else:
            out = driver.Create(outname, ysize, xsize, predarr.shape[-1], gdal.GDT_Float32)
            out.SetGeoTransform(reffile.GetGeoTransform())
            out.SetProjection(srs.ExportToWkt())
            for band in range(1, predarr.shape[-1]+1):
                    exportband = predarr[:, :, band - 1].astype(np.float32)
                    out.GetRasterBand(band).WriteArray(exportband)
            del out
        counter+=1

def clip_raster_to_shp(raster,shp_filename,clipped_filename, \
                        shp_srs=None,\
                        dst_srs=None, \
                        dst_format='GTiff'):
    ds = gpd.read_file(shp_filename)
    if shp_srs is not None:
        ds.crs = shp_srs
    if type(raster) != osgeo.gdal.Dataset:
        raster = gdal.Open(raster)
    #raster_srs = raster.GetProjection()
    if dst_srs is not None:
        ds = ds.to_crs(dst_srs)
    ulx, lry, lrx, uly = ds.bounds.values[0]
    clip_bounds = [ulx, uly, lrx, lry]
    clipped = gdal.Translate(clipped_filename, raster, \
                            projWin=clip_bounds, format=dst_format)
    return clipped
 
def convert_to_ascii_deprecated(filenames,reflcp):
    for filename in filenames:
        lcp=gdal.Open(reflcp)
        ulx, xres, xskew, uly, yskew, yres  = lcp.GetGeoTransform()
        lrx = ulx + (lcp.RasterXSize * xres)
        lry = uly + (lcp.RasterYSize * yres)
        print(ulx,uly,lrx,lry)
        src_ds=gdal.Open(filename)
        #src_ds = gdal.Open(filename)
        # Open source dataset
        out_ds=gdal.Translate(filename[:-4] + '.asc', src_ds, format='AAIGrid', noData= -9999,outputSRS=lcp.GetProjection(),outputBounds = [ulx, uly, lrx, lry],width=lcp.RasterXSize, height=lcp.RasterYSize, outputType=gdal.GDT_Float32)

def convert_to_ascii(filenames,reflcp):
    for filename in filenames:
        lcp = gdal.Open(reflcp)
        nrows = lcp.RasterYSize
        ncols = lcp.RasterXSize
        gt = lcp.GetGeoTransform()
        xmin, ymin, xmax, ymax = gt[0], gt[3]+(gt[-1]*nrows), gt[0]+(gt[1]*ncols), gt[3]
        res = gt[1] 
        extent = [xmin, ymin, xmax, ymax]
        prj = lcp.GetProjection()
        srs = osr.SpatialReference(wkt=prj)
        #src_ds = gdal.Open(filename)
        # Open source dataset
        # Reproject the data first...
        dst_file = '/tmp/tmp.tif'
        print(f'warping {filename} and generating ascii..')
        src_ds = gdal.Warp(srcDSOrSrcDSTab=filename, \
            destNameOrDestDS=dst_file, \
            format='GTiff',\
            xRes=res, yRes=res,\
            dstSRS=srs,\
            dstNodata=-9999,\
            outputBoundsSRS=srs,\
            outputBounds=[xmin, ymin, xmax, ymax])
        asc_file = f'{os.path.splitext(filename)[0]}.asc'
        out_ds = gdal.Translate(asc_file, src_ds, format='AAIGrid', noData=-9999)

def mod_headers(asciis):
    for filename in asciis:
        with open(filename, 'r') as file:
            # read a list of lines into data
            data = file.readlines()
            del data[4:6]
            data.insert(4, 'cellsize     30\n')
            file.close()
        with open(filename, 'w') as file:
            file.writelines(data)

def export_tif(image, ref_tif, outname, bands=None, dtype=gdal.GDT_Float32, metadata=None, bandmeta=None, rasterdriver='GTiff', verbose=True):
    '''
    Input a numpy array image and a reference geotif 
    to convert image to geotiff of same geotransform
    and projection. Note, if alpha_mask is not None,
    creates a 4 channel geotiff (alpha as last channel)

    Parameters:
    -----------
    image - 3D <numpy array>
    ref_tif - Geotiff reference <str> filename or <gdal object> of same dimensions
    outname - <str> file name to output to (use .tif extension)
    dtype - <str> denoting data type for GeoTiff. Defaults to 8 bit image,
    but can use gdal.GDT_Float32
    '''
    if type(ref_tif) is not osgeo.gdal.Dataset:
        ref_tif = gdal.Open(ref_tif)
    gt = ref_tif.GetGeoTransform()
    proj = ref_tif.GetProjection()
    xsize = np.shape(image)[1] 
    ysize = np.shape(image)[0] 
    if bands is None:
        bands = ref_tif.RasterCount 
    driver = gdal.GetDriverByName(rasterdriver)
    out = driver.Create(outname, xsize, ysize, bands, dtype)
    out.SetGeoTransform(gt)
    out.SetProjection(proj)
    if metadata is not None:
        out.SetMetadata(metadata)
    if bands == 1:
        band = out.GetRasterBand(1)
        band.WriteArray(image)
        if bandmeta is not None:
            band.SetMetadata(bandmeta) 
            band.SetDescription(bandmeta)
    else:
        for i in range(bands):
            band = out.GetRasterBand(i+1)
            band.WriteArray(image[:,:,i]) #if we want a red image for a 4 channel
            if bandmeta is not None:
                band.SetMetadata(bandmeta[i]) 
                band.SetDescription(bandmeta[i])
    out = None
    if verbose:
        return(print('created %s'%(outname)))

def tif_to_ascii(tif_file, asc_file):
    out_ds = gdal.Translate(asc_file, tif_file, format='AAIGrid', noData=-9999)

def histogram_match_lcp(lcp_var, reflcp, srcfile, scale=None, outname=None, verbose=False, plotfolder=None):
    '''
    input a structure variable by name
    function will:
        1) open the lcp, get the band
        2) open the corresponding chimera asc
        3) scale and histogram match chimera
        4) write out match to tif 
        5) return matched array
    '''    
    var_map = {'Canopy cover':'canopycover', 'Canopy height':'standheight', \
              'Canopy base height':'cbh', 'Canopy bulk density':'cbd'}
    lcp_bands = [n.split('\n')[0][1:] \
                 for n in gdal.Info(reflcp).split('Description =')][1:]
    if scale is None:
        lcp_scale = [int(n.split('\n')[0].split('x')[1]) \
                 if len(n.split('\n')[0].split('x'))==2 else 1 \
                 for n in gdal.Info(reflcp).split('UNIT_NAME=')]
        scale_map = {key:value for (key,value) in zip(lcp_bands, lcp_scale)}
        scale = scale_map[lcp_var]
    i = [l for l in range(len(lcp_bands)) if lcp_var in lcp_bands[l]]
    lcp_arr = reflcp.ReadAsArray()[i][0]
    src = gdal.Open(srcfile).ReadAsArray()*scale
    matched =  exposure.match_histograms(src, lcp_arr)
    if verbose:
        print(f'Opening {srcfile}')
        print(f'Matched, {lcp_var} with scale {scale}!')
    #write it out with export
    if outname is not None:
        tif_file = f'{os.path.splitext(outname)[0]}.tif'
        export_tif(matched, reflcp, tif_file, bands=1, verbose=False)
        tif_to_ascii(tif_file, outname)
        print(f'Wrote {outname}')
    if plotfolder is not None:
        if not os.path.exists(plotfolder):
            os.makedirs(plotfolder)
        labels = ['Source', 'Ref', 'Matched']
        f, ax = plt.subplots(2,3,figsize=(20,13))
        for i, img in enumerate((src, lcp_arr, matched)):
            ax[0,i].imshow(img)
            ax[0,i].set_title(labels[i], fontsize=22)
            ax[0,i].tick_params(axis='both', which='major', labelsize=20)
            img_hist, bins = exposure.histogram(img)
            ax[1,i].plot(bins, img_hist/img_hist.max())
            img_cdf, bins = exposure.cumulative_distribution(img)
            ax[1,i].plot(bins, img_cdf)
            ax[1,i].tick_params(axis='both', which='major', labelsize=20)
        print(f'Writing {plotfolder}/{lcp_var}_match.png')
        outprefix = os.path.splitext(os.path.basename(srcfile))[0]
        f.savefig(f'{plotfolder}/{outprefix}_match.png')
    return(matched) 

def replace_lcp_layers(image, ref_lcp, outname, bands_to_replace):
    '''
    Input a numpy array image and a reference lcp
    at the respective bands same geotransform
    and projection. 

    Parameters:
    -----------
    image - 3D <numpy array> with first dimension as bands
    ref_lcp - LCP reference <gdal object> of same dimensions
    outname - <str> file name to output to (use .lcp extension)
    bands_to_replace - <list> of bands to replace in lcp
    '''
    gt = ref_lcp.GetGeoTransform()
    proj = ref_lcp.GetProjection()
    driver = gdal.GetDriverByName('GTiff')
    lcp_array = ref_lcp.ReadAsArray()
    xsize = np.shape(lcp_array)[2] 
    ysize = np.shape(lcp_array)[1] 
    bands = np.shape(lcp_array)[0]

    out = driver.Create(outname, xsize, ysize, bands, gdal.GDT_Int16)
    out.SetGeoTransform(gt)
    out.SetProjection(proj)
    for i in range(bands):
        band = out.GetRasterBand(i+1)
        band.WriteArray(lcp_array[i]) 
    for b in range(len(bands_to_replace)):
        band = out.GetRasterBand(bands_to_replace[b])
        band.WriteArray(image[b])
    out.FlushCache()
    out = None
    return(print('created %s'%(outname)))

'''            
def build_shapepoints(refasc,plotshape):
    outproj=gdal.Open(refasc).GetProjection()
    print(outproj)
    shapegeom=eet.convert_geometry(plotshape,outgeo=False,outproj=outproj)
    xlist=[]
    ylist=[]
    for i in range(0, shapegeom.GetGeometryCount()):
        basegeom=shapegeom.GetGeometryRef(i)
        for j in range(0, basegeom.GetGeometryCount()):
            subgeom=basegeom.GetGeometryRef(j)
            if subgeom is not None:
                for k in range(0,subgeom.GetPointCount()):
                    pt = subgeom.GetPoint(k)
                    xlist.append(pt[0])
                    ylist.append(pt[1])
            else:
                pt = basegeom.GetPoint(j)
                xlist.append(pt[0])
                ylist.append(pt[1])
    nonzeros=np.where(np.array(xlist) != 0)
    xlist,ylist=np.array(xlist)[nonzeros],np.array(ylist)[nonzeros]
    return(xlist,ylist)
'''
def stylize_outputs(xlists,ylists,orderednames,savefig=False, cmap='YlGn'):
    titlecounter=0
    titlelist=['Pre-Treatment CBD (kg/m^3)', 'Post-Treatment CBD (kg/m^3)','Pre-Treatment CC (percent)','Post-Treatment CC (percent)',
              'Pre-Treatment CFA (class)', 'Post-Treatment CFA (class)','Pre-Treatment FLI (Kw/m)', 'Post-Treatment FLI (Kw/m)',
              'Pre-Treatment MWS (mph)', 'Post-Treatment MWS (mph)']
    for fname in orderednames:
        image=np.swapaxes(np.swapaxes(gdal.Open(fname).ReadAsArray(),0,-1),0,-1)
        if 'cfa' in fname:
            image=np.clip(image,1,3)
            cmap='Wistia'
        #find the xmin, ymin, xmax, ymax
        counter = 0
        for xs, ys in zip(xlists,ylists):
            if counter == 0:
                xmin = np.min(xs)
                xmax = np.max(xs)
                ymin = np.min(ys)
                ymax = np.max(ys)
            else:
                xmin = np.min(np.append(xmin,xs))
                xmax = np.max(np.append(xmax,xs))
                ymin = np.min(np.append(ymin,ys))
                ymax = np.max(np.append(ymax,ys))
            counter+=1
        print(xmin,xmax,ymin,ymax)
        plt.xlim((xmin,xmax))
        plt.ylim((ymin,ymax))
        #plot all list of shapes
        plt.imshow(X=image,cmap=cmap,extent=[xmin,xmax,ymin,ymax])
        for xs, ys in zip(xlists,ylists):
            plt.plot(xs,ys,'-',color='blue',linewidth=2)
        plt.title(titlelist[titlecounter])
        print('Mean of area is ' + str(np.mean(image.flatten()[image.flatten() > 0])),
              'Std of area is ' + str(np.std(image.flatten()[image.flatten() > 0])))
        yspace=(ymax-ymin)/4
        plt.yticks([ymin+yspace,ymin+yspace*2,ymax-yspace])
        xspace=(xmax-xmin)/4
        plt.xticks([xmin+xspace,xmin+xspace*2,xmax-xspace])
        cb=plt.colorbar(shrink=.4)
        tick_locator = ticker.MaxNLocator(nbins=4)
        cb.locator = tick_locator
        cb.update_ticks()
        if savefig:
            plt.savefig(orderednames[titlecounter][:-4]+'figure.png',bbox_inches='tight',dpi=300)
        plt.show()
        titlecounter+=1
        plt.close()
            
    
    
        
