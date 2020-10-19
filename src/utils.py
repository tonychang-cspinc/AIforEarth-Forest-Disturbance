import rasterio
import rasterio.warp
from shapely.geometry import Polygon
import geopandas as gpd
import osgeo
import osr
import gdal
import yaml
import os

def add_suffix(filename,suffix):
    ext = os.path.splitext(filename)[1]
    prefix = os.path.splitext(filename)[0]
    final = f'{prefix}_{suffix}{ext}'
    return final
    
def dict_merge(dict1, dict2): 
    res = {**dict1, **dict2} 
    return res 

def write_pickle(filename,data):
    print(f'Creating {filename}')
    f = bz2.open(f'{filename}','wb') 
    pickle.dump(data,f)    
    f.close()
    
def parameter_loader(yamlfile=f'/contents/parameters.yaml'):
    with open(yamlfile) as parameter_yaml:
        params = yaml.safe_load(parameter_yaml)
        wd = params['working_directory'] 
        site = params['site']
        exp_dir = params['export_directory']
        localpath = f'{wd}/{exp_dir}/{site}'
        if not os.path.exists(localpath):
            os.makedirs(localpath) 
        reflcp = params['reference_lcp']
        bbox_name = params['target_bounding_box']
        shp = params['target_shapefile']
        if not reflcp: 
            reflcp = f'{wd}/{site}/{site}.lcp'
        if not bbox_name: 
            bbox_name = f'{wd}/{site}/{site}_bbox.shp'
        if not shp: 
            shp = f'{wd}/{site}/{site}.shp'
        idealpredate = params['pre_treatment_date']
        idealpostdate = params['post_treatment_date']
        max_cloud_thr = params['max_cloud_threshold']
        outbucket = params['gcp_bucket']
        outpath = params['gcp_bucket_subfolder']
    return wd, site, localpath, reflcp, bbox_name, shp, idealpredate, idealpostdate,\
          max_cloud_thr, outbucket, outpath

def make_bounding_box_from_bounds(bounds):
    bbox = [(bounds.left,bounds.top),\
              (bounds.right,bounds.top),\
              (bounds.right,bounds.bottom),\
              (bounds.left,bounds.bottom),\
              (bounds.left,bounds.top)]
    return bbox

def get_points_list(bbox):
    xs = [i[0] for i in bbox]
    ys = [i[1] for i in bbox]
    return(xs,ys)

def make_bbox_polygon(bounds):
    bbox = make_bounding_box(bounds)
    xs,ys = get_point_list(bbox)
    return Polygon(zip(xs,ys))

def generate_bbox_shapefile(ref, dst_file='name', dst_crs={'init':'EPSG:4326'}):
    ds = rasterio.open(ref)
    bbox = make_bounding_box_from_bounds(ds.bounds)
    xs, ys = get_points_list(bbox)
    tfm_x, tfm_y = rasterio.warp.transform(src_crs=ds.crs,\
                                            dst_crs=dst_crs,\
                                            xs=xs, ys=ys)
    bbox_poly = Polygon(zip(tfm_x, tfm_y))
    out_df = gpd.GeoDataFrame({'name':dst_file, 'geometry':[bbox_poly]})
    out_df.crs = dst_crs
    out_df.to_file(dst_file)
    return out_df

def make_warp_parameters(ds):
    ref = gdal.Open(ds)
    res = ref.GetGeoTransform()[1]
    proj = ref.GetProjection()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(proj)   
    return(srs, res)

def make_extent(ds):
    if type(ds) != osgeo.gdal.Dataset:
        ds = gdal.Open(ds)
    gt = ds.GetGeoTransform()
    r = ds.RasterYSize
    c = ds.RasterXSize
    extent = {'xmin':gt[0],'xmax':gt[0]+gt[1]*c,'ymax':gt[3],'ymin':gt[3]+gt[-1]*r}
    ext_list = [extent['xmin'], extent['xmax'], extent['ymin'], extent['ymax']]
    return extent, ext_list
