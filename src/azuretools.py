import tempfile
import warnings
import urllib
import shutil
import os
import io
import requests

# Workaround for a problem in older rasterio versions
os.environ["CURL_CA_BUNDLE"] = "/etc/ssl/certs/ca-certificates.crt" 

# Less standard, but still pip- or conda-installable
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import gdal
import osr
import rasterio
import rtree
import shapely
import pickle
import math

# pip install progressbar2, not progressbar
import progressbar

from geopy.geocoders import Nominatim
from rasterio.windows import Window 
from tqdm import tqdm
from azure.storage.blob import ContainerClient

crs = "EPSG:4326"

# Storage locations are documented at http://aka.ms/ai4edata-naip
blob_root = 'https://naipblobs.blob.core.windows.net/naip'

index_files = ["tile_index.dat", "tile_index.idx", "tiles.p"]
index_blob_root = 'https://naipblobs.blob.core.windows.net/naip-index/rtree/'
# Storage locations are documented at http://aka.ms/ai4edata-nasadem
nasadem_account_name = 'nasadem'
nasadem_container_name = 'nasadem-nc'
nasadem_account_url = 'https://' + nasadem_account_name + '.blob.core.windows.net'
nasadem_blob_root = nasadem_account_url + '/' + nasadem_container_name + '/v001/'

# A full list of files is available at:
#
# https://nasadem.blob.core.windows.net/nasadem-nc/v001/index/file_list.txt
nasadem_file_index_url = nasadem_blob_root + 'index/nasadem_file_list.txt'
nasadem_content_extension = '.nc'
nasadem_file_prefix = 'NASADEM_NC_'

temp_dir = os.path.join(tempfile.gettempdir(),'naip')
os.makedirs(temp_dir,exist_ok=True)

# Storage locations are documented at http://aka.ms/ai4edata-hls
hls_container_name = 'hls'
hls_account_name = 'hlssa'
hls_account_url = 'https://' + hls_account_name + '.blob.core.windows.net/'
hls_blob_root = hls_account_url + hls_container_name

# This file is provided by NASA; it indicates the lat/lon extents of each
# hls tile.
#
# The file originally comes from:
#
# https://hls.gsfc.nasa.gov/wp-content/uploads/2016/10/S2_TilingSystem2-1.txt
#
# ...but as of 8/2019, there is a bug with the column names in the original file, so we
# access a copy with corrected column names.

hls_tile_extents_url = 'https://ai4edatasetspublicassets.blob.core.windows.net/assets/S2_TilingSystem2-1.txt?st=2019-08-23T03%3A25%3A57Z&se=2028-08-24T03%3A25%3A00Z&sp=rl&sv=2018-03-28&sr=b&sig=KHNZHIJuVG2KqwpnlsJ8truIT5saih8KrVj3f45ABKY%3D'

# Load this file into a table, where each row is:
#
# Tile ID, Xstart, Ystart, UZ, EPSG, MinLon, MaxLon, MinLon, MaxLon
#print('Reading tile extents...')
s = requests.get(hls_tile_extents_url).content
hls_tile_extents = pd.read_csv(io.StringIO(s.decode('utf-8')),delimiter=r'\s+')
#print('Read tile extents for {} tiles'.format(len(hls_tile_extents)))

# Read-only shared access signature (SAS) URL for the hls container
hls_sas_token = 'st=2019-08-07T14%3A54%3A43Z&se=2050-08-08T14%3A54%3A00Z&sp=rl&sv=2018-03-28&sr=c&sig=EYNJCexDl5yxb1TxNH%2FzILznc3TiAnJq%2FPvCumkuV5U%3D'

hls_container_client = ContainerClient(account_url=hls_account_url, 
container_name=hls_container_name,
credential=None)


class DownloadProgressBar():
    """
    https://stackoverflow.com/questions/37748105/how-to-use-progressbar-module-with-urlretrieve
    """
    
    def __init__(self):
        self.pbar = None

    def __call__(self, block_num, block_size, total_size):
        if not self.pbar:
            #self.pbar = progressbar.ProgressBar(max_value=total_size)
            self.pbar = progressbar.ProgressBar(maxval=total_size)
            self.pbar.start()
            
        downloaded = block_num * block_size
        if downloaded < total_size:
            self.pbar.update(downloaded)
        else:
            self.pbar.finish()
            
##########################        
#####NAIP tools###########        
##########################        

class NAIPTileIndex:
    """
    Utility class for performing NAIP tile lookups by location.
    """
    
    tile_rtree = None
    tile_index = None
    base_path = None
    
    def __init__(self, base_path=None):
        
        if base_path is None:
            
            base_path = temp_dir
            os.makedirs(base_path,exist_ok=True)
            
            for file_path in index_files:
                download_url(index_blob_root + file_path, base_path + '/' + file_path,
                             progress_updater=DownloadProgressBar())
                
        self.base_path = base_path
        self.tile_rtree = rtree.index.Index(base_path + "/tile_index")
        self.tile_index = pickle.load(open(base_path  + "/tiles.p", "rb"))
      
    
    def lookup_tile(self, lat, lon):
        """"
        Given a lat/lon coordinate pair, return the list of NAIP tiles that contain
        that location.

        Returns an array containing [mrf filename, idx filename, lrc filename].abs        """

        point = shapely.geometry.Point(float(lon),float(lat))
        intersected_indices = list(self.tile_rtree.intersection(point.bounds))

        intersected_files = []
        tile_intersection = False

        for idx in intersected_indices:

            intersected_file = self.tile_index[idx][0]
            intersected_geom = self.tile_index[idx][1]
            if intersected_geom.contains(point):
                tile_intersection = True
                intersected_files.append(intersected_file)

        if not tile_intersection and len(intersected_indices) > 0:
            print('''Error: there are overlaps with tile index, 
                      but no tile completely contains selection''')   
            return None
        elif len(intersected_files) <= 0:
            print("No tile intersections")
            return None
        else:
            return intersected_files
        
            
def download_url(url, destination_filename=None, progress_updater=None,\
        force_download=False, quiet=True):
    """
    Download a URL to a temporary file
    """
    
    # This is not intended to guarantee uniqueness, we just know it happens to guarantee
    # uniqueness for this application.
    if destination_filename is None:
        url_as_filename = url.replace('://', '_').replace('/', '_')    
        destination_filename = \
            os.path.join(temp_dir,url_as_filename)
    if (not force_download) and (os.path.isfile(destination_filename)):
        if not quiet:
            print('Bypassing download of already-downloaded file {}'.format(os.path.basename(url)))
        return destination_filename
    print('Downloading file {} to {}'.format(os.path.basename(url),destination_filename),end='')
    urllib.request.urlretrieve(url, destination_filename, progress_updater)  
    assert(os.path.isfile(destination_filename))
    nBytes = os.path.getsize(destination_filename)
    print('...done, {} bytes.'.format(nBytes))
    return destination_filename
    

def display_naip_tile(filename):
    """
    Display a NAIP tile using rasterio.
    
    For .mrf-formatted tiles (which span multiple files), 'filename' should refer to the 
    .mrf file.
    """
    
    # NAIP tiles are enormous; downsize for plotting in this notebook
    dsfactor = 10
    
    with rasterio.open(filename) as raster:

        # NAIP imagery has four channels: R, G, B, IR
        #
        # Stack RGB channels into an image; we won't try to render the IR channel
        #
        # rasterio uses 1-based indexing for channels.
        h = int(raster.height/dsfactor)
        w = int(raster.width/dsfactor)
        print('Resampling to {},{}'.format(h,w))
        r = raster.read(1, out_shape=(1, h, w))
        g = raster.read(2, out_shape=(1, h, w))
        b = raster.read(3, out_shape=(1, h, w))        
    
    rgb = np.dstack((r,g,b))
    fig = plt.figure(figsize=(7.5, 7.5), dpi=100, edgecolor='k')
    plt.imshow(rgb)
    raster.close()
    
    
def get_coordinates_from_address(address):
    """
    Look up the lat/lon coordinates for an address.
    """
    
    geolocator = Nominatim(user_agent="NAIP")
    location = geolocator.geocode(address)
    print('Retrieving location for address:\n{}'.format(location.address))
    return location.latitude, location.longitude

##############################        
######NASADEM tools###########        
##############################        


def lat_lon_to_nasadem_tile(lat,lon):
    """
    Get the NASADEM file name for a specified latitude and longitude
    """

    # A tile name looks like:
    #
    # NASADEM_NUMNC_n00e016.nc
    #
    # The translation from lat/lon to that string is represented nicely at:
    #
    # https://dwtkns.com/srtm30m/

    # Force download of the file list
    get_nasadem_file_list()

    ns_token = 'n' if lat >=0 else 's'
    ew_token = 'e' if lon >=0 else 'w'

    lat_index = abs(math.floor(lat))
    lon_index = abs(math.floor(lon))

    lat_string = ns_token + '{:02d}'.format(lat_index)
    lon_string = ew_token + '{:03d}'.format(lon_index)

    filename =  nasadem_file_prefix + lat_string + lon_string + \
        nasadem_content_extension

    if filename not in nasadem_file_list:
        print('Lat/lon {},{} not available'.format(lat,lon))
        filename = None

    return filename


def get_nasadem_file_list():
    """
    Retrieve the full list of NASADEM tiles
    """

    global nasadem_file_list
    if nasadem_file_list is None:
        nasadem_file = download_url(nasadem_file_index_url)
        with open(nasadem_file) as f:
            nasadem_file_list = f.readlines()
            nasadem_file_list = [f.strip() for f in nasadem_file_list]
            nasadem_file_list = [f for f in nasadem_file_list if \
                                     f.endswith(nasadem_content_extension)]
    return nasadem_file_list

def lat_lon_to_nasadem_tile(lat,lon,current_dem_list=None):
    """
    Get the NASADEM file name for a specified latitude and longitude
    """

    # A tile name looks like:
    #
    # NASADEM_NUMNC_n00e016.nc
    #
    # The translation from lat/lon to that string is represented nicely at:
    #
    # https://dwtkns.com/srtm30m/

    # Force download of the file list
    nasadem_file_list = get_nasadem_file_list(current_dem_list)

    ns_token = 'n' if lat >=0 else 's'
    ew_token = 'e' if lon >=0 else 'w'

    lat_index = abs(math.floor(lat))
    lon_index = abs(math.floor(lon))

    lat_string = ns_token + '{:02d}'.format(lat_index)
    lon_string = ew_token + '{:03d}'.format(lon_index)

    filename =  nasadem_file_prefix + lat_string + lon_string + \
        nasadem_content_extension

    if filename not in nasadem_file_list:
        print('Lat/lon {},{} not available'.format(lat,lon))
        filename = None

    return filename

def get_nasadem_file_list(current_dem_list):
    """
    Retrieve the full list of NASADEM tiles
    """
    if current_dem_list is None:
        nasadem_file = download_url(nasadem_file_index_url)
        with open(nasadem_file) as f:
            nasadem_file_list = f.readlines()
            nasadem_file_list = [f.strip() for f in nasadem_file_list]
            nasadem_file_list = [f for f in nasadem_file_list if \
                f.endswith(nasadem_content_extension)]
    return nasadem_file_list

##############################        
##########HLS tools###########        
##############################        

def get_hls_tile(blob_url):
    """
    Given a URL pointing to an HLS image in blob storage, load that image via GDAL
    and return both data and metadata.
    """    

    formatted_gdal_bloburl='/{}/{}'.format('vsicurl',blob_url)

    tile_open = gdal.Open(formatted_gdal_bloburl)
    data = tile_open.GetRasterBand(1)
    ndv,xsize,ysize = data.GetNoDataValue(),tile_open.RasterXSize,tile_open.RasterYSize

    projection = osr.SpatialReference()
    projection.ImportFromWkt(tile_open.GetProjectionRef())

    datatype = data.DataType
    datatype = gdal.GetDataTypeName(datatype)  
    data_array = data.ReadAsArray()

    return ndv,xsize,ysize,projection,data_array


def list_available_tiles(prefix):
    """
    List all blobs in an Azure blob container matching a prefix.  

    We'll use this to query tiles by location and year.
    """

    files = []
    generator = hls_container_client.list_blobs(name_starts_with=prefix)
    for blob in generator:
        files.append(blob.name)
    return files


def lat_lon_to_hls_tile_id(lat,lon):
    """
    Get the hls tile ID for a given lat/lon coordinate pair
    """  
    found_matching_tile = False

    for i_row,row in hls_tile_extents.iterrows():
        found_matching_tile = lat >= row.MinLat and lat <= row.MaxLat \
                and lon >= row.MinLon and lon <= row.MaxLon
        if found_matching_tile:
            break

    if not found_matching_tile:
        return None
    else:
        return row.TilID

