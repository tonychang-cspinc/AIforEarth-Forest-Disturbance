import tempfile
import warnings
import urllib
import shutil
import os

# Workaround for a problem in older rasterio versions
os.environ["CURL_CA_BUNDLE"] = "/etc/ssl/certs/ca-certificates.crt" 

# Less standard, but still pip- or conda-installable
import matplotlib.pyplot as plt
import numpy as np
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

