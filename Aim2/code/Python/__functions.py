"""
Helper functions for aspen intensity/severity work
Python library imports
maxwell.cook@colorado.edu
"""

import gc, time, os, sys, glob
import shutil
import numpy as np
import pandas as pd
import rioxarray as rxr
import matplotlib.pyplot as plt
import pyproj

from itertools import combinations
from collections import Counter
from datetime import datetime
from zoneinfo import ZoneInfo
from shapely.geometry import box
from shapely.geometry import Polygon, MultiPolygon
from rasterstats import zonal_stats

import warnings
warnings.filterwarnings("ignore")  # suppresses annoying geopandas warning


def list_files(path, ext, recursive):
    """
    List files of a specific type in a directory or subdirectories
    """
    if recursive is True:
        return glob.glob(os.path.join(path, '**', '*{}'.format(ext)), recursive=True)
    else:
        return glob.glob(os.path.join(path, '*{}'.format(ext)), recursive=False)


def save_zip(gdf, zip_path, temp_dir):
    """
    Save a GeoDataFrame as a zipped shapefile, optionally using a specific temporary directory.

    Parameters:
    - gdf (GeoDataFrame): The GeoDataFrame to save.
    - zip_path (str): Path to the ZIP archive to create.
    - temp_dir (str, optional): Path to the directory to use for temporary shapefile storage.
                                If None, a temporary directory will be created.

    Returns:
    - str: Path to the ZIP archive.
    """
    try:
        os.makedirs(temp_dir, exist_ok=True)

        # Define the shapefile name and path within the temp_dir
        shp_name = os.path.splitext(os.path.basename(zip_path))[0]
        shp_path = os.path.join(temp_dir, shp_name + '.shp')

        # Save the shapefile
        gdf.to_file(shp_path)

        # Ensure all associated files are present
        base_name = os.path.join(temp_dir, shp_name)
        if not all(os.path.exists(base_name + ext) for ext in ['.shp', '.shx', '.dbf']):
            raise FileNotFoundError("One or more shapefile components are missing.")

        # Compress the shapefile into a ZIP archive
        shutil.make_archive(os.path.splitext(zip_path)[0], 'zip', root_dir=temp_dir, base_dir='.')

        return zip_path

    finally:
        # Clean up the temp directory
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)


def convert_datetime(acq_date, acq_time, zone=None):
    """ Function to convert ACQ_DATE and ACQ_TIME to a datetime object in UTC """
    # Ensure ACQ_TIME is in HHMM format
    acq_time = str(acq_time).zfill(4)  # Pad to 4 digits if necessary

    # Convert acq_date to a string in the format '%Y-%m-%d'
    if isinstance(acq_date, str):
        acq_date = datetime.strptime(acq_date, '%m/%d/%Y')  # Adjust based on input format

    acq_date_str = acq_date.strftime('%Y-%m-%d')
    dt = datetime.strptime(acq_date_str + acq_time, '%Y-%m-%d%H%M')

    # Localize the datetime object to the specified timezone
    dt_utc = dt.replace(tzinfo=ZoneInfo("UTC"))  # localize the UTC timezone
    if zone is not None:
        dt_zone = dt_utc.astimezone(ZoneInfo(zone))
        return dt_zone
    else:
        return dt_utc


def get_coords(geom, buffer, crs='EPSG:4326'):
    """ Returns the bounding box coordinates for a given geometry(ies) and buffer """
    _geom = geom.copy()
    _geom['geometry'] = _geom.geometry.buffer(buffer)
    bounds = _geom.to_crs(crs).unary_union.envelope
    coords = list(bounds.exterior.coords)

    # Calculate extent
    min_lon = min([p[0] for p in coords])
    max_lon = max([p[0] for p in coords])
    min_lat = min([p[1] for p in coords])
    max_lat = max([p[1] for p in coords])

    extent = [min_lon, max_lon, min_lat, max_lat]

    del _geom, bounds
    return coords, extent


def compute_band_stats(geoms, image_da, id_col, attr=None, stats=None, ztype='categorical'):
    """
    Function to compute band statistics for geometries and a single raster band.
    Args:
        geoms: the geometries for which to calculate zonal statistics
        image_da: categorical raster image array
        id_col: the unique identifier for geometries
        attr: the attribute to calculate (example, 'CBH_') for naming
        stats: statistics to calculate (if continuous input data) as list of strings
        ztype: whether to treat raster data as categorical or continuous
    """
    affine = image_da.rio.transform()
    nodataval = image_da.rio.nodata
    arr = image_da.values

    if ztype == 'categorical':

        if attr is None:
            attr = 'evt'

        zs = zonal_stats(
            vectors=geoms[[id_col, 'geometry']],
            raster=arr,
            affine=affine,
            nodata=nodataval,
            categorical=True,
            all_touched=True,
            geojson_out=True
        )

        # Extract the results (properties)
        stats_df = pd.DataFrame(zs)
        stats_df[id_col] = stats_df['properties'].apply(lambda x: x.get(id_col))
        stats_df['properties'] = stats_df['properties'].apply(
            lambda x: {key: val for key, val in x.items() if key != id_col})
        stats_df['props_list'] = stats_df['properties'].apply(lambda x: list(x.items()))

        # Explode the properties to column
        props = stats_df.explode('props_list').reset_index(drop=True)
        props[[attr,'count']] = pd.DataFrame(props['props_list'].tolist(), index=props.index)

        # Handle NaN values
        props.dropna(subset=[attr], inplace=True)

        # Tidy the columns.
        props[attr] = props[attr].astype(int)
        props = props[[id_col,attr,'count']].reset_index(drop=True)

        # Calculate the total pixels and percent cover
        total_pixels = props.groupby(props[id_col])['count'].transform('sum')
        props['total_pixels'] = total_pixels
        props['pct_cover'] = (props['count'] / props['total_pixels']) * 100

        del arr, stats, stats_df  # clean up
        gc.collect()

        return props

    elif ztype == 'continuous':
        # Make sure 'stats' is defined
        if stats is None:
            print("! Please provide list of statistics to calculate !")
            return None
        else:
            zs = zonal_stats(
                vectors=geoms[[id_col, 'geometry']],
                raster=arr,
                affine=affine,
                nodata=nodataval,
                stats=stats,
                categorical=False,
                all_touched=True,
                geojson_out=True
            )

            # Extract the dataframe
            stats_df = pd.DataFrame(zs)
            stats_df[id_col] = stats_df['properties'].apply(lambda x: x.get(id_col))
            for stat in stats:
                stats_df[stat] = stats_df['properties'].apply(lambda x: x.get(stat))
                stats_df.rename(columns={stat: f'{attr}_{stat}'}, inplace=True)

            # Tidy the columns
            cols_to_keep = [id_col] + [f'{attr}_{stat}' for stat in stats]
            stats_df = stats_df[cols_to_keep]

            return stats_df


def create_bounds(gdf, buffer=None, method='bounds'):
    """
    Calculate a bounding rectangle for a given geometry and buffer
    Args:
        gdf: perimeter geometry
        buffer: buffer distance to be applied
        method: one of ['bounds','convex_hull','exact']
    """
    if method == 'bounds':
        geom = gdf.geometry.apply(lambda geo: box(*geo.bounds))
        if buffer is not None:
            geom = geom.buffer(buffer)
        # Apply the new geometry
        gdf_ = gdf.copy()
        gdf_.geometry = geom.geometry.apply(
            lambda geo: Polygon(geo) if geo.geom_type == 'Polygon' else MultiPolygon([geo]))
        return gdf_
    elif method == 'convex_hull':
        gdf_ = gdf.copy()
        gdf_['geometry'] = gdf_.geometry.convex_hull
        if buffer is not None:
            gdf_['geometry'] = gdf_.geometry.buffer(buffer)
        return gdf_
    elif method == 'exact':
        if buffer is not None:
            gdf_ = gdf.geometry.buffer(buffer)
            return gdf_
        else:
            gdf_ = gdf.copy()
            gdf_.geometry = gdf_.geometry.buffer(buffer)
            return gdf_


def weighted_variance(values, weights):
    """ Calculate weighted variance. """
    average = np.average(values, weights=weights)
    variance = np.average((values - average) ** 2, weights=weights)
    return variance


def find_best_match(perim, neighbors, max_size_diff):
    """
    Identifies the 'best match' based on spatial distance and size difference
    Args:
        perim: the polygon geometries
        neighbors: nearest neighbors identified (see 'find_nearest' function)
        max_size_diff: the maximum allowed size difference (in % difference)
    """
    # Get the fire size from the perimeter data
    perim_size = perim['Final_Acres']

    # Initialize best score and match
    best_score = float('inf')
    best_match = None

    for _, point in neighbors.iterrows():
        # Calculate the size difference
        ics_size = point.get('FINAL_ACRES', np.nan)

        if perim_size != 0:
            size_diff = abs((ics_size - perim_size) / perim_size) * 100
        else:
            size_diff = float('inf')  # If perim_size is 0, treat it as infinite difference

        if max_size_diff is not None:
            if size_diff > max_size_diff:
                continue  # Skip if size difference exceeds max_size_diff

        # Calculate the spatial distance (assuming it's precomputed)
        spatial_dist = point.get('spatial_dist', 0)

        # Composite score (you can adjust the weights as needed)
        score = spatial_dist + size_diff

        # Check if this is the best match
        if score < best_score:
            best_score = score
            best_match = point

    return best_match


def find_nearest(perims, points, nn, max_dist=50000, max_size_diff=150):
    """
    Finds the nearest points based on spatial proximity, size, and temporal alignment.
    """

    out_nns = []  # storing the resulting nearest neighbors for each perimeter
    no_matches = []  # to store fires with no matches

    for _, perim in perims.iterrows():
        fire_id = perim['Fire_ID']
        perim_geom = perim.geometry
        fire_year = perim['Fire_Year']

        # Filter incident points to the fire year (filtered once per perimeter)
        inci_points = points[points['START_YEAR'] == fire_year]

        # Early check if no points match the fire year
        if inci_points.empty:
            print(f"No matching points for fire year and fire id: {fire_year} / {fire_id}")
            no_matches.append(perim.to_frame().T)  # Append as DataFrame
            continue

        # Calculate distances from the fire perimeter to the incident points
        distances = inci_points.geometry.apply(lambda x: perim_geom.distance(x))

        # Filter by the maximum distance if provided
        if max_dist is not None:
            inci_points = inci_points[distances <= max_dist]
            distances = distances[distances <= max_dist]

        # Check if there are still points left after filtering
        if inci_points.empty:
            no_matches.append(perim.to_frame().T)  # Convert row to DataFrame and append
            continue

        # Sort by distance and retain the nearest NN points
        nearest_points = inci_points.iloc[distances.argsort()[:nn]].copy()

        # Calculate the best match based on size and distance
        best_ = find_best_match(perim, nearest_points, max_size_diff=max_size_diff)

        if best_ is not None:
            best_['Fire_ID'] = perim['Fire_ID']
            best_['Fire_Name'] = perim['Fire_Name']
            best_['Final_Acres'] = perim['Final_Acres']
            best_['Source'] = perim['Source']
            best_['Start_Date'] = perim['Start_Date']
            best_['Aspen_Pct'] = perim['pct_aspen']
            out_nns.append(best_.to_frame().T)  # Convert best match to DataFrame before appending

    # Concatenate the no_matches
    print(f"There were [{len(no_matches)}/{len(perims)}] fires with no matches.")
    if len(no_matches) > 0:
        no_matches = pd.concat(no_matches, ignore_index=True)
    else:
        no_matches = pd.DataFrame()

    # Concatenate the matches
    if len(out_nns) > 0:
        out_nns = pd.concat(out_nns, ignore_index=True)
    else:
        out_nns = pd.DataFrame()

    return out_nns, no_matches


def monitor_export(task, timeout=30):
    """ Monitors EE export task """
    while task.active():
        print('Waiting for export to finish..\n\tPatience young padawan.')
        time.sleep(timeout)  # Check every 30 seconds

    # Get the status of the task
    status = task.status()

    # Check if the task failed or succeeded
    if status['state'] == 'COMPLETED':
        print("Export completed successfully !!!!")
    elif status['state'] == 'FAILED':
        print(f"Export failed! Bummer. Reason: {status.get('error_message', 'Unknown error')}")
    else:
        print(f"Export ended with state: {status['state']}")


def get_spp_coo(df, spps_list=None, grid_col='grid_index', sp_col='fortypnm_gp', wt_col=None):
    """
    Analyzes species co-occurrence within grid cells.
    Parameters:
        df (pd.DataFrame): Input DataFrame containing species and grid information.
        spps_list (list): List of species to filter and analyze.
        grid_col (str): Column name representing the spatial unit (default: 'grid_index').
        sp_col (str): Column name representing the species (default: 'fortypnm_gp').
        wt_col (str): Column to use for weighting by abundance or dominance
    Returns:
        pd.DataFrame: DataFrame with co-occurrence counts and percentages.
    """

    # Filter the DataFrame to only include relevant species
    df[sp_col] = df[sp_col].str.lower()  # force to lower case
    if spps_list is not None:
        df_ = df[df[sp_col].isin(spps_list)].copy()
    else:
        df_ = df.copy()
    print(f"\nSpecies occurrence counts:\n{df_[sp_col].value_counts()}\n")

    # If not weight is provided, assume binary presence
    if wt_col is None:
        df_['weight'] = 1
    else:
        df_['weight'] = df_[wt_col]

    # Create species-grid abundance pivot table
    spp_grid = df_.pivot_table(
        index=grid_col, columns=sp_col, values='weight', aggfunc='sum', fill_value=0
    )

    # Compute co-occurrence by multiplying the species-grid matrix with its transpose
    coo_mat = spp_grid.T.dot(spp_grid)

    # Compute co-occurrence percentage: (shared grids) / (grids where either species appears)
    species_totals = np.diag(coo_mat)  # Total occurrences per species
    coo_pct = coo_mat.div(species_totals[:, None] + species_totals - coo_mat, axis=0)

    # Set self-co-occurrence to 1
    np.fill_diagonal(coo_pct.values, 1)

    return coo_mat, coo_pct


def resample_bilinear(in_img, scale_factor, crs='EPSG:5070', match_img=None):
    """
    :param in_img:
    :param to_img:
    :param scale_factor:
    :param proj4:
    :return:
    """

    resampled = in_img.rio.reproject(
        in_img.rio.crs,
        shape=(
            int(in_img.rio.height * scale_factor),
            int(in_img.rio.width * scale_factor),
        ),
        resampling=Resampling.bilinear
    )

    if match_img is not None:
        match_img = rxr.open_rasterio(match_img, masked=True, cache=False).squeeze()
        resampled = resampled.rio.reproject_match(match_img)
        del match_img
        return resampled
    else:
        return resampled


def rasterize_it(zones, ref_img, zone_col, open=False, crs='EPSG:5070'):
    """
    :param shp:
    :param crs:
    :param to_img:
    :return:
    """

    if open is True:
        ref = rxr.open_rasterio(ref_img, masked=True, lock=False, chunks=True).squeeze()
    else:
        ref = ref_img

    # make sure the CRS matches
    zones = zones.to_crs(crs)

    # Using GeoCube to rasterize the Vector
    rastered = make_geocube(
        vector_data=zones,
        measurements=[zone_col],
        resolution=(-1000, 1000),
        fill=np.nan
    )

    # Convert from Dataset to DataArray (extract first variable)
    arr = rastered[zone_col]
    # Ensure it matches the reference raster spatially
    arr = arr.rio.write_crs(crs)  # Assign CRS
    matched = arr.rio.reproject_match(ref, resampling=Resampling.bilinear)

    del rastered, arr

    return matched


def plot_raster(raster, title="", cmap="viridis", legend_lab="", save_file=False, out_png=None):
    """
    Plot a raster dataset using matplotlib.

    :param raster: xarray.DataArray containing the raster data.
    :param title: Title of the plot.
    :param cmap: Colormap for visualization.
    """

    # Get spatial bounds
    bounds = raster.rio.bounds()
    xmin, ymin, xmax, ymax = bounds

    fig, ax = plt.subplots(figsize=(10, 6))

    # Mask out no-data values for better visualization
    im = ax.imshow(
        raster.squeeze(),
        cmap=cmap,
        extent=[xmin, xmax, ymin, ymax],  # Keep original projection extent
        origin="upper"
    )

    # Convert tick locations from EPSG:5070 to EPSG:4326
    transformer = pyproj.Transformer.from_crs("EPSG:5070", "EPSG:4326", always_xy=True)

    def transform_ticks(ticks, axis="x"):
        if axis == "x":
            lat, lon = transformer.transform(ticks, np.full_like(ticks, ymin))
        else:
            lat, lon = transformer.transform(np.full_like(ticks, xmin), ticks)
        return lat if axis == "x" else lon

    # Set new x-axis ticks
    xticks = np.linspace(xmin, xmax, num=6)
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{lon:.2f}" for lon in transform_ticks(xticks, "x")])

    # Set new y-axis ticks
    yticks = np.linspace(ymin, ymax, num=6)
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{lat:.2f}" for lat in transform_ticks(yticks, "y")])

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.7)
    cbar.set_label(legend_lab)

    # Labels and title
    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

    # save if specific
    if save_file is True:
        if out_png is not None:
            plt.savefig(out_png, dpi=500, bbox_inches='tight')
        else:
            print("File not saved, please specify an output path.")

    plt.show()


def download_file(img_url, outname, chunk_size=1024):
    try:
        with requests.get(img_url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(outname, 'wb') as f:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if chunk:  # filter out keep-alive chunks
                        f.write(chunk)
    except Exception as e:
        print(f"Error downloading {img_url}: {e}")
        raise


def get_image_service_array(url, ply, out_prefix, res=30, outSR=3857, outdir=None, cleanup=True, plot=False):
    """
    Downloads raster imagery from an ArcGIS REST ImageService for a given polygon.

    Parameters:
    ----------
    url : str
        REST endpoint of the ImageServer
    ply : GeoDataFrame or GeoSeries
        Area of interest (single polygon or feature)
    out_prefix : str
        Prefix for naming downloaded image tiles
    res : int
        Pixel resolution (in meters)
    outSR : str or int
        Output spatial reference (e.g., 5070); default uses service SR
    output_dir : str or None
        Directory to save tiles. If None, uses a temporary directory.
    cleanup : bool
        Whether to delete files after processing (if using temp directory)

    Returns:
    --------
    xarray.DataArray
        Combined raster from all tiles (or single tile)
    """

    # Query the metadata from the service (url)
    meta = requests.get(url + '?f=pjson').json()

    if 'spatialReference' not in meta:
        # Try parent layer for spatial reference
        print(f"[WARNING] No spatialReference found for {url}, checking parent...")
        parent_url = '/'.join(url.strip('/').split('/')[:-1])
        parent_meta = requests.get(parent_url + '?f=pjson').json()
        if 'spatialReference' not in parent_meta:
            raise ValueError(f"No spatial reference found in metadata for {url} or its parent.")
        spr = parent_meta['spatialReference']
    else:
        spr = meta['spatialReference']

    epsg = spr.get('latestWkid') or spr.get('wkid')
    if epsg is None:
        raise ValueError("Unable to determine EPSG code from spatialReference.")

    # Use default max dimensions if missing
    max_w = meta.get('maxImageWidth', 1024)
    max_h = meta.get('maxImageHeight', 1024)

    # Reproject the region of interest to the correct CRS
    # Extract the bounding geometry (optional buffer)
    if ply.crs is None:
        raise ValueError("Input geometry must have a defined CRS.")
    ply_proj = ply.to_crs(epsg)
    xmin, ymin, xmax, ymax = ply_proj.total_bounds

    # Computing the tiling extents
    wcells = int((xmax - xmin) / res)
    hcells = int((ymax - ymin) / res)
    tile_w = min(wcells, max_w)
    tile_h = min(hcells, max_h)
    wcells_l = np.arange(0, wcells, tile_w)
    hcells_l = np.arange(0, hcells, tile_h)

    # Prepare output directory
    temp_dir = None
    if outdir is None:
        temp_dir = tempfile.mkdtemp()
        outdir = temp_dir
    os.makedirs(outdir, exist_ok=True)

    # Download the tiles
    rasters = []
    tile = 1
    for w in wcells_l:
        for h in hcells_l:
            xmin2 = xmin + w * res
            xmax2 = min(xmin + (w + tile_w) * res, xmax)
            ymin2 = ymin + h * res
            ymax2 = min(ymin + (h + tile_h) * res, ymax)

            qry = url + '/exportImage?'
            parm = {
                'f': 'json',
                'where': 'Year = 2020',
                'bbox': f"{xmin2},{ymin2},{xmax2},{ymax2}",
                'size': f"{tile_w},{tile_h}",
                'imageSR': outSR,
                'format': 'tiff'
            }
            resp = requests.get(qry, parm, timeout=60)
            if resp.status_code == 200:
                result = resp.json()
                img_url = result.get('href') or result.get('imageUrl')

                if img_url is None:
                    continue  # Skip to next tile

                out_fp = os.path.join(outdir, f"{out_prefix}_{tile}.tif")
                download_file(img_url, out_fp)
                rasters.append(rxr.open_rasterio(out_fp, masked=True).squeeze())
                tile += 1
            else:
                print(f"[WARNING] Tile {tile} failed: HTTP {resp.status_code}")

    if not rasters:
        raise ValueError("No tiles were successfully downloaded.")

    # 6. Merge tiles
    result = rasters[0] if len(rasters) == 1 else xr.combine_by_coords(rasters)

    # 7. Cleanup
    if cleanup and temp_dir is not None:
        shutil.rmtree(temp_dir)

    # 8. Plot
    if plot is True:
        result.plot(cmap='viridis')
        # plt.title(f"({key})")
        plt.show()

    return result