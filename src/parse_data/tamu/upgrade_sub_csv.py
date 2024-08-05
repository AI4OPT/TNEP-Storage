import geopandas as gpd
import pandas as pd
from shapely.geometry import Point

GEOJSON_PATH = "../../../data/geojson/world.geojson"
SUB_CSV_PATH = "../../tamu/base_grid/sub.csv"

df = pd.read_csv(SUB_CSV_PATH)

gdf_points = gpd.GeoDataFrame(df, geometry=[Point(xy) for xy in zip(df.lon, df.lat)], crs="EPSG:4326")
gdf_polygons = gpd.read_file(GEOJSON_PATH)
gdf_polygons = gdf_polygons[gdf_polygons['countryName'] == "United States"]


counter = 0

# Define a function to find the zoneName for each point
def find_zone_name(point, polygons):
    global counter
    if (counter % 10000 == 0):
        print(counter)
    counter += 1
    
    for _, poly in polygons.iterrows():
        if point.within(poly.geometry):
            return poly['zoneName']
    return find_closest_zone_name(point, polygons)

def find_closest_zone_name(point, polygons):
    min_distance = float('inf')
    closest_zone_name = None
    
    # Iterate through each polygon to find the closest one
    for _, poly in polygons.iterrows():
        distance = point.distance(poly.geometry)
        if distance < min_distance:
            min_distance = distance
            closest_zone_name = poly['zoneName']
            
    return closest_zone_name


# Apply the function to find the zoneName for each point
gdf_points['zoneName'] = gdf_points.apply(lambda x: find_zone_name(x.geometry, gdf_polygons), axis=1)

# Convert back to a regular DataFrame and drop the geometry column if no longer needed
df_final = pd.DataFrame(gdf_points.drop(columns=['geometry']))

# Display the first few rows to verify
print(df_final.head())

df_final.to_csv('../data/geojson/upgraded_sub.csv', index=False)