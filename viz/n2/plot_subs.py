import plotly.graph_objects as go
import pandas as pd
import plotly.express as px
import numpy as np

SUB_PATH = "../../data/geojson/upgraded_sub.csv"
df = pd.read_csv(SUB_PATH)

# Get unique zone names and sort them alphabetically
unique_zones = sorted(df['zoneName'].unique())

# Calculate the number of zones
num_zones = len(unique_zones)
print(num_zones)

# Use Plotly Express continuous color scale
color_scale = px.colors.sequential.Viridis  # You can choose any other continuous color scale

# Create a color mapping
color_mapping = {zone: color_scale[i * len(color_scale) // num_zones] for i, zone in enumerate(unique_zones)}

# Apply color mapping to create a new column 'marker_color'
df['marker_color'] = df['zoneName'].map(color_mapping)

fig = go.Figure()

fig.add_trace(go.Scattergeo(
    lon=df['lon'],
    lat=df['lat'],
    mode='markers',
    marker=dict(
        color=df['marker_color'],  # Use the new 'marker_color' column
        size=1,
        colorscale='Viridis',  # Specify the color scale
        cmin=0,
        cmax=num_zones - 1,
    )
))

fig.update_layout(
    title='Nodes',
    geo_scope='usa',
)

fig.write_image("../../data/geojson/zones.png")
