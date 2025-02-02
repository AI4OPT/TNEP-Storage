import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import sys
import json

simdir = sys.argv[1]
date = sys.argv[2]

df_storage = pd.read_csv(f"../../../{simdir}/output/storage_investments.csv")
df_lines = pd.read_csv(f"../../../{simdir}/output/line_investments.csv")
df_energy = pd.read_csv(f"../../../{simdir}/output/{date}/energy.csv")
production_cols = [col for col in df_energy.columns if col.endswith('_production')]

# Calculate curtailment for each production type
for prod_col in production_cols:
    original_col = prod_col.replace('_production', '')
    curtailment_col = f'Curtailment_{original_col}'
    df_energy[curtailment_col] = df_energy[prod_col] - df_energy[original_col]

# Calculate total curtailment across all production types for each node and hour
curtailment_cols = [f'Curtailment_{col.replace("_production", "")}' for col in production_cols]
df_energy['Total_Curtailment'] = df_energy[curtailment_cols].sum(axis=1)

# Plot for Total Curtailment
fig = px.scatter_mapbox(
    df_energy,
    lat='Lat',
    lon='Lon',
    size='Total_Curtailment',               # Size based on total curtailment
    color='Total_Curtailment',              # Color based on total curtailment
    color_continuous_scale='Burg_r',     # Use a reversed viridis color scale (smaller curtailments are darker)
    center=dict(lat=31.5, lon=-99.5),       # Center the map on Texas
    zoom=4,                                 # Set zoom level
    labels={'Total_Curtailment': 'Total Curtailment (100 MWh)'},  # Label for the color legend
    mapbox_style="carto-positron",          # Map style
)

# Function to add lines to the map
def add_lines(df, color="black", width=0.2):
    traces = []
    for _, row in df.iterrows():
        traces.append(go.Scattermapbox(
            lat=[row['Lat1'], row['Lat2']],
            lon=[row['Lon1'], row['Lon2']],
            mode='lines',
            showlegend=False,
            line=dict(color=color, width=width)
        ))
    return traces

# Add all lines (black, width 0.2)
fig.add_traces(add_lines(df_lines, color="black", width=0.2))

fig.write_image(f"../../../{simdir}/visual/curtailment_{date}.png")



