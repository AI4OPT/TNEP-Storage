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

df_energy['Imbalance_Size'] = df_energy['Energy_Imbalance'].apply(lambda x: abs(x) if x < 0 else 0)
df_energy = df_energy[df_energy['Imbalance_Size'] > 0]

fig = px.scatter_mapbox(
    df_energy,
    lat='Lat',
    lon='Lon',
    size='Imbalance_Size',               # Size based on solar capacity
    color='Imbalance_Size',              # Color based on solar capacity
    color_continuous_scale='solar',    # Use viridis color scale
    center=dict(lat=31.5, lon=-99.5),    # Center the map
    zoom=4,                              # Set zoom level
    labels={'Imbalance_Size': 'Unserved Energy (100 MW)'},
    mapbox_style="carto-positron",       # Map style
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

fig.write_image(f"../../../{simdir}/visual/ue_{date}.png")