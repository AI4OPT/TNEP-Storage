import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import sys
import json

simdir = sys.argv[1]
cand_file = sys.argv[2]
df_storage = pd.read_csv("../../../{}/output/storage_investments.csv".format(simdir))
df_lines = pd.read_csv("../../../{}/output/line_investments.csv".format(simdir))

with open(f"../../../{cand_file}", "r") as file:
    cands = json.load(file)

df_storage = df_storage[df_storage['Node_Index'].isin([int(c) for c in cands])]

fig = px.scatter_mapbox(
    df_storage, 
    lat='Lat', 
    lon='Lon',
    center=dict(lat=31.5, lon=-99.5), 
    zoom=4,
    mapbox_style="carto-positron",
)

# Add filtered points with fixed size and red color using go.Scattermapbox
fig.add_trace(go.Scattermapbox(
    lat=df_storage['Lat'],
    lon=df_storage['Lon'],
    mode='markers',
    marker=dict(size=10, color='red', opacity=0.8),
    name='Candidate Nodes'
))

# Function to add lines to the map
def add_lines(df, color, width_func):
    traces = []
    for _, row in df.iterrows():
        traces.append(go.Scattermapbox(
            lat=[row['Lat1'], row['Lat2']],
            lon=[row['Lon1'], row['Lon2']],
            mode='lines',
            showlegend=False,
            line=dict(color=color, width=width_func(row))
        ))
    return traces

# Add regular lines (thin black lines)
fig.add_traces(add_lines(df_lines, "black", lambda row: 0.2))

fig.write_image(f"../../../{simdir}/visual/storage_candidates.png")