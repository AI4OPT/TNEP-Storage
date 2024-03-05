import plotly.graph_objects as go
import pandas as pd
import sys

simdir = sys.argv[1]
# Read CSV file
df_storage = pd.read_csv("../../{}/output/storage_investments.csv".format(simdir))
df_lines = pd.read_csv("../../{}/output/line_investments.csv".format(simdir))

fig = go.Figure()

# Add lines
for _, row in df_lines.iterrows():
    fig.add_trace(go.Scattergeo(
        lon = [row['Lon1'], row['Lon2']],
        lat = [row['Lat1'], row['Lat2']],
        mode = 'lines',
        line=dict(width=0.5, color='blue'),
    ))

# Add nodes
fig.add_trace(go.Scattergeo(
    lon = df_storage['Lon'],
    lat = df_storage['Lat'],
    mode = 'markers',
    marker_color = 'red',
    marker_size = 1
    ))

fig.update_geos(
    lonaxis_range=[-106, -94],  # Min and max longitude for the USA
    lataxis_range=[25, 37]      # Min and max latitude for the USA
)

fig.update_layout(showlegend=False)

fig.write_image("../../{}/visual/nodes.png".format(simdir))