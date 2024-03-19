import plotly.graph_objects as go
import pandas as pd
import sys

simdir = sys.argv[1]
# Read CSV file
df_storage = pd.read_csv("../../{}/output/storage_investments.csv".format(simdir))
df_lines = pd.read_csv("../../{}/output/line_investments.csv".format(simdir))


fig = go.Figure(data=go.Scattergeo(
        lon = df_storage['Lon'],
        lat = df_storage['Lat'],
        mode = 'markers',
        marker = dict(
            color = 'red',
            size = df_storage['Storage_Power'],  # Adjust marker size based on 'Storage_Power'
            sizemode = 'area',  # This ensures the size of the marker represents the area
            sizeref = 2.*max(df_storage['Storage_Power'])/(25.**2),  # Adjusts the scale of the marker sizes
            sizemin = 4  # Minimum marker size
            )
        ))

# Add lines
for _, row in df_lines.iterrows():
    fig.add_trace(go.Scattergeo(
        lon = [row['Lon1'], row['Lon2']],
        lat = [row['Lat1'], row['Lat2']],
        mode = 'lines',
        line=dict(width=row['Upgrade_Lvl'], color='blue'),
    ))

fig.update_layout(
        title = 'Nodes',
        geo_scope='usa',
    )
fig.write_image("../../{}/visual/upgrades.png".format(simdir))