import plotly.graph_objects as go
import pandas as pd
import plotly.express as px
import numpy as np
import json
import sys

FILE_PATH = "/storage/home/hcoda1/1/kwu381/TNEP-Storage/data/topology/tamu/n2/tamu_aggregate_ps_data.json"

# Open the JSON file and load its contents into a Python object
with open(FILE_PATH, 'r') as file:
    data = json.load(file)

fig = go.Figure()

for i in data["branch"]:
    branch = data["branch"][i]
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    lon1 = data["bus"][str(f_bus)]["lon"]
    lon2 = data["bus"][str(t_bus)]["lon"]
    lat1 = data["bus"][str(f_bus)]["lat"]
    lat2 = data["bus"][str(t_bus)]["lat"]

    fig.add_trace(go.Scattergeo(
        lon = [lon1, lon2],
        lat = [lat1, lat2],
        mode = 'lines',
        line=dict(width=1, color='red'),
    ))

for i in data["bus"]:
    bus = data["bus"][i]
    fig.add_trace(go.Scattergeo(
        lon = [bus["lon"]],
        lat = [bus["lat"]], 
        mode = 'markers',  
        marker=dict(size=4, color='blue'),  
    ))


fig.update_layout(
    title='Nodes',
    geo_scope='usa',
    showlegend=False
)

fig.write_image("../data/topology/tamu/n2/topology.png")

    