import plotly.express as px
import pandas as pd
import sys

simdir = sys.argv[1]
# Read CSV file
df_energy = pd.read_csv("../../{}/output/energy.csv".format(simdir))
df_energy['Imbalance_Sign'] = df_energy['Energy_Imbalance'].apply(lambda x: 'Positive' if x >= 0 else 'Negative')
df_energy['Imbalance_Size'] = df_energy['Energy_Imbalance'].abs()

fig = px.scatter_geo(df_energy, 
                     lon='Lon',
                     lat='Lat',
                     hover_name='Node_Name', 
                     size='Imbalance_Size',
                     color='Imbalance_Sign',
                     color_discrete_map={'Positive': 'blue', 'Negative': 'red'},
                     animation_frame='Hour',
                     scope='usa')

fig.write_html("../../{}/visual/energy.html".format(simdir))