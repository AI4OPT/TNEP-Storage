import plotly.express as px
import pandas as pd
import sys

simdir = sys.argv[1]
# Read CSV file
df_energy = pd.read_csv(f"../../{simdir}/output/energy.csv")

# Create a new column for plotting (size)
df_energy['Imbalance_Size'] = df_energy['Energy_Imbalance'].apply(lambda x: abs(x) if x < 0 else 0)

# Generate the plot focusing only on negatives, shown in red
fig = px.scatter_geo(df_energy, 
                     lon='Lon',
                     lat='Lat',
                     hover_name='Node_Name', 
                     size='Imbalance_Size',
                     color_discrete_sequence=['red'],  # Use a single color for clarity
                     animation_frame='Hour',
                     scope='usa')

# Save the plot as an HTML file
fig.write_html(f"../../{simdir}/visual/ue_energy.html")
