using JSON

SUMS_DATA_FILE = "data/topology/tamu/n2/tamu_aggregate_ps_data.json"
MAXES_DATA_FILE = "data/topology/tamu/n2/tamu_aggregate_ps_data.json"

data_sums = JSON.parsefile(SUMS_DATA_FILE)
data_maxes = JSON.parsefile(MAXES_DATA_FILE)