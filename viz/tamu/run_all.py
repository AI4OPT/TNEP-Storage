import subprocess
import sys

def run_script(script_name, simdir):
    subprocess.run(["python3", script_name, simdir], check=True)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 run_all.py simdir")
        sys.exit(1)

    simdir = sys.argv[1]

    # List of your scripts
    scripts = ["plot_nodes.py", "plot_ue_oe.py", "plot_hourly_generation.py", "plot_investments.py", "plot_curtailment.py", "plot_congested_lines.py"]

    for script in scripts:
        run_script(script, simdir)