import tkinter as tk
from tkinter import ttk, filedialog
import subprocess

import datetime

# Function to get a formatted timestamp
def get_timestamp():
    return datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

log_filename = f"m2d_{get_timestamp()}.log"

# Function to run the m2d init command
def run_m2d_init_command():
    command = ["m2d", "init"]
    # Add arguments to the command based on the user's input for m2d init
    if config_var.get():
        command.extend(["-c", config_var.get()])
    if approach_var.get():
        command.extend(["-a", approach_var.get()])
    if maxs_var.get():
        command.extend(["-m", maxs_var.get()])
    if number_var.get():
        command.extend(["-n", number_var.get()])
    if direction_var.get():
        command.extend(["-d"] + direction_var.get().split())
    if ranges_var.get():
        command.extend(["-r"] + ranges_var.get().split())
    if properties_var.get():
        command.extend(["-p", properties_var.get()])
    if verbose_var.get():
        command.append("-v")
    if back_var.get():
        command.append("-b")
    
    # Run the command
    with open(log_filename, "a") as log_file:
        subprocess.run(command, stdout=log_file, stderr=subprocess.STDOUT)


# Function to run the m2d run command
def run_m2d_run_command():
    command = ["m2d", "run"]
    # Add arguments to the command based on the user's input for m2d run
    if approach_run_var.get():
        command.extend(["-a", approach_run_var.get()])
    if properties_run_var.get():
        command.extend(["-p", properties_run_var.get()])
    if manual_var.get():
        command.append("--manual")
    if verbose_run_var.get():
        command.append("-v")
    if input_var.get():
        command.append(input_var.get())
    
    # Run the command
    with open(log_filename, "a") as log_file:
        subprocess.run(command, stdout=log_file, stderr=subprocess.STDOUT)

# Function to run the m2d post command
def run_m2d_post_command():
    command = ["m2d", "post"]
    # Add arguments to the command based on the user's input for m2d post
    if approach_post_var.get():
        command.extend(["-a", approach_post_var.get()])
    if inputfile_var.get():
        command.extend(["-i", inputfile_var.get()])
    if properties_post_var.get():
        command.extend(["-p", properties_post_var.get()])
    if skip_var.get():
        command.append("--skip")
    if order_var.get():
        command.extend(["-o", order_var.get()])
    if fmt_var.get():
        command.extend(["-f", fmt_var.get()])
    if dpi_var.get():
        command.extend(["-d", dpi_var.get()])
    if plot_var.get():
        command.append("--plot")
    if verbose_post_var.get():
        command.append("-v")
    
    # Run the command
    with open(log_filename, "a") as log_file:
        subprocess.run(command, stdout=log_file, stderr=subprocess.STDOUT)

# Function to open a file dialog and update the entry with the selected file path
def browse_file(entry_widget):
    filename = filedialog.askopenfilename()
    entry_widget.delete(0, tk.END)
    entry_widget.insert(0, filename)

# Function to update default values based on approach
def update_approach_init_var(*args):
    properties=properties_var.get()
    approach = approach_var.get()
    if properties=="elc":
       direction_var.set(value="")       
       ranges_var.set(value="")
       if approach == 'energy':
           maxs_var.set("0.030")  # Update to default max strain for energy approach
           # Update other default values specific to the energy approach if needed
       elif approach == 'stress':
           maxs_var.set("0.0050")  # Update to default max strain for stress approach
           # Update other default values specific to the stress approach if needed
    else:
       direction_var.set(value="xx")       
       ranges_var.set(value="0.0 0.2")


def update_approach_post_var(*args):
    approach = approach_post_var.get()
    if approach == 'energy':
        order_var.set(4)
    elif approach == 'stress':
        order_var.set(3)

# Create the main window
root = tk.Tk()
root.title("Mech2D")

# Create the tab control
tab_control = ttk.Notebook(root)

# Create the tab for m2d init
tab_init = ttk.Frame(tab_control)
tab_control.add(tab_init, text='m2d init')

# Create the tab for m2d run
tab_run = ttk.Frame(tab_control)
tab_control.add(tab_run, text='m2d run')

# Create the tab for m2d post
tab_post = ttk.Frame(tab_control)
tab_control.add(tab_post, text='m2d post')

# Define variables for each parameter for m2d init
config_var = tk.StringVar(value="")  # Assuming no default file
approach_var = tk.StringVar(value="stress")  # Default approach
maxs_var = tk.StringVar(value="0.0050")  # Default max strain for stress approach
number_var = tk.StringVar(value="5")  # Default odd number > 4
direction_var = tk.StringVar(value="")  # Default direction
ranges_var = tk.StringVar(value="")  # Default strain range
properties_var = tk.StringVar(value="elc")  # Default property to calculate
verbose_var = tk.BooleanVar(value=False)  # Default verbose is off
back_var = tk.BooleanVar(value=False)  # Default back is off

# Bind the update_defaults function to the approach_var
properties_var.trace('w',update_approach_init_var)
approach_var.trace('w', update_approach_init_var)

# Define variables for each parameter for m2d run
approach_run_var = tk.StringVar(value="stress")
properties_run_var = tk.StringVar(value="elc")
manual_var = tk.BooleanVar(value=False)
verbose_run_var = tk.BooleanVar(value=False)
input_var = tk.StringVar(value="")

# Define variables for each parameter for m2d post
approach_post_var = tk.StringVar(value="stress")
inputfile_var = tk.StringVar(value="")
properties_post_var = tk.StringVar(value="elc")
skip_var = tk.BooleanVar(value=False)
order_var = tk.StringVar(value=3)
fmt_var = tk.StringVar(value="png")
dpi_var = tk.StringVar(value=100)
plot_var = tk.BooleanVar(value=True)
verbose_post_var = tk.BooleanVar(value=False)

approach_post_var.trace('w', update_approach_post_var)
#-----------------------------------------
# GUI elements for m2d init inside tab_init
# Labels, entries, checkboxes, and buttons for each parameter
tk.Label(tab_init, text="Config File:").grid(row=0, column=0, sticky='w')
config_entry = tk.Entry(tab_init, textvariable=config_var)
config_entry.grid(row=0, column=1)
tk.Button(tab_init, text="Browse", command=lambda: browse_file(config_entry)).grid(row=0, column=2)

tk.Label(tab_init, text="Approach:").grid(row=1, column=0, sticky='w')
tk.OptionMenu(tab_init, approach_var, "stress", "energy").grid(row=1, column=1, sticky='w')

tk.Label(tab_init, text="Properties:").grid(row=2, column=0, sticky='w')
tk.OptionMenu(tab_init, properties_var, "elc", "ssc").grid(row=2, column=1, sticky='w')

tk.Label(tab_init, text="Max Strain:").grid(row=3, column=0, sticky='w')
maxs_entry = tk.Entry(tab_init, textvariable=maxs_var)
maxs_entry.grid(row=3, column=1)

tk.Label(tab_init, text="Number of Structures:").grid(row=4, column=0, sticky='w')
number_entry = tk.Entry(tab_init, textvariable=number_var)
number_entry.grid(row=4, column=1)

tk.Label(tab_init, text="Direction:").grid(row=5, column=0, sticky='w')
direction_entry = tk.Entry(tab_init, textvariable=direction_var)
direction_entry.grid(row=5, column=1)

tk.Label(tab_init, text="Strain Ranges:").grid(row=6, column=0, sticky='w')
ranges_entry = tk.Entry(tab_init, textvariable=ranges_var)
ranges_entry.grid(row=6, column=1)

tk.Checkbutton(tab_init, text="Verbose", variable=verbose_var).grid(row=7, column=0, sticky='w')
tk.Checkbutton(tab_init, text="Backup", variable=back_var).grid(row=7, column=1, sticky='w')

tk.Button(tab_init, text="Run m2d init", command=run_m2d_init_command).grid(row=8, column=0, columnspan=3)


# GUI elements for m2d run inside tab_run
# Labels, entries, checkboxes, and buttons for each parameter
tk.Label(tab_run, text="Input File:").grid(row=0, column=0, sticky='w')
input_entry = tk.Entry(tab_run, textvariable=input_var)
input_entry.grid(row=0, column=1)
tk.Button(tab_run, text="Browse", command=lambda: browse_file(input_entry)).grid(row=0, column=2)

tk.Label(tab_run, text="Approach:").grid(row=1, column=0, sticky='w')
tk.OptionMenu(tab_run, approach_run_var, "stress", "energy").grid(row=1, column=1, sticky='w')

tk.Label(tab_run, text="Properties:").grid(row=2, column=0, sticky='w')
tk.OptionMenu(tab_run, properties_run_var, "elc", "ssc").grid(row=2, column=1, sticky='w')

tk.Checkbutton(tab_run, text="Manual", variable=manual_var).grid(row=3, column=0, sticky='w')
tk.Checkbutton(tab_run, text="Verbose", variable=verbose_run_var).grid(row=3, column=1, sticky='w')

tk.Button(tab_run, text="Run m2d run", command=run_m2d_run_command).grid(row=4, column=0, columnspan=3)


# GUI elements for m2d post inside tab_post
tk.Label(tab_post, text="Approach:").grid(row=1, column=0, sticky='w')
tk.OptionMenu(tab_post, approach_post_var, "stress", "energy").grid(row=1, column=1, sticky='w')

tk.Label(tab_post, text="Input File:").grid(row=0, column=0, sticky='w')
inputfile_entry = tk.Entry(tab_post, textvariable=inputfile_var)
inputfile_entry.grid(row=0, column=1)
tk.Button(tab_post, text="Browse", command=lambda: browse_file(inputfile_entry)).grid(row=0, column=2)

tk.Label(tab_post, text="Properties:").grid(row=2, column=0, sticky='w')
tk.OptionMenu(tab_post, properties_post_var, "elc", "ssc").grid(row=2, column=1, sticky='w')

tk.Checkbutton(tab_post, text="Skip Data Parsing", variable=skip_var).grid(row=3, column=0, sticky='w')

tk.Label(tab_post, text="Polynomial Order:").grid(row=4, column=0, sticky='w')
order_entry = tk.Entry(tab_post, textvariable=order_var)
order_entry.grid(row=4, column=1)

tk.Label(tab_post, text="Output Format:").grid(row=5, column=0, sticky='w')
fmt_entry = tk.Entry(tab_post, textvariable=fmt_var)
fmt_entry.grid(row=5, column=1)

tk.Label(tab_post, text="Output DPI:").grid(row=6, column=0, sticky='w')
dpi_entry = tk.Entry(tab_post, textvariable=dpi_var)
dpi_entry.grid(row=6, column=1)

tk.Checkbutton(tab_post, text="Plot Figures", variable=plot_var).grid(row=7, column=0, sticky='w')

tk.Checkbutton(tab_post, text="Verbose", variable=verbose_post_var).grid(row=7, column=1, sticky='w')

tk.Button(tab_post, text="Run m2d post", command=run_m2d_post_command).grid(row=8, column=0, columnspan=3)

# Pack the tab control
tab_control.pack(expand=1, fill='both')

# Run the main loop
root.mainloop()

