#---------------------------------------------------------------------------------
# IMPORTS
#---------------------------------------------------------------------------------
from scipy.io import loadmat        # For extracting data from matlab .fig/.mat files
import matplotlib.pyplot as plt     # For plotting
import math                         # For basic math functions
import argparse                     # For parsing command line arguments
import ast                          

#---------------------------------------------------------------------------------
# PLOT
#---------------------------------------------------------------------------------
def plot(filename, save=False):
    """
        Grabs a MATLAB .fig file and converts it to a matplotlib format. 
        The save flag allows saving of this data into a png. 

    Args:
        filename (str): 
            Contains the path to the .fig file to be converted and displayed
        save (bool): 
            Flag to determine if the plot should be saved as an image. 

    Returns:
        int: The sum of the two integers.

    Raises:
        ValueError: If either argument is not an integer.
    """
    # Loads the .fig filenames one at a time. 
    d = loadmat(filename,squeeze_me=True, struct_as_record=False)
    
    # Grabs the distinct mat_structs (one per subplot in the .fig, usually)
    # Unfortunately, this requires some guess work to determine what is inside
    # the .fig if you haven't seen it before. So the exact information in the
    # mat_struct can be a little bit annoying to figure out. 
    ax = d['hgS_070000'].children

    # Figure out how many rows and columns in the .fig file. Assume
    # that we are trying to distribute the columns in a square grid. 
    nrows = math.floor(math.sqrt(len(ax)))
    ncols = math.ceil(math.sqrt(len(ax)))
    fig, subplt = plt.subplots(nrows=nrows, ncols=ncols)

    # Use this counter to figure out when we label text and then 
    # take appropriate acctions to populate the title, x-axis and y-axis
    # labels. This can break if the fig is not so simplistic. 
    counter = 0

    # Go through each individual mat_struct and start to figure out 
    # how to graph this.
    for i, a in enumerate(ax):
        # Grab the contents of each mat_struct and assume it is a line
        for line in a.children:
            # Check whether the type of the mat_struct contents is actually a line
            if line.type == 'graph2d.lineseries' or line.type=='line':
                # Check if there is a Marker property in the mat_struct datatype
                # If there is, grab it.
                if hasattr(line.properties,'Marker'):
                    marker = "%s" % line.properties.Marker
                    if (marker == 'diamond'): marker = 'D'
                    elif (marker == 'square'): marker = 's'
                else:
                    marker = None
                    
                # Check if there is a linestyle property in mat_struct datatype
                # If there is, grab it. 
                if hasattr(line.properties,'LineStyle'):
                    linestyle = "%s" % line.properties.LineStyle
                    if (linestyle == 'none'): linestyle=None
                else:
                    linestyle = None
                
                # Check if there is a color property in mat_struct datatype
                # If there is, grab it. 
                if hasattr(line.properties, 'Color'):
                    r,g,b =  line.properties.Color
                    color = (r,g,b)
                else:
                    color = None

                # Create a dictionary of the above parameters
                # We'll pass this dictionary filtered to remove None
                # values as parameters into the plot functions
                params = {
                    'marker': marker,
                    'linestyle': linestyle,
                    'color': color
                }

                # Filter out None values from our parameters
                filtered_params = {k: v for k, v in params.items() if v is not None}

                # Grab the x and y coordinates
                # This may break if the fig doesn't contain this data.
                x = line.properties.XData
                y = line.properties.YData
                
                # Use the filtered parameters to plot. If there is no linestyle
                # then we assume a scatter plot. If there is then we assume a line
                # plot. 
                if (not linestyle and marker):
                    subplt.flat[i].scatter(x, y, **filtered_params)
                else:
                    subplt.flat[i].plot(x,y,**filtered_params)

            # Check if the datatype is actually text. This will be label
            # data most likely. 
            elif line.type == 'text':
                # If the data is a string then we know that it is likely 
                # to be label data. 
                if hasattr(line.properties,'String'):

                    print (line.properties.String)
                    print(len(line.properties.String))
                    # print (counter)
                    
                    # check 
                    if (line.properties.String!=" "):
                        # If its the first label of a set of three encountered, 
                        # then assume that it is the title of the sub plot
                        if (counter == 0):
                            subplt.flat[i].set_title(line.properties.String)
                            counter = 1
                        # If its the second label of a set of three, then 
                        # assume it is an x-axis label
                        elif (counter == 1):
                            subplt.flat[i].set_xlabel(line.properties.String)
                            counter = 2
                        # If its the third label of a set of three, then
                        # assume it is an y-axis label 
                        elif (counter == 2):
                            subplt.flat[i].set_ylabel(line.properties.String)
                            counter = 0
                
    # Adjust layout to be tight
    plt.tight_layout()

    # Show the figure. The figure will need to be closed to terminate
    # the code. 
    plt.show()

# Save the plot as a PNG file with 600 DPI
# plt.savefig(filename+'.png', format='png', dpi=1000)


if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="This script allows display of a MATLAB .fig file without having MATLAB. It might be a bit buggy since there may be many different structures to a .fig file and this assumes a grid of line subplots. The function would need to be customized for other kinds of plots. But this can be a guessing game unless the user has an idea of what the original MATLAB plot looked like.")

    # Add arguments
    parser.add_argument('-f', '--filename', type=str, help='the path to the .fig file', required=True)

    # Parse the arguments
    args = parser.parse_args()

    # Grab the filename
    filename=args.filename

    # Display the plot
    plot(filename)
