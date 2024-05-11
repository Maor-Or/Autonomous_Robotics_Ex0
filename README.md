# Autonomous Robotics: Ex0 - GNSS Raw Measurements

This project involves implementing a basic positioning algorithm for Global Navigation Satellite Systems (GNSS), focusing on Root Mean Square (RMS) of weighted pseudo-ranges. The objective is to convert GNSS raw measurement files into CSV format, compute X, Y, Z coordinates based on GPS time, and convert them to Latitude, Longitude, and Altitude. The solution integrates these components to generate a KML file depicting the computed path and a CSV file with comprehensive positional information. Testing is conducted on provided and custom datasets, with the final solution hosted on GitHub.

## How to run:

Clone this repository to your local machine.
Prepare your GNSS raw measurement files or use provided datasets.<br/>

1: In the terminal, while in the project directory, run the command "pip install -r requirements.txt" to install all requirements.<br/><br/>
2: Optional - To choose your own GNSS txt file to convert to kml, load the file into the data/sample/ directory, then in Running_files/Main.py, in the main function, update the input_filepath variable with your own file's name<br/>  
3: In the terminal, in the Running_files directory, run the command "python Main.py", it will generate three output files (found in the data directory):<br/>
&nbsp;a - ex0_part2.csv - the required output for part 2 of the assignment <br/>
&nbsp;b - ex0_part5.csv - the required output for part 5 of the assignment (like part 2, but includes user's X,Y,Z, Lat,Lon,Alt) <br/>
&nbsp;c - satellite_coordinates.kml
