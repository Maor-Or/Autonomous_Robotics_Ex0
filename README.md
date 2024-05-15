# Autonomous Robotics: Ex0 - GNSS Raw Measurements

Submitted by Raz Saad, Maor Or

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

## Tests
Our tests can be found at data/boaz_data_results directory, and at data/our_data_results. <br/>
An example of a kml output file from data/boaz_data_results/satellite_coordinates_boaz_driving.kml in Google Earth: <br/><br/>
![boaz_driving_output_kml_googleearth](https://github.com/Raz-Saad/Autonomous_Robotics_Ex0/assets/118377261/b3776ea6-a3e5-4885-b5d5-b92970450e51)


## About the algorithm

### Least Squares Trilateration Algorithm

This repository contains an implementation of the least squares trilateration algorithm in Python. Trilateration is a technique used in navigation systems to determine the position of a receiver based on the distances to multiple satellites with known positions.

#### Description

The algorithm takes as input:

- `satellite_positions`: Positions of satellites in Earth-Centered, Earth-Fixed (ECEF) coordinate system.
- `measured_pseudorange`: Measured pseudorange (distance from the receiver to each satellite).
- `initial_position_guess`: Initial guess of the receiver's position.
- `initial_clock_bias`: Initial guess of the clock bias.

It iteratively refines the initial guess of the receiver's position and clock bias using least squares optimization until a satisfactory solution is found.

#### Algorithm Steps

1. **Initialization**: Initialize position and clock bias correction vectors, and set up the G matrix.
2. **Iteration**: Enter a loop until position correction becomes sufficiently small.
   - Calculate range estimates.
   - Predict pseudoranges.
   - Calculate residuals.
   - Update G matrix.
   - Calculate corrections using least squares estimation.
   - Update position and clock bias.
3. **Conversion**: Convert optimized ECEF coordinates to geodetic coordinates (latitude, longitude, altitude).
4. **Output**: Return optimized ECEF and geodetic coordinates of the receiver.
