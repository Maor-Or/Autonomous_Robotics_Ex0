from datetime import datetime, timezone, timedelta
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
import navpy
from gnssutils import EphemerisManager
import sys, os, csv
import pyproj
from scipy.optimize import least_squares
import simplekml



def export_gnss_data_from_file(input_filepath):
    with open(input_filepath) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if row[0][0] == '#':
                if 'Fix' in row[0]:
                    android_fixes = [row[1:]]
                elif 'Raw' in row[0]:
                    measurements = [row[1:]]
            else:
                if row[0] == 'Fix':
                    android_fixes.append(row[1:])
                elif row[0] == 'Raw':
                    measurements.append(row[1:])

    android_fixes = pd.DataFrame(android_fixes[1:], columns=android_fixes[0])
    measurements = pd.DataFrame(measurements[1:], columns=measurements[0])

    # Format satellite IDs
    measurements.loc[measurements['Svid'].str.len() == 1, 'Svid'] = '0' + measurements['Svid']
    measurements.loc[measurements['ConstellationType'] == '1', 'Constellation'] = 'G'
    measurements.loc[measurements['ConstellationType'] == '3', 'Constellation'] = 'R'
    measurements['SvName'] = measurements['Constellation'] + measurements['Svid']

    # Remove all non-GPS measurements
    measurements = measurements.loc[measurements['Constellation'] == 'G']

    # Convert columns to numeric representation
    measurements['Cn0DbHz'] = pd.to_numeric(measurements['Cn0DbHz'])
    measurements['TimeNanos'] = pd.to_numeric(measurements['TimeNanos'])
    measurements['FullBiasNanos'] = pd.to_numeric(measurements['FullBiasNanos'])
    measurements['ReceivedSvTimeNanos'] = pd.to_numeric(measurements['ReceivedSvTimeNanos'])
    measurements['PseudorangeRateMetersPerSecond'] = pd.to_numeric(measurements['PseudorangeRateMetersPerSecond'])
    measurements['ReceivedSvTimeUncertaintyNanos'] = pd.to_numeric(measurements['ReceivedSvTimeUncertaintyNanos'])

    # A few measurement values are not provided by all phones
    # We'll check for them and initialize them with zeros if missing
    if 'BiasNanos' in measurements.columns:
        measurements['BiasNanos'] = pd.to_numeric(measurements['BiasNanos'])
    else:
        measurements['BiasNanos'] = 0
    if 'TimeOffsetNanos' in measurements.columns:
        measurements['TimeOffsetNanos'] = pd.to_numeric(measurements['TimeOffsetNanos'])
    else:
        measurements['TimeOffsetNanos'] = 0

    # print(measurements.columns)

    measurements['GpsTimeNanos'] = measurements['TimeNanos'] - (
                measurements['FullBiasNanos'] - measurements['BiasNanos'])
    gpsepoch = datetime(1980, 1, 6, 0, 0, 0)
    measurements['UnixTime'] = pd.to_datetime(measurements['GpsTimeNanos'], utc=True, origin=gpsepoch)
    measurements['UnixTime'] = measurements['UnixTime']

    # Split data into measurement epochs
    measurements['Epoch'] = 0
    measurements.loc[
        measurements['UnixTime'] - measurements['UnixTime'].shift() > timedelta(milliseconds=200), 'Epoch'] = 1
    measurements['Epoch'] = measurements['Epoch'].cumsum()

    WEEKSEC = 604800
    LIGHTSPEED = 2.99792458e8

    # This should account for rollovers since it uses a week number specific to each measurement

    measurements['tRxGnssNanos'] = measurements['TimeNanos'] + measurements['TimeOffsetNanos'] - (
                measurements['FullBiasNanos'].iloc[0] + measurements['BiasNanos'].iloc[0])
    measurements['GpsWeekNumber'] = np.floor(1e-9 * measurements['tRxGnssNanos'] / WEEKSEC)
    measurements['tRxSeconds'] = 1e-9 * measurements['tRxGnssNanos'] - WEEKSEC * measurements['GpsWeekNumber']
    measurements['tTxSeconds'] = 1e-9 * (measurements['ReceivedSvTimeNanos'] + measurements['TimeOffsetNanos'])
    # Calculate pseudorange in seconds
    measurements['prSeconds'] = measurements['tRxSeconds'] - measurements['tTxSeconds']

    # Conver to meters
    measurements['PrM'] = LIGHTSPEED * measurements['prSeconds']
    measurements['PrSigmaM'] = LIGHTSPEED * 1e-9 * measurements['ReceivedSvTimeUncertaintyNanos']

    # Initialize with NaN values for updating the df at the end of the parsing
    measurements['Sat.X'] = np.nan
    measurements['Sat.Y'] = np.nan
    measurements['Sat.Z'] = np.nan
    measurements['Pseudo-Range'] = np.nan

    return measurements


def calculate_satellite_position(ephemeris, transmit_time):
    mu = 3.986005e14
    OmegaDot_e = 7.2921151467e-5
    F = -4.442807633e-10
    sv_position = pd.DataFrame()
    sv_position['sv'] = ephemeris.index
    sv_position.set_index('sv', inplace=True)
    sv_position['t_k'] = transmit_time - ephemeris['t_oe']
    A = ephemeris['sqrtA'].pow(2)
    n_0 = np.sqrt(mu / A.pow(3))
    n = n_0 + ephemeris['deltaN']
    M_k = ephemeris['M_0'] + n * sv_position['t_k']
    E_k = M_k
    err = pd.Series(data=[1] * len(sv_position.index))
    i = 0
    while err.abs().min() > 1e-8 and i < 10:
        new_vals = M_k + ephemeris['e'] * np.sin(E_k)
        err = new_vals - E_k
        E_k = new_vals
        i += 1

    sinE_k = np.sin(E_k)
    cosE_k = np.cos(E_k)
    delT_r = F * ephemeris['e'].pow(ephemeris['sqrtA']) * sinE_k
    delT_oc = transmit_time - ephemeris['t_oc']
    sv_position['delT_sv'] = ephemeris['SVclockBias'] + ephemeris['SVclockDrift'] * delT_oc + ephemeris[
        'SVclockDriftRate'] * delT_oc.pow(2)

    v_k = np.arctan2(np.sqrt(1 - ephemeris['e'].pow(2)) * sinE_k, (cosE_k - ephemeris['e']))

    Phi_k = v_k + ephemeris['omega']

    sin2Phi_k = np.sin(2 * Phi_k)
    cos2Phi_k = np.cos(2 * Phi_k)

    du_k = ephemeris['C_us'] * sin2Phi_k + ephemeris['C_uc'] * cos2Phi_k
    dr_k = ephemeris['C_rs'] * sin2Phi_k + ephemeris['C_rc'] * cos2Phi_k
    di_k = ephemeris['C_is'] * sin2Phi_k + ephemeris['C_ic'] * cos2Phi_k

    u_k = Phi_k + du_k

    r_k = A * (1 - ephemeris['e'] * np.cos(E_k)) + dr_k

    i_k = ephemeris['i_0'] + di_k + ephemeris['IDOT'] * sv_position['t_k']

    x_k_prime = r_k * np.cos(u_k)
    y_k_prime = r_k * np.sin(u_k)

    Omega_k = ephemeris['Omega_0'] + (ephemeris['OmegaDot'] - OmegaDot_e) * sv_position['t_k'] - OmegaDot_e * ephemeris[
        't_oe']

    sv_position['x_k'] = x_k_prime * np.cos(Omega_k) - y_k_prime * np.cos(i_k) * np.sin(Omega_k)
    sv_position['y_k'] = x_k_prime * np.sin(Omega_k) + y_k_prime * np.cos(i_k) * np.cos(Omega_k)
    sv_position['z_k'] = y_k_prime * np.sin(i_k)
    return sv_position


def update_df_with_xyz_and_pseudorange(measurements, ephemeris_data_directory):
    manager = EphemerisManager(ephemeris_data_directory)
    LIGHTSPEED = 2.99792458e8

    for epoch in measurements['Epoch'].unique():
        one_epoch = measurements.loc[(measurements['Epoch'] == epoch) & (measurements['prSeconds'] < 0.1)]
        one_epoch = one_epoch.drop_duplicates(subset='SvName').set_index('SvName')
        timestamp = one_epoch.iloc[0]['UnixTime'].to_pydatetime(warn=False)
        sats = one_epoch.index.unique().tolist()
        ephemeris = manager.get_ephemeris(timestamp, sats)
        sv_position = calculate_satellite_position(ephemeris, one_epoch['tTxSeconds'])

        xs = sv_position[['x_k', 'y_k', 'z_k']].to_numpy()
        pr = one_epoch['PrM'] + LIGHTSPEED * sv_position['delT_sv']
        pr = pr.to_numpy()

        for idx, (x, y, z), pr_val in zip(range(len(xs)), xs, pr):
            svid = sats[idx]  # Get the satellite ID
            epoch_mask = (measurements['Epoch'] == epoch) & (measurements['SvName'] == svid)
            measurements.loc[epoch_mask, 'Sat.X'] = x
            measurements.loc[epoch_mask, 'Sat.Y'] = y
            measurements.loc[epoch_mask, 'Sat.Z'] = z
            measurements.loc[epoch_mask, 'Pseudo-Range'] = pr_val

    return measurements


def export_updated_df_to_csv(measurements, parent_directory, csv_output_filename):
    columns_to_keep = ['UnixTime', 'SvName', 'Sat.X', 'Sat.Y', 'Sat.Z', 'Pseudo-Range', 'Cn0DbHz']

    measurements_to_csv = measurements[columns_to_keep].copy()
    measurements_to_csv.rename(columns={'UnixTime': 'GPS time', 'SvName': 'SatPRN (ID)', 'Cn0DbHz': 'CN0'},
                               inplace=True)
    measurements_to_csv.reset_index(drop=True, inplace=True)

    output_filepath = os.path.join(parent_directory, 'data', csv_output_filename)
    measurements_to_csv.to_csv(output_filepath, index=False)
    print(f"{csv_output_filename} has been saved to:", output_filepath)

    return measurements_to_csv


def trilateration(satellite_data):
    # Define the coordinate system transformations
    ecef_to_lla = pyproj.Transformer.from_crs("EPSG:4978", "EPSG:4326", always_xy=True).transform

    # Initial guess for position (using mean of satellite positions)
    x_initial = satellite_data['Sat.X'].mean()
    y_initial = satellite_data['Sat.Y'].mean()
    z_initial = satellite_data['Sat.Z'].mean()

    # Define function for optimization
    def residuals(position):
        x, y, z = position
        predicted_range = np.sqrt((satellite_data['Sat.X'] - x) ** 2 + (satellite_data['Sat.Y'] - y) ** 2 + (
                    satellite_data['Sat.Z'] - z) ** 2)
        difference = predicted_range - satellite_data['Pseudo-Range']
        return difference

    # Perform least squares optimization
    result = least_squares(residuals, [x_initial, y_initial, z_initial])

    # Extract optimized position
    x_final, y_final, z_final = result.x

    # Convert optimized ECEF coordinates to geodetic (lat, lon, alt)
    lon, lat, alt = ecef_to_lla(x_final, y_final, z_final)
    alt = result.fun.mean()  # Use the mean of the residuals as altitude

    return x_final, y_final, z_final, lat, lon, alt


# calculate user's location for each epoch
def calculate_locations_in_df(satellite_data):
    result_coordinates = {}
    # Group the data by 'GPS time' and iterate over each group
    for time, group in satellite_data.groupby('GPS time'):
        # Call the trilateration function for each group
        x_final, y_final, z_final, lat, lon, alt = trilateration(group)
        result_coordinates[time] = (x_final, y_final, z_final, lat, lon, alt)
    return result_coordinates



def update_satellite_df_with_user_location_and_export(satellite_data, result_coordinates, parent_directory, csv_output_filename):
    # Update the original DataFrame with the calculated coordinates
    x_final_list = []
    y_final_list = []
    z_final_list = []
    lat_list = []
    lon_list = []
    alt_list = []

    for index, row in satellite_data.iterrows():
        time = row['GPS time']
        x_final, y_final, z_final, lat, lon, alt = result_coordinates[time]
        x_final_list.append(x_final)
        y_final_list.append(y_final)
        z_final_list.append(z_final)
        lat_list.append(lat)
        lon_list.append(lon)
        alt_list.append(alt)

    # Add new columns for the coordinates to the original DataFrame
    satellite_data['x_final'] = x_final_list
    satellite_data['y_final'] = y_final_list
    satellite_data['z_final'] = z_final_list
    satellite_data['lat'] = lat_list
    satellite_data['lon'] = lon_list
    satellite_data['alt'] = alt_list

    output_filepath = os.path.join(parent_directory, 'data', csv_output_filename)
    satellite_data.to_csv(output_filepath, index=False)
    print(f"{csv_output_filename} has been saved to:", output_filepath)
    return satellite_data


def export_to_kml(satellite_data, parent_directory , kml_output_filename):
    # Assuming `data` is your DataFrame containing satellite data
    # Group the data by 'GPS time' and iterate over each group
    kml = simplekml.Kml()

    for time, group in satellite_data.groupby('GPS time'):
        # Extract coordinates for the group
        lat = group['lat'].iloc[0]
        lon = group['lon'].iloc[0]
        alt = group['alt'].iloc[0]

        # Create a point for the coordinates at this time
        point = kml.newpoint(name=str(time), coords=[(lon, lat, alt)])

    output_filepath = os.path.join(parent_directory, 'data', kml_output_filename)
    # Save the KML file
    kml.save(output_filepath)
    print(f"{kml_output_filename} has been saved to:", output_filepath)


if __name__ == "__main__":
    ### Parsing phase ###
    
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Get the parent directory
    parent_directory = os.path.dirname(script_dir)
    ephemeris_data_directory = os.path.join(parent_directory, 'data')
    sys.path.insert(0, parent_directory)
    # Get path to sample file in data directory, which is located in the parent directory of this notebook
    input_filepath = os.path.join(parent_directory, 'data', 'sample', 'gnss_log_2024_04_13_19_51_17_boaz_fixed.txt')
    df = export_gnss_data_from_file(input_filepath)
    df = update_df_with_xyz_and_pseudorange(df, ephemeris_data_directory)
    satellite_df = export_updated_df_to_csv(df, parent_directory, "ex0_part2.csv")

    ### Location Calculation phase ###
    result_coordinates = calculate_locations_in_df(satellite_df)
    updated_satellite_df = update_satellite_df_with_user_location_and_export(satellite_df, result_coordinates,
                                                                             parent_directory, 'ex0_part5.csv')
    export_to_kml(satellite_df, parent_directory, 'satellite_coordinates.kml')
