from geopy.geocoders import Nominatim
from skyfield.api import load, Topos, Star
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from skyfield.data import hipparcos
import plotly.graph_objects as go

# Function to get the coordinates for a given city
def get_coordinates(city):
    geolocator = Nominatim(user_agent="geo_app")
    location = geolocator.geocode(city)
    if location:
        return location.latitude, location.longitude
    else:
        raise ValueError(f"Could not find coordinates for {city}.")

# Function to get visible stars
def get_visible_stars(latitude, longitude, date_time):
    ts = load.timescale()
    eph = load('de421.bsp')
    observer = eph['earth'] + Topos(latitude_degrees=latitude, longitude_degrees=longitude)
    observation_time = ts.utc(date_time.year, date_time.month, date_time.day, date_time.hour, date_time.minute)

    # Load star catalog
    with load.open(hipparcos.URL) as f:
        stars = hipparcos.load_dataframe(f)

    hip_id = 677
    sirius_data = stars.loc[hip_id]
    sirius_star = Star.from_dataframe(sirius_data)

    # Compute visibility
    visible_stars = []
    for hip_id, star_data in stars.iterrows():
        star = Star.from_dataframe(star_data)
        apparent = observer.at(observation_time).observe(star).apparent()
        alt, az, _ = apparent.altaz()
        if alt.degrees > 0:
            visible_stars.append((hip_id, alt.degrees, az.degrees))
    return visible_stars

# Function to plot star map
def plot_star_map(visible_stars):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, polar=True)

    for hip_id, alt, az in visible_stars:
        ax.scatter(np.radians(az), 90 - alt, s=10, label=f"HIP {hip_id}")

    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 90)
    ax.set_title("Visible Stars", va='bottom')
    plt.legend(loc='upper right')
    plt.show()

# Function to plot star map
def plot_star_map_3d(visible_stars):
    # Convert altitude and azimuth to Cartesian coordinates for 3D plotting
    x_coords = []
    y_coords = []
    z_coords = []
    star_names = []

    for hip_id, alt, az in visible_stars:
        # Convert spherical (alt, az) to Cartesian
        r = 1  # Assume all stars are at a fixed distance for plotting
        x = r * np.cos(np.radians(alt)) * np.cos(np.radians(az))
        y = r * np.cos(np.radians(alt)) * np.sin(np.radians(az))
        z = r * np.sin(np.radians(alt))

        x_coords.append(x)
        y_coords.append(y)
        z_coords.append(z)
        star_names.append(f"HIP {hip_id}")

    # Create 3D scatter plot
    fig = go.Figure()

    fig.add_trace(go.Scatter3d(
        x=x_coords,
        y=y_coords,
        z=z_coords,
        mode='markers',
        marker=dict(
            size=5,
            color=z_coords,  # Use z-coordinate for color mapping
            colorscale='Viridis',  # Colorscale for depth
            opacity=0.8
        ),
        text=star_names,  # Display HIP ID on hover
        hoverinfo='text'
    ))

    # Add labels and layout adjustments
    fig.update_layout(
        title="3D Star Map",
        scene=dict(
            xaxis_title="X (Azimuth)",
            yaxis_title="Y (Altitude)",
            zaxis_title="Z (Elevation)",
            aspectmode="cube"  # Keep proportions equal
        )
    )

    # Show the plot
    fig.show()


if __name__ == "__main__":
    try:
        # Ask user for input
        city = input("Enter the city you are viewing from: ")
        date_input = input("Enter the date (YYYY-MM-DD): ")
        time_input = input("Enter the time (HH:MM in 24-hour format): ")

        # Validate and parse the date and time
        try:
            date_time = datetime.strptime(f"{date_input} {time_input}", "%Y-%m-%d %H:%M")
        except ValueError:
            raise ValueError("Invalid date or time format. please use  YYYY-MM-DD for date and HH:MM for time.")
        
        # Get coordinates for the city
        latitude, longitude = get_coordinates(city)
        print(f"Coordinates for {city}: Latitude {latitude}, Longitude {longitude}")

        # Get visible stars
        visible_stars = get_visible_stars(latitude, longitude, date_time)
        print(f"Found {len(visible_stars)} visible stars.")

        # Plot star map
        # plot_star_map(visible_stars)

        # Plot star map
        plot_star_map_3d(visible_stars)

    except Exception as e:
        print(f"Error: {e}")

# This is a test comment.