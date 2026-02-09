import numpy as np
import matplotlib.pyplot as plt
import os
from rocket import Rocket, RocketStage, R_EARTH
import orbit

def run_simulation():
    # Define Rocket Stages (Approximate Falcon 9)
    # Stage 1
    stage1 = RocketStage(
        dry_mass=25600, 
        fuel_mass=395700, 
        thrust=7607000, 
        isp=300, 
        diameter=3.7, 
        drag_coeff=0.3
    )
    
    # Stage 2
    stage2 = RocketStage(
        dry_mass=3900, 
        fuel_mass=92670, 
        thrust=934000, 
        isp=348, 
        diameter=3.7, 
        drag_coeff=0.3
    )
    
    payload_mass = 10000  # Payload to LEO
    
    rocket = Rocket("Falcon 9-ish", [stage1, stage2], payload_mass)
    
    dt = 0.1
    max_time = 3000  # seconds
    
    print("Starting simulation...")
    
    # Simulation Loop
    while rocket.time < max_time:
        alt = rocket.step(dt)
        
        # Stop if we crash
        if alt < -100:
            print("Rocket crashed!")
            break
            
        # Stop if we are in stable orbit (simplified check)
        # v_orbit = sqrt(GM/r)
        r = np.linalg.norm(rocket.position)
        v = np.linalg.norm(rocket.velocity)
        
        # If we are coasting and have been flying for a while, stop
        if rocket.current_stage_index >= len(rocket.stages) and rocket.time > 900:
            print("Mission complete (fuel exhausted or target reached)")
            break

    # Calculate Final Orbit
    final_r = rocket.position
    final_v = rocket.velocity
    elements = orbit.get_orbital_elements(final_r, final_v)
    
    print("\nFinal Orbital Elements:")
    for key, val in elements.items():
        print(f"{key}: {val:.4f}")
        
    # Hohmann Transfer Calculation
    # Assume we are in a circular orbit at the final altitude (approx)
    r_parking = np.linalg.norm(final_r)
    alt_parking = r_parking - R_EARTH
    r_geo = R_EARTH + 35786000 # GEO altitude
    
    print(f"\nCalculating Hohmann Transfer from {alt_parking/1000:.2f} km to GEO ({35786} km)...")
    dv1, dv2, total_dv, t_transfer = orbit.hohmann_transfer(r_parking, r_geo)
    
    print(f"Delta-V 1 (LEO burn): {dv1:.2f} m/s")
    print(f"Delta-V 2 (GEO circularization): {dv2:.2f} m/s")
    print(f"Total Delta-V: {total_dv:.2f} m/s")
    print(f"Transfer Time: {t_transfer/3600:.2f} hours")
    
    # Plotting
    plot_results(rocket, elements, r_parking, r_geo)

def plot_results(rocket, elements, r1, r2):
    # Ensure images directory exists
    if not os.path.exists("images"):
        os.makedirs("images")
        
    history = rocket.history
    t = np.array(history['time'])
    alt = np.array(history['altitude']) / 1000 # km
    vel = np.array(history['velocity'])
    q = np.array(history['dynamic_pressure']) / 1000 # kPa
    acc = np.array(history['acceleration']) / 9.81 # g
    
    # Figure 1: Ascent Profile
    fig1, axs = plt.subplots(2, 2, figsize=(12, 10))
    fig1.suptitle(f"Ascent Profile: {rocket.name}")
    
    axs[0, 0].plot(t, alt)
    axs[0, 0].set_title("Altitude vs Time")
    axs[0, 0].set_xlabel("Time (s)")
    axs[0, 0].set_ylabel("Altitude (km)")
    axs[0, 0].grid(True)
    
    axs[0, 1].plot(t, vel)
    axs[0, 1].set_title("Velocity vs Time")
    axs[0, 1].set_xlabel("Time (s)")
    axs[0, 1].set_ylabel("Velocity (m/s)")
    axs[0, 1].grid(True)
    
    axs[1, 0].plot(t, q)
    axs[1, 0].set_title("Dynamic Pressure (Q) vs Time")
    axs[1, 0].set_xlabel("Time (s)")
    axs[1, 0].set_ylabel("Q (kPa)")
    axs[1, 0].grid(True)
    
    axs[1, 1].plot(t, acc)
    axs[1, 1].set_title("Acceleration vs Time")
    axs[1, 1].set_xlabel("Time (s)")
    axs[1, 1].set_ylabel("Acceleration (g)")
    axs[1, 1].grid(True)
    
    plt.tight_layout()
    plt.savefig("images/ascent_profile.png")
    print("Saved images/ascent_profile.png")
    
    # Figure 2: Trajectory
    plt.figure(figsize=(10, 6))
    downrange = np.array(history['downrange_distance']) / 1000 # km
    plt.plot(downrange, alt)
    plt.title("Ascent Trajectory")
    plt.xlabel("Downrange Distance (km)")
    plt.ylabel("Altitude (km)")
    plt.grid(True)
    plt.savefig("images/trajectory.png")
    print("Saved images/trajectory.png")
    
    # Figure 3: Orbit Visualization
    # Draw Earth
    theta = np.linspace(0, 2*np.pi, 100)
    x_earth = R_EARTH * np.cos(theta)
    y_earth = R_EARTH * np.sin(theta)
    
    plt.figure(figsize=(10, 10))
    plt.plot(x_earth/1000, y_earth/1000, 'b-', label='Earth')
    
    # Parking Orbit (simplified as circular)
    x_parking = r1 * np.cos(theta)
    y_parking = r1 * np.sin(theta)
    plt.plot(x_parking/1000, y_parking/1000, 'g--', label='Parking Orbit')
    
    # Target Orbit (GEO)
    x_geo = r2 * np.cos(theta)
    y_geo = r2 * np.sin(theta)
    plt.plot(x_geo/1000, y_geo/1000, 'r--', label='GEO')
    
    # Transfer Orbit
    # Ellipse with a = (r1+r2)/2, e = (r2-r1)/(r2+r1)
    a_t = (r1 + r2) / 2
    e_t = (r2 - r1) / (r2 + r1)
    # Focus is at (0,0). Periapsis at r1 (right side for visualization).
    # r = a(1-e^2)/(1+e cos theta)
    r_transfer = a_t * (1 - e_t**2) / (1 + e_t * np.cos(theta))
    x_transfer = r_transfer * np.cos(theta)
    y_transfer = r_transfer * np.sin(theta)
    
    # We only need the half ellipse from periapsis to apoapsis
    # But full ellipse is fine to show
    plt.plot(x_transfer/1000, y_transfer/1000, 'y-', label='Hohmann Transfer')
    
    plt.title("PySpaceFlight: Orbital Maneuvers Visualization")
    plt.xlabel("X (km)")
    plt.ylabel("Y (km)")
    plt.axis('equal')
    plt.legend()
    plt.grid(True)
    plt.savefig("images/orbit_viz.png")
    print("Saved images/orbit_viz.png")

if __name__ == "__main__":
    run_simulation()
