import numpy as np

# Constants
G = 6.67430e-11  # Gravitational constant, m^3 kg^-1 s^-2
M_EARTH = 5.972e24  # Mass of Earth, kg
MU = G * M_EARTH  # Standard gravitational parameter
R_EARTH = 6371000  # Radius of Earth, m

def get_orbital_elements(r_vec, v_vec):
    """
    Calculates orbital elements from state vectors.
    
    Args:
        r_vec (np.array): Position vector (m) relative to central body.
        v_vec (np.array): Velocity vector (m/s).
        
    Returns:
        dict: Orbital elements (a, e, i, omega, w, nu).
    """
    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)
    
    # Specific angular momentum vector h = r x v
    # In 2D (x, y), h is in z direction
    # r_vec is [x, y], v_vec is [vx, vy]
    # h = x*vy - y*vx
    
    # Let's treat vectors as 3D for general calculation
    if len(r_vec) == 2:
        r_vec_3d = np.array([r_vec[0], r_vec[1], 0.0])
        v_vec_3d = np.array([v_vec[0], v_vec[1], 0.0])
    else:
        r_vec_3d = r_vec
        v_vec_3d = v_vec
        
    h_vec = np.cross(r_vec_3d, v_vec_3d)
    h = np.linalg.norm(h_vec)
    
    # Node vector n = k x h (k is unit vector in z)
    # n points towards ascending node
    k = np.array([0, 0, 1])
    n_vec = np.cross(k, h_vec)
    n = np.linalg.norm(n_vec)
    
    # Eccentricity vector e_vec
    # e = (1/mu) * ((v^2 - mu/r)*r - (r.v)*v)
    term1 = (v**2 - MU/r) * r_vec_3d
    term2 = np.dot(r_vec_3d, v_vec_3d) * v_vec_3d
    e_vec = (term1 - term2) / MU
    e = np.linalg.norm(e_vec)
    
    # Specific mechanical energy E = v^2/2 - mu/r
    E = (v**2)/2 - MU/r
    
    # Semi-major axis a = -mu / (2*E)
    if abs(E) < 1e-10:
        a = float('inf') # Parabolic
    else:
        a = -MU / (2*E)
        
    # Inclination i = acos(h_z / h)
    i_rad = np.arccos(h_vec[2] / h)
    i_deg = np.degrees(i_rad)
    
    # Longitude of ascending node (Omega)
    # For 2D equatorial orbits, this is undefined/0
    if n == 0:
        Omega_deg = 0.0
    else:
        Omega_deg = np.degrees(np.arccos(n_vec[0] / n))
        if n_vec[1] < 0:
            Omega_deg = 360 - Omega_deg
            
    # Argument of periapsis (omega)
    if n == 0:
        # Equatorial orbit
        # Angle between I (x-axis) and e_vec
        if e == 0:
            omega_deg = 0.0
        else:
            omega_deg = np.degrees(np.arccos(e_vec[0] / e))
            if e_vec[1] < 0:
                omega_deg = 360 - omega_deg
    else:
        omega_deg = np.degrees(np.arccos(np.dot(n_vec, e_vec) / (n * e)))
        if e_vec[2] < 0:
            omega_deg = 360 - omega_deg

    # True anomaly (nu)
    # Angle between e_vec and r_vec
    if e == 0:
        # Circular orbit
        # Angle between I (x-axis) and r_vec? No, usually argument of latitude u
        # Let's just use angle of r_vec
        nu_deg = 0.0 
    else:
        cos_nu = np.dot(e_vec, r_vec_3d) / (e * r)
        cos_nu = np.clip(cos_nu, -1.0, 1.0)
        nu_deg = np.degrees(np.arccos(cos_nu))
        if np.dot(r_vec_3d, v_vec_3d) < 0:
            nu_deg = 360 - nu_deg

    return {
        'a': a,
        'e': e,
        'i': i_deg,
        'Omega': Omega_deg,
        'omega': omega_deg,
        'nu': nu_deg,
        'period': 2 * np.pi * np.sqrt(a**3 / MU) if a > 0 else 0
    }

def hohmann_transfer(r1, r2):
    """
    Calculates Delta-V for a Hohmann transfer between two circular orbits.
    
    Args:
        r1 (float): Radius of initial orbit (m).
        r2 (float): Radius of target orbit (m).
        
    Returns:
        tuple: (dv1, dv2, total_dv, transfer_time)
    """
    # Velocity of initial orbit
    v1 = np.sqrt(MU / r1)
    
    # Velocity of target orbit
    v2 = np.sqrt(MU / r2)
    
    # Semi-major axis of transfer orbit
    at = (r1 + r2) / 2
    
    # Velocity at periapsis of transfer orbit (at r1)
    vp_t = np.sqrt(MU * (2/r1 - 1/at))
    
    # Velocity at apoapsis of transfer orbit (at r2)
    va_t = np.sqrt(MU * (2/r2 - 1/at))
    
    # Delta-V 1 (at r1)
    dv1 = abs(vp_t - v1)
    
    # Delta-V 2 (at r2)
    dv2 = abs(v2 - va_t)
    
    total_dv = dv1 + dv2
    
    # Transfer time (half period of transfer orbit)
    t_transfer = np.pi * np.sqrt(at**3 / MU)
    
    return dv1, dv2, total_dv, t_transfer

def orbit_state_at_time(a, e, dt):
    """
    Calculates the position in orbit after time dt (simplified for circular/elliptical).
    Does not account for initial true anomaly yet (assumes starts at periapsis).
    """
    n = np.sqrt(MU / a**3) # Mean motion
    M = n * dt # Mean anomaly
    
    # Solve Kepler's Equation M = E - e*sin(E) for Eccentric Anomaly E
    # Newton-Raphson
    E_anom = M
    for _ in range(10):
        f = E_anom - e * np.sin(E_anom) - M
        df = 1 - e * np.cos(E_anom)
        E_anom -= f / df
        
    # Calculate True Anomaly nu
    # tan(nu/2) = sqrt((1+e)/(1-e)) * tan(E/2)
    nu = 2 * np.arctan(np.sqrt((1+e)/(1-e)) * np.tan(E_anom/2))
    
    # Radius r
    r = a * (1 - e * np.cos(E_anom))
    
    return r, nu
