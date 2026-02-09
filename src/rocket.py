import numpy as np

# Constants
G = 6.67430e-11  # Gravitational constant, m^3 kg^-1 s^-2
M_EARTH = 5.972e24  # Mass of Earth, kg
R_EARTH = 6371000  # Radius of Earth, m
RHO_0 = 1.225  # Sea level air density, kg/m^3
H_SCALE = 8500  # Scale height for atmosphere, m
G0 = 9.80665  # Standard gravity, m/s^2


class RocketStage:
    def __init__(self, dry_mass, fuel_mass, thrust, isp, diameter, drag_coeff=0.5):
        """
        Initializes a rocket stage.

        Args:
            dry_mass (float): Mass of the stage structure without fuel (kg).
            fuel_mass (float): Mass of the fuel (kg).
            thrust (float): Thrust force (N).
            isp (float): Specific impulse (s).
            diameter (float): Diameter of the stage (m).
            drag_coeff (float): Drag coefficient (dimensionless).
        """
        self.dry_mass = dry_mass
        self.fuel_mass = fuel_mass
        self.thrust = thrust
        self.isp = isp
        self.diameter = diameter
        self.area = np.pi * (diameter / 2)**2
        self.drag_coeff = drag_coeff
        self.total_mass = dry_mass + fuel_mass

    def burn(self, dt):
        """
        Simulates fuel burn for a time step dt.
        Returns the thrust force and mass flow rate.
        """
        if self.fuel_mass <= 0:
            return 0.0, 0.0

        mdot = self.thrust / (self.isp * G0)
        fuel_burned = mdot * dt

        if fuel_burned > self.fuel_mass:
            fuel_burned = self.fuel_mass
            mdot = fuel_burned / dt # adjust mdot for last step
            self.fuel_mass = 0
        else:
            self.fuel_mass -= fuel_burned
        
        self.total_mass = self.dry_mass + self.fuel_mass
        return self.thrust, mdot


class Rocket:
    def __init__(self, name, stages, payload_mass):
        """
        Initializes a multi-stage rocket.

        Args:
            name (str): Name of the rocket.
            stages (list): List of RocketStage objects, ordered from bottom to top.
            payload_mass (float): Mass of the payload (kg).
        """
        self.name = name
        self.stages = stages
        self.payload_mass = payload_mass
        self.current_stage_index = 0
        self.time = 0.0
        self.position = np.array([0.0, R_EARTH])  # x, y (Earth centered, y is up at launch)
        self.velocity = np.array([0.0, 0.0])      # vx, vy
        self.history = {
            'time': [],
            'altitude': [],
            'velocity': [],
            'acceleration': [],
            'mass': [],
            'dynamic_pressure': [],
            'downrange_distance': []
        }

    @property
    def current_mass(self):
        mass = self.payload_mass
        for i in range(self.current_stage_index, len(self.stages)):
            mass += self.stages[i].total_mass
        return mass

    def get_atmospheric_density(self, altitude):
        if altitude < 0:
            return RHO_0
        return RHO_0 * np.exp(-altitude / H_SCALE)

    def get_gravity(self, position):
        r = np.linalg.norm(position)
        if r == 0:
            return np.array([0.0, -G0]) # Should not happen
        
        # Gravity vector points towards Earth center (0,0)
        direction = -position / r
        magnitude = G * M_EARTH / (r**2)
        return direction * magnitude

    def get_drag(self, velocity, density):
        v_mag = np.linalg.norm(velocity)
        if v_mag == 0:
            return np.array([0.0, 0.0])
        
        # If coasting (no active stage), use the last stage's geometry or payload geometry
        # simplified: use last stage geometry if available
        idx = min(self.current_stage_index, len(self.stages) - 1)
        stage = self.stages[idx]
        
        drag_force_mag = 0.5 * density * (v_mag**2) * stage.drag_coeff * stage.area
        drag_direction = -velocity / v_mag
        return drag_direction * drag_force_mag

    def step(self, dt):
        """
        Advances the simulation by time step dt.
        """
        if self.current_stage_index >= len(self.stages):
            # Coasting phase or mission complete
            thrust_force = 0.0
            mdot = 0.0
        else:
            stage = self.stages[self.current_stage_index]
            thrust_force, mdot = stage.burn(dt)
            
            # Simple staging logic: if fuel runs out, separate stage immediately
            if stage.fuel_mass <= 0:
                self.current_stage_index += 1
                if self.current_stage_index < len(self.stages):
                    print(f"Stage separation! Ignition of stage {self.current_stage_index + 1} at t={self.time:.2f}s")

        # Physics calculation
        r = np.linalg.norm(self.position)
        altitude = r - R_EARTH
        
        density = self.get_atmospheric_density(altitude)
        gravity_force = self.get_gravity(self.position) * self.current_mass
        drag_force = self.get_drag(self.velocity, density)
        
        # Thrust direction (simplified: always tangent to velocity or vertical initially)
        # For a simple gravity turn, we need a pitch program.
        # Let's implement a basic pitch maneuver.
        if np.linalg.norm(self.velocity) < 1.0:
             # Vertical launch
             thrust_direction = np.array([0.0, 1.0]) # Assuming launch from (0, R_EARTH) means up is +y
             # But wait, position is (0, R_EARTH). So local up is (0, 1).
             # Let's refine the coordinate system: 2D plane through Earth center.
             # Launch site is at (0, R_EARTH).
             thrust_direction = self.position / r
        else:
             # Gravity turn: thrust vector follows velocity vector, but maybe pitch down slightly?
             # A real gravity turn starts with a small pitch over maneuver, then follows gravity.
             # Here we will just follow velocity vector for simplicity after initial pitch over.
             
             # Simple Pitch Program:
             # Vertical until 1km, then pitch over by 2 degree, then gravity turn (thrust aligns with velocity)
             v_unit = self.velocity / np.linalg.norm(self.velocity)
             
             if altitude < 500:
                 thrust_direction = self.position / r # Vertical
             elif altitude < 1500:
                  # Initiate pitch over (nudge velocity vector slightly horizontal)
                  # We can't change velocity directly, we steer thrust.
                  # Let's rotate the "vertical" vector slightly.
                  angle = np.deg2rad(5.0) # 5 degree pitch over
                  # Rotate (0,1) by 1 deg to the right -> (sin(a), cos(a))
                  # Actually we need to rotate relative to local vertical.
                  # Local vertical is position/r.
                  # Let's say we launch from North Pole (0, R) for simplicity of math visualization?
                  # No, let's stick to (0, R) as launch site.
                  # Tangent is (1, 0).
                  # We want to tip towards (1, 0).
                  vertical = self.position / r
                  tangent = np.array([1.0, 0.0]) # East
                  thrust_direction = np.cos(angle)*vertical + np.sin(angle)*tangent
                  thrust_direction /= np.linalg.norm(thrust_direction)
             else:
                  # Gravity turn: thrust aligns with velocity
                  thrust_direction = v_unit

        total_force = (thrust_direction * thrust_force) + gravity_force + drag_force
        acceleration = total_force / self.current_mass

        # Update state (Euler integration)
        self.velocity += acceleration * dt
        self.position += self.velocity * dt
        self.time += dt

        # Calculate dynamic pressure
        q = 0.5 * density * (np.linalg.norm(self.velocity)**2)
        
        # Calculate downrange distance (arc length on Earth surface)
        # Angle between initial position (0, R) and current position
        # theta = arccos( dot(p0, p) / (|p0| |p|) )
        p0 = np.array([0.0, R_EARTH])
        # clip dot product to avoid numerical errors
        dot_prod = np.dot(p0, self.position) / (R_EARTH * r)
        dot_prod = np.clip(dot_prod, -1.0, 1.0)
        theta = np.arccos(dot_prod)
        downrange = theta * R_EARTH

        # Log history
        self.history['time'].append(self.time)
        self.history['altitude'].append(altitude)
        self.history['velocity'].append(np.linalg.norm(self.velocity))
        self.history['acceleration'].append(np.linalg.norm(acceleration))
        self.history['mass'].append(self.current_mass)
        self.history['dynamic_pressure'].append(q)
        self.history['downrange_distance'].append(downrange)

        return altitude
