from .orbit import Orbit
from .frame import Frame
from .trajectory import Trajectory
import matplotlib.pyplot as plt
import numpy as np

class Plotter:
    def __init__(self, plot3d=True):
        """
        Initializes the Plotter class

        :param plot3d: Boolean to plot in 3D
        """
        self.plot3d = plot3d
    
        fig = plt.figure(figsize=(16, 9))
        if plot3d:
            ax = fig.add_subplot(projection='3d')
        else:
            ax = fig.add_subplot()
        self.ax = ax


    ########## Color functions ##########
    @staticmethod
    def color_gradient(value, colors=[(1,0.5,0,1), (0,0,1,1)]):
        """
        Returns an interpolated color based on a number between 0 and 1.

        :param colors: List of RGBA colors to interpolate (each color as tuple of 4 values 0-1)
        :param value: Value between 0 and 1 to interpolate the color 
        :return: Interpolated RGBA color tuple
        """
        if value < 0 or value > 1:
            raise ValueError("Value must be between 0 and 1.")
            
        if len(colors) < 2:
            raise ValueError("At least 2 colors are needed for interpolation.")
            
        # Find the two colors to interpolate between
        num_segments = len(colors) - 1
        segment = value * num_segments
        i = int(segment)
        
        # Handle last segment edge case
        if i >= num_segments:
            return colors[-1]
            
        # Interpolate between the two colors
        t = segment - i
        c1 = np.array(colors[i])
        c2 = np.array(colors[i+1])
        
        return tuple(c1 * (1-t) + c2 * t)

    @staticmethod
    def lighten(color, amount=0.5):
        """
        Lightens a color by mixing it with white

        :param color: RGBA tuple (values between 0 and 1)
        :param amount: Amount of white to add (0 to 1)
        :return: Lightened color
        """
        white = (1, 1, 1, color[3])  # Keep the original alpha
        return tuple(c1 * (1-amount) + c2 * amount for c1, c2 in zip(color, white))

    @staticmethod
    def darken(color, amount=0.5):
        """
        Darkens a color by mixing it with black

        :param color: RGBA tuple (values between 0 and 1)
        :param amount: Amount of black to add (0 to 1)
        :return: Darkened color
        """
        black = (0, 0, 0, color[3])  # Keep the original alpha
        return tuple(c1 * (1-amount) + c2 * amount for c1, c2 in zip(color, black))

    @staticmethod
    def derive_color(color, object = "trajectory"):
        if object == "trajectory":
            return (color[0]*0.1, color[1]*0.1, color[2]*0.1, 0.1)
        elif object == "position":
            return (color[0]*0.5, color[1]*0.5, color[2]*0.5, 0.5)
        elif object == "h_vec":
            return Plotter.darken(color)
        elif object == "e_vec":
            return color
        elif object == "n_vec":
            return Plotter.lighten(color)


    ########## Draw functions ##########
    def draw_body(self, body, r_vec_bc, t_clock_bc, equator=True, label=None, color=(1,0,0,0.4)):        
        # Body surface
        u = np.linspace(0, 2 * np.pi, 36)
        v = np.linspace(0, np.pi, 36)
        sphere_vec_bc = np.array([r_vec_bc[0] + body.radius * np.outer(np.cos(u), np.sin(v)),
                               r_vec_bc[1] + body.radius * np.outer(np.sin(u), np.sin(v)),
                               r_vec_bc[2] + body.radius * np.outer(np.ones(np.size(u)), np.cos(v))])
        sphere_vec = self.frame.transform_bc_to_frame(sphere_vec_bc, t_clock_bc)
        if equator:
            # Equator line
            phi = np.linspace(0, 2 * np.pi, 36)
            equator_vec_bc = np.array([r_vec_bc[0] + body.radius * np.cos(phi),
                                       r_vec_bc[1] + body.radius * np.sin(phi),
                                       r_vec_bc[2] * np.ones(np.size(phi))])
            equator_vec = self.frame.transform_bc_to_frame(equator_vec_bc, t_clock_bc)
            #plot equator line
            if self.plot3d:
                self.ax.plot(equator_vec[0], equator_vec[1], equator_vec[2], color='gray')
            else:
                self.ax.plot(equator_vec[0], equator_vec[1], color='gray')
        
        #plot body surface
        if label is None:
            label = body.name
        if self.plot3d:
            self.ax.plot_surface(sphere_vec[0], sphere_vec[1], sphere_vec[2], color=color, label=label, zorder=2)   # type: ignore
        else:
            center_vec = self.frame.transform_bc_to_frame(np.array([0, 0, 0]), t_clock_bc)
            circle = plt.Circle((center_vec[0], center_vec[1]), body.radius, color=color, label=label, zorder=2) # type: ignore
            self.ax.add_patch(circle)

    def draw_points(self, orbit):
        if orbit.points_in_orbital_plane:
            for i, point in enumerate(orbit.points_in_orbital_plane):
                label = point['name']
                # List of colors
                colors = ['red', 'purple', 'brown', 'deeppink', 'olive', 'deepskyblue', 'green']
                color = colors[i % len(colors)]
                r_vec_bc = orbit.state_vectors_at_theta(point['theta'], "bodycentric")[0]
                r_vec = self.frame.transform_bc_to_frame(r_vec_bc, point['t_clock'])
                if self.plot3d:
                    self.ax.plot(r_vec[0], r_vec[1], r_vec[2], 'o', color=color, label=label, zorder=2)
                else:
                    self.ax.plot(r_vec[0], r_vec[1], 'o', color=color, label=label, zorder=2)
    
    def draw_positions(self, orbit, velocities=True, groundtrack=False, color=(1,0.5,0,1), v_scale=1):
        if orbit.orbital_positions:
            t_min = min(position['t_clock'] for position in orbit.orbital_positions if position['t_clock'] is not None)
            t_max = max(position['t_clock'] for position in orbit.orbital_positions if position['t_clock'] is not None)
            
            for i, position in enumerate(orbit.orbital_positions):
                label = position.get('name', f'Position {i+1}')
                r_vec_bc, v_vec_bc = orbit.state_vectors_at_theta(position['theta'], "bodycentric")
                r_vec = self.frame.transform_bc_to_frame(r_vec_bc, position['t_clock'])
                v_vec = self.frame.transform_bc_to_frame(v_vec_bc, position['t_clock'])
                label += f'({position["t_clock"]:.2f}s)'
                
                if t_max == t_min:
                    t_norm = 1
                else:
                    t_norm = (position['t_clock'] - t_min) / (t_max - t_min)
                point_color = self.color_gradient(t_norm, [self.derive_color(color, object="position"), color])
                
                # Plot a marker at the position
                if self.plot3d:
                    self.ax.plot(r_vec[0], r_vec[1], r_vec[2], 'o', color=point_color, label=label)
                else:
                    self.ax.plot(r_vec[0], r_vec[1], 'o', color=point_color, label=label, zorder=2)
                
                # Plot velocity vector
                if velocities:
                    v_vec = v_vec * orbit.main_body.radius / 20 * v_scale  # Velocity vector scale
                    if groundtrack:
                        v_r = np.dot(v_vec, r_vec/np.linalg.norm(r_vec))
                        v_vec = v_vec - v_r * r_vec/np.linalg.norm(r_vec)
                    if self.plot3d:
                        self.ax.quiver(r_vec[0], r_vec[1], r_vec[2], v_vec[0], v_vec[1], v_vec[2], color=point_color, zorder=5)
                    else:
                        self.ax.quiver(r_vec[0], r_vec[1], v_vec[0], v_vec[1], color=point_color, scale_units='xy', scale=1, width=0.003, zorder=5)

    def draw_orbit(self, orbit, t_clock_bc=0, groundtrack=False, h_arrow=True, e_arrow=True, n_arrow=True, style="solid", color=(0,0,1,1), label=None):
        """
        Draws the orbit

        :param orbit: Orbit object
        :param t_clock_bc: t_clock instant when plotting origin and auxiliary vectors
        :param groundtrack: Boolean to plot groundtrack
        """

        h_vec = self.frame.transform_bc_to_frame(orbit.h_vec_bc, t_clock_bc)
        e_vec = self.frame.transform_bc_to_frame(orbit.e_vec_bc, t_clock_bc)
        n_vec = self.frame.transform_bc_to_frame(orbit.n_vec_bc, t_clock_bc)
        
        ### Scaling vectors for plotting ###
        if orbit.h == 0:
            h_vec = [0, 0, 0]
        else:
            h_vec = h_vec/orbit.h * orbit.rp
        if orbit.n == 0:
            n_vec = [0, 0, 0]
        else:
            theta_n = orbit.theta_ascending_node()
            r_n = orbit.r_at_theta(theta_n)
            n_vec = n_vec/orbit.n * r_n
        if orbit.e < 1e-7:
            e_vec = [0, 0, 0]
        else:
            e_vec = e_vec/orbit.e * orbit.rp

        # Create array of true anomaly
        if orbit.e == 1:
            theta_min = -120
            theta_max = 120
            if orbit.orbital_positions:
                for pos in orbit.orbital_positions:
                    if pos['theta'] < theta_min:
                        theta_min = pos['theta']
                    if pos['theta'] > theta_max:
                        theta_max = pos['theta']
            theta_list = np.linspace(theta_min, theta_max, 360)
        elif orbit.e > 1:
            epsilon = 15
            theta_min = -orbit.theta_infinity() + epsilon
            theta_max = orbit.theta_infinity() - epsilon
            if orbit.orbital_positions:
                for pos in orbit.orbital_positions:
                    if pos['theta'] < theta_min:
                        theta_min = pos['theta']
                    if pos['theta'] > theta_max:
                        theta_max = pos['theta']
            theta_list = np.linspace(theta_min, theta_max, 360)
        else:
            theta_list = np.linspace(0, 360, 360)
        
        r_vecs = []
        for theta in theta_list:
            t_clock = orbit.t_clock_at_theta(theta)
            r_vec_bc = orbit.state_vectors_at_t_clock(t_clock, "bodycentric")[0]
            r_vec = self.frame.transform_bc_to_frame(r_vec_bc, t_clock)
            r_vecs.append(r_vec)
        r_vecs = np.array(r_vecs).T
        
        if groundtrack:
            r_vecs = r_vecs / np.linalg.norm(r_vecs, axis=1, keepdims=True) * orbit.main_body.radius
        
        # Plot orbit
        center_vec = self.frame.transform_bc_to_frame(np.array([0, 0, 0]), t_clock_bc)
        if style == "solid":
            style = '-'
        elif style == "dashed": 
            style = '--'
        if label is None:
            label = orbit.name
        if self.plot3d:
            self.ax.plot(r_vecs[0], r_vecs[1], r_vecs[2], style, color=color, label=label, linewidth=2, zorder=1)
            if h_arrow:
                self.ax.quiver(center_vec[0], center_vec[1], center_vec[2], h_vec[0], h_vec[1], h_vec[2], color=self.derive_color(color, object="h_vec"), label='h - Angular Momentum')
            if e_arrow:
                self.ax.quiver(center_vec[0], center_vec[1], center_vec[2], e_vec[0], e_vec[1], e_vec[2], color=self.derive_color(color, object="e_vec"), label='e - Eccentricity')
            if n_arrow:
                self.ax.quiver(center_vec[0], center_vec[1], center_vec[2], n_vec[0], n_vec[1], n_vec[2], color=self.derive_color(color, object="n_vec"), label='n - Ascending Node')
        else:
            self.ax.plot(r_vecs[0], r_vecs[1], style, color=color, label=label, linewidth=2, zorder=1)
            if h_arrow:
                self.ax.quiver(center_vec[0], center_vec[1], h_vec[0], h_vec[1], color=self.derive_color(color, object="h_vec"), label='h', scale_units='xy', scale=1, width=0.003)
            if e_arrow:
                self.ax.quiver(center_vec[0], center_vec[1], e_vec[0], e_vec[1], color=self.derive_color(color, object="e_vec"), label='e', scale_units='xy', scale=1, width=0.003)
            if n_arrow:
                self.ax.quiver(center_vec[0], center_vec[1], n_vec[0], n_vec[1], color=self.derive_color(color, object="n_vec"), label='n', scale_units='xy', scale=1, width=0.003)

    def draw_trajectory(self, orbit_number, traj_idx=0, time_step=300, positions=True, velocities=True, groundtrack=False, color=(1,0.5,0,1), v_scale=1):
        """
        Draws the trajectory

        :param orbit_number: Orbit number
        :param traj_idx: Trajectory index
        :param groundtrack: Boolean to plot groundtrack
        """

        orbit = self.trajectories[traj_idx].orbits[orbit_number]['orbit']
        t1_clock = self.trajectories[traj_idx].get_trajectory_position(orbit_number, position_index="last")['t_clock']
        t0_clock = self.trajectories[traj_idx].get_trajectory_position(orbit_number, position_index="first")['t_clock']
        annotations = []  # List to store all annotations
        if positions:
            for position in self.trajectories[traj_idx].orbits[orbit_number]['trajectory_positions']:
                r_vec, v_vec = orbit.state_vectors_at_t_clock(position['t_clock'], self.frame)
                if t1_clock == t0_clock:
                    t_norm = 0
                else:
                    t_norm = (position['t_clock'] - t0_clock) / (t1_clock - t0_clock)
                point_color = self.color_gradient(t_norm, colors=[self.derive_color(color, object="position"), color])
                if groundtrack:
                    r_vec = r_vec / np.linalg.norm(r_vec) * orbit.main_body.radius
                if self.plot3d:
                    point_plot = self.ax.plot(r_vec[0], r_vec[1], r_vec[2], 'o', color=point_color,
                                        picker=True, pickradius=5, label=position['name']+f' ({position["t_clock"]:.2f}s)', zorder=4)[0]
                else:
                    point_plot = self.ax.plot(r_vec[0], r_vec[1], 'o', color=point_color,
                                        picker=True, pickradius=5, label=position['name']+f' ({position["t_clock"]:.2f}s)', zorder=4)[0]
                    # Create annotation
                    annotation = self.ax.annotate(f't = {position["t_clock"]:.2f}s',
                        xy=(r_vec[0], r_vec[1]), xytext=(10, 10),
                        textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                        arrowprops=dict(arrowstyle='->'),
                        visible=False,
                        zorder=10)
                    annotations.append((point_plot, annotation))
                # Plot velocity vector
                if velocities:
                    v_vec = v_vec * orbit.main_body.radius / 20 * v_scale # Velocity vector scale
                    if groundtrack:
                        v_r = np.dot(v_vec, r_vec/np.linalg.norm(r_vec))
                        v_vec = v_vec - v_r * r_vec/np.linalg.norm(r_vec)
                    if self.plot3d:
                        self.ax.quiver(r_vec[0], r_vec[1], r_vec[2], v_vec[0], v_vec[1], v_vec[2], color=point_color, zorder=5)
                    else:
                        self.ax.quiver(r_vec[0], r_vec[1], v_vec[0], v_vec[1], color=point_color, scale_units='xy', scale=1, width=0.003, zorder=5)

    
        # Hover function that manages all annotations
        def hover(event):
            if event.inaxes == plt.gca():
                for point_plot, annotation in annotations:
                    cont, _ = point_plot.contains(event)
                    annotation.set_visible(cont)
                plt.draw()
        
        plt.gcf().canvas.mpl_connect('motion_notify_event', hover)
        
        # Generate high resolution points for the continuous line
        times = np.arange(t0_clock, t1_clock, time_step)
        # Add the final point if it's not already included
        if t1_clock not in times:
            times = np.append(times, t1_clock)
        
        r_vecs = []
        colors = []
        for t in times:
            r_vec = orbit.state_vectors_at_t_clock(t, frame=self.frame)[0]
            if groundtrack:
                r_vec = r_vec / np.linalg.norm(r_vec) * orbit.main_body.radius
                
            if t0_clock == t1_clock:
                t_norm = 1
            else:
                t_norm = (t - t0_clock) / (t1_clock - t0_clock)
            r_vecs.append(r_vec)
            colors.append(self.color_gradient(t_norm, colors=[self.derive_color(color, object="trajectory"), color]))

        # Convert to numpy array for better performance
        r_vecs = np.array(r_vecs)
        for i in range(len(r_vecs)-1):
            if self.plot3d:
                self.ax.plot(r_vecs[i:i+2,0], r_vecs[i:i+2,1], r_vecs[i:i+2,2], 
                        '-', color=colors[i], alpha=0.6, linewidth=2, zorder=3)
            else:
                self.ax.plot(r_vecs[i:i+2,0], r_vecs[i:i+2,1], 
                            '-', color=colors[i], alpha=0.6, linewidth=2, zorder=3)

    def draw_maneuver(self, orbit, maneuver, color=(0.7,0,0), v_scale=1):
        """
        Draws the maneuver
        """
        label = maneuver.name
        r_vec_bc = orbit.state_vectors_at_t_clock(maneuver.t_clock, "bodycentric")[0]
        r_vec = self.frame.transform_bc_to_frame(r_vec_bc, maneuver.t_clock)
        label += f' ({maneuver.t_clock:.2f}s)'
        
        if self.plot3d:
            self.ax.plot(r_vec[0], r_vec[1], r_vec[2], 'H', color=color, label=label, zorder=1, markersize=15)
        else:
            self.ax.plot(r_vec[0], r_vec[1], 'H', color=color, label=label, zorder=1, markersize=15)
        
        if maneuver.delta_v_vec_bc is not None:
            delta_v_vec = self.frame.transform_bc_to_frame(maneuver.delta_v_vec_bc, maneuver.t_clock)
            delta_v_vec = delta_v_vec * orbit.main_body.radius/20 * v_scale  # Velocity vector scale
            if self.plot3d:
                self.ax.quiver(r_vec[0], r_vec[1], r_vec[2], delta_v_vec[0], delta_v_vec[1], delta_v_vec[2], color=color, zorder=12)
            else:
                self.ax.quiver(r_vec[0], r_vec[1], delta_v_vec[0], delta_v_vec[1], color=color, scale_units='xy', scale=1, width=0.003, zorder=12)


    ########## Plot functions ##########
    def plot_orbit(self, orbit, frame="bodycentric", points=True, positions=True, velocities=True, groundtrack=False, v_scale=1, h_arrow=True, e_arrow=True, n_arrow=True):
        """
        Plots the orbit
        """
        if not isinstance(orbit, Orbit):
            raise ValueError("The 'orbit' parameter must be an instance of the Orbit class")
        
        if isinstance(frame, str):
            if frame == "perifocal":
                frame = Frame(name="perifocal",
                              Omega0_bc_frame=orbit.Omega0,
                              Omega_dot_bc_frame=orbit.Omega_dot,
                              omega0_bc_frame=orbit.omega0,
                              omega_dot_bc_frame=orbit.omega_dot,
                              i0_bc_frame=orbit.i,
                              t_clock_bc_frame=orbit.t0_clock)
            elif frame == "perifocal_t0":
                frame = Frame(name="perifocal_t0",
                              Omega0_bc_frame=orbit.Omega0,
                              Omega_dot_bc_frame=0,
                              omega0_bc_frame=orbit.omega0,
                              omega_dot_bc_frame=0,
                              i0_bc_frame=orbit.i,
                              t_clock_bc_frame=orbit.t0_clock)
            elif frame == "bodycentric":
                frame = Frame.bodycentric()
            elif frame == "rotating_bodycentric":
                frame = Frame.rotating_bodycentric(orbit.main_body)
            else:
                raise ValueError("Invalid frame name. Please use 'perifocal', 'perifocal_t0', 'bodycentric', 'rotating_bodycentric' or a custom Frame object.")
        if not isinstance(frame, Frame):
            raise ValueError("The 'frame' parameter must be an instance of the Frame class")
        self.frame = frame

        self.draw_body(orbit.main_body, r_vec_bc=np.array([0, 0, 0]), t_clock_bc=orbit.t0_clock)
        self.draw_orbit(orbit, t_clock_bc=orbit.t0_clock, groundtrack=groundtrack, h_arrow=h_arrow, e_arrow=e_arrow, n_arrow=n_arrow)
        if points:
            self.draw_points(orbit)
        if positions:
            self.draw_positions(orbit, velocities=velocities, groundtrack=groundtrack, v_scale=v_scale)

        self.ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        self.ax.set_xlabel('X (km)')
        self.ax.set_ylabel('Y (km)')
        if self.plot3d:
            self.ax.set_zlabel('Z (km)') # type: ignore
        self.ax.axis('equal')
        self.ax.grid(True)
        plt.title(f'Frame: {self.frame.name}')
        plt.tight_layout()
        plt.get_current_fig_manager().window.showMaximized() # type: ignore
        plt.show()

    def plot_trajectory(self, trajectory, time_step=300, frame="bodycentric", orbits=True, orbit_arrows=False, points=True, positions=True, velocities=True, maneuvers=True, groundtrack=False, v_scale=1):
        """
        Plots the trajectory
        """
        if not isinstance(trajectory, Trajectory):
            raise ValueError("The 'trajectory' parameter must be an instance of the Trajectory class")
        self.trajectories = [trajectory]

        # Use first orbit to define frame
        orbit0 = trajectory.orbits[0]['orbit']
        if isinstance(frame, str):
            if frame == "perifocal":
                frame = Frame(name="perifocal",
                              Omega0_bc_frame=orbit0.Omega0,
                              Omega_dot_bc_frame=orbit0.Omega_dot,
                              omega0_bc_frame=orbit0.omega0,
                              omega_dot_bc_frame=orbit0.omega_dot,
                              i0_bc_frame=orbit0.i,
                              t_clock_bc_frame=orbit0.t0_clock)
            elif frame == "perifocal_t0":
                frame = Frame(name="perifocal_t0",
                              Omega0_bc_frame=orbit0.Omega0,
                              Omega_dot_bc_frame=0,
                              omega0_bc_frame=orbit0.omega0,
                              omega_dot_bc_frame=0,
                              i0_bc_frame=orbit0.i,
                              t_clock_bc_frame=orbit0.t0_clock)
            elif frame == "bodycentric":
                frame = Frame.bodycentric()
            elif frame == "rotating_bodycentric":
                frame = Frame.rotating_bodycentric(orbit0.main_body)
            else:
                raise ValueError("Invalid frame name. Please use 'perifocal', 'perifocal_t0', 'bodycentric', 'rotating_bodycentric' or a custom Frame object.")
        if not isinstance(frame, Frame):
            raise ValueError("The 'frame' parameter must be an instance of the Frame class")
        self.frame = frame
        self.draw_body(orbit0.main_body, r_vec_bc=np.array([0, 0, 0]), t_clock_bc=orbit0.t0_clock)

        colors = [(1.0, 0.5, 0.0, 1),
                  (0.0, 0.0, 1.0, 1),
                  (1.0, 0.0, 0.2, 1),
                  (0.0, 0.5, 0.8, 1),
                  (0.9, 0.9, 0.0, 1)]
        
        for orbit_number in trajectory.orbits:
            orbit = trajectory.orbits[orbit_number]['orbit']
            if orbits:
                if orbit_arrows:
                    h_arrow = True
                    e_arrow = True
                    n_arrow = True
                else:
                    h_arrow = False
                    e_arrow = False
                    n_arrow = False
                self.draw_orbit(orbit, t_clock_bc=orbit.t0_clock, groundtrack=groundtrack, h_arrow=h_arrow, e_arrow=e_arrow, n_arrow=n_arrow, style="dashed", color=colors[orbit_number%len(colors)], label=f"{orbit_number} - {orbit.name}")
            if points:
                self.draw_points(orbit)
            self.draw_trajectory(orbit_number, time_step=time_step, positions=positions, velocities=velocities, groundtrack=groundtrack, color=colors[orbit_number%len(colors)], v_scale=v_scale)
            if maneuvers:
                for maneuver in trajectory.orbits[orbit_number]['maneuvers']:
                    self.draw_maneuver(orbit, maneuver, v_scale=v_scale)
        
        self.ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        self.ax.set_xlabel('X (km)')
        self.ax.set_ylabel('Y (km)')
        if self.plot3d:
            self.ax.set_zlabel('Z (km)') # type: ignore
        self.ax.axis('equal')
        self.ax.grid(True)
        plt.title(f'Frame: {self.frame.name}')
        plt.tight_layout()
        plt.get_current_fig_manager().window.showMaximized() # type: ignore
        plt.show()

    def plot_groundtrack(self, trajectory, time_step=300, frame="rotating_bodycentric"):
        """
        Plots the groundtrack
        """
        if not isinstance(trajectory, Trajectory):
            raise ValueError("The 'trajectory' parameter must be an instance of the Trajectory class")

        # Use first orbit to define frame
        orbit0 = trajectory.orbits[0]['orbit']
        if isinstance(frame, str):
            if frame == "perifocal":
                frame = Frame(name="perifocal",
                              Omega0_bc_frame=orbit0.Omega0,
                              Omega_dot_bc_frame=orbit0.Omega_dot,
                              omega0_bc_frame=orbit0.omega0,
                              omega_dot_bc_frame=orbit0.omega_dot,
                              i0_bc_frame=orbit0.i,
                              t_clock_bc_frame=orbit0.t0_clock)
            elif frame == "perifocal_t0":
                frame = Frame(name="perifocal_t0",
                              Omega0_bc_frame=orbit0.Omega0,
                              Omega_dot_bc_frame=0,
                              omega0_bc_frame=orbit0.omega0,
                              omega_dot_bc_frame=0,
                              i0_bc_frame=orbit0.i,
                              t_clock_bc_frame=orbit0.t0_clock)
            elif frame == "bodycentric":
                frame = Frame.bodycentric()
            elif frame == "rotating_bodycentric":
                frame = Frame.rotating_bodycentric(orbit0.main_body)
            else:
                raise ValueError("Invalid frame name. Please use 'perifocal', 'perifocal_t0', 'bodycentric', 'rotating_bodycentric' or a custom Frame object.")
        if not isinstance(frame, Frame):
            raise ValueError("The 'frame' parameter must be an instance of the Frame class")
        
        fig_gt = plt.figure(figsize=(16, 9))
        ax_gt = fig_gt.add_subplot()

        annotations = []  # List to store all annotations        
        
        for orbit_number in trajectory.orbits:
            orbit = trajectory.orbits[orbit_number]['orbit']
            t1_clock = trajectory.get_trajectory_position(orbit_number, position_index="last")['t_clock']
            t0_clock = trajectory.get_trajectory_position(orbit_number, position_index="first")['t_clock']
            for position in trajectory.orbits[orbit_number]['trajectory_positions']:
                r_vec = orbit.state_vectors_at_t_clock(position['t_clock'], frame=frame)[0]
                ra, dec = frame.convert_cartesian_to_ra_dec(r_vec)
                t_clock = position['t_clock']
                t_norm = (t_clock - t0_clock) / (t1_clock - t0_clock)
                color = self.color_gradient(t_norm)
                point_plot = ax_gt.plot(ra, dec, 'o', color=color,
                                        picker=True, pickradius=5,
                                        markersize=5, zorder=4)[0]
                    
                # Create annotation
                annotation = ax_gt.annotate(f't = {t_clock:.2f}s',
                    xy=(ra, dec), xytext=(10, 10),
                    textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                    arrowprops=dict(arrowstyle='->'),
                    visible=False,
                    zorder=10)
                
                annotations.append((point_plot, annotation))
            
            # Hover function that manages all annotations
            def hover(event):
                if event.inaxes == plt.gca():
                    for point_plot, annotation in annotations:
                        cont, _ = point_plot.contains(event)
                        annotation.set_visible(cont)
                    plt.draw()
            
            fig_gt.canvas.mpl_connect('motion_notify_event', hover)
            
            ##### Continuous line #####
            samples = int((t1_clock - t0_clock)/time_step)
            times = np.linspace(t0_clock, t1_clock, samples)
            for t in times:
                r_vec = orbit.state_vectors_at_t_clock(t, frame=frame)[0]
                ra, dec = frame.convert_cartesian_to_ra_dec(r_vec)
                t_norm = (t - t0_clock) / (t1_clock - t0_clock)
                color = self.color_gradient(t_norm)
                ax_gt.plot(ra, dec, '.', color=color,
                        alpha=0.6, markersize=3, zorder=3)

        ax_gt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        ax_gt.set_xlabel('Right Ascension (°)')
        ax_gt.set_ylabel('Declination (°)')
        ax_gt.set_xlim(-180, 180)
        ax_gt.set_ylim(-90, 90)
        ax_gt.grid(True)
        ax_gt.axhline(y=0, color='k', linestyle='-', linewidth=1)
        ax_gt.axvline(x=0, color='k', linestyle='-', linewidth=1)

        # Add text for coordinates
        coord_text = ax_gt.text(0.02, 0.98, '', transform=ax_gt.transAxes, 
                            bbox=dict(facecolor='white', alpha=0.7),
                            verticalalignment='top')

        # Show coordinates on hover
        def update_coords(event):
            if event.inaxes == ax_gt:
                ra = event.xdata
                dec = event.ydata
                coord_text.set_text(f'RA: {ra:.2f}°\nDec: {dec:.2f}°')
                plt.draw()

        plt.gcf().canvas.mpl_connect('motion_notify_event', update_coords)

        plt.title(f'Groundtrack')
        plt.tight_layout()
        plt.get_current_fig_manager().window.showMaximized() # type: ignore
        plt.show()

    def plot_orbits(self, orbits, frame="bodycentric", points=True, positions=True, velocities=True, groundtrack=False, v_scale=1, h_arrow=True, e_arrow=True, n_arrow=True):
        """
        Plots the orbits
        """
        # Convert single input to list for uniform processing
        if isinstance(orbits, Orbit):
            orbits = [orbits]
        
        # Use first orbit of first trajectory to define frame
        orbit0 = orbits[0]
        if isinstance(frame, str):
            if frame == "perifocal":
                frame = Frame(name="perifocal",
                              Omega0_bc_frame=orbit0.Omega0,
                              Omega_dot_bc_frame=orbit0.Omega_dot,
                              omega0_bc_frame=orbit0.omega0,
                              omega_dot_bc_frame=orbit0.omega_dot,
                              i0_bc_frame=orbit0.i,
                              t_clock_bc_frame=orbit0.t0_clock)
            elif frame == "perifocal_t0":
                frame = Frame(name="perifocal_t0",
                              Omega0_bc_frame=orbit0.Omega0,
                              Omega_dot_bc_frame=0,
                              omega0_bc_frame=orbit0.omega0,
                              omega_dot_bc_frame=0,
                              i0_bc_frame=orbit0.i,
                              t_clock_bc_frame=orbit0.t0_clock)
            elif frame == "bodycentric":
                frame = Frame.bodycentric()
            elif frame == "rotating_bodycentric":
                frame = Frame.rotating_bodycentric(orbit0.main_body)
            else:
                raise ValueError("Invalid frame name. Please use 'perifocal', 'perifocal_t0', 'bodycentric', 'rotating_bodycentric' or a custom Frame object.")
        if not isinstance(frame, Frame):
            raise ValueError("The 'frame' parameter must be an instance of the Frame class")
        self.frame = frame

        # Draw central body once
        self.draw_body(orbit0.main_body, r_vec_bc=np.array([0, 0, 0]), t_clock_bc=orbit0.t0_clock)

        # Different colors for each orbit
        orbit_colors = [(0.0, 0.0, 1.0, 1),
                        (1.0, 0.8, 0.0, 1),
                        (0.5, 0.0, 0.8, 1),
                        (1.0, 0.0, 0.0, 1),
                        (0.5, 0.5, 0.5, 1)]
        
        # For each orbit
        for orbit_number, orbit in enumerate(orbits):
            self.draw_orbit(orbit, t_clock_bc=orbit.t0_clock, groundtrack=groundtrack, color=orbit_colors[orbit_number%len(orbit_colors)], label=f"{orbit_number} - {orbit.name}", h_arrow=h_arrow, e_arrow=e_arrow, n_arrow=n_arrow)
            if points:
                self.draw_points(orbit)
            if positions:
                self.draw_positions(orbit, velocities=velocities, groundtrack=groundtrack, v_scale=v_scale)

        self.ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        self.ax.set_xlabel('X (km)')
        self.ax.set_ylabel('Y (km)')
        if self.plot3d:
            self.ax.set_zlabel('Z (km)') # type: ignore
        self.ax.axis('equal')
        self.ax.grid(True)
        plt.title(f'Frame: {self.frame.name}')
        plt.tight_layout()
        plt.get_current_fig_manager().window.showMaximized() # type: ignore
        plt.show()

    def plot_trajectories(self, trajectories, time_step=300, frame="bodycentric", orbits=True, points=True, positions=True, velocities=True, maneuvers=True, groundtrack=False, v_scale=1):
        """
        Plots a single or multiple trajectories
        
        Args:
            trajectories: A single trajectory or a list of trajectories to plot
            time_step: Time step for continuous line
            frame: Reference frame for plotting ("bodycentric", "perifocal", etc)
            orbits: If True, plots the orbits
            points: If True, plots the points of interest
            positions: If True, plots the positions along the trajectory
            velocities: If True, plots the velocity vectors
            maneuvers: If True, plots the maneuvers
            groundtrack: If True, plots the groundtrack
            v_scale: Scale factor for the velocity vectors
        """
        # Convert single input to list for uniform processing
        if isinstance(trajectories, Trajectory):
            trajectories = [trajectories]
            
        if not all(isinstance(traj, Trajectory) for traj in trajectories):
            raise ValueError("All elements must be instances of the Trajectory class")

        self.trajectories = trajectories

        # Use first orbit of first trajectory to define frame
        orbit0 = trajectories[0].orbits[0]['orbit']
        if isinstance(frame, str):
            if frame == "perifocal":
                frame = Frame(name="perifocal",
                              Omega0_bc_frame=orbit0.Omega0,
                              Omega_dot_bc_frame=orbit0.Omega_dot,
                              omega0_bc_frame=orbit0.omega0,
                              omega_dot_bc_frame=orbit0.omega_dot,
                              i0_bc_frame=orbit0.i,
                              t_clock_bc_frame=orbit0.t0_clock)
            elif frame == "perifocal_t0":
                frame = Frame(name="perifocal_t0",
                              Omega0_bc_frame=orbit0.Omega0,
                              Omega_dot_bc_frame=0,
                              omega0_bc_frame=orbit0.omega0,
                              omega_dot_bc_frame=0,
                              i0_bc_frame=orbit0.i,
                              t_clock_bc_frame=orbit0.t0_clock)
            elif frame == "bodycentric":
                frame = Frame.bodycentric()
            elif frame == "rotating_bodycentric":
                frame = Frame.rotating_bodycentric(orbit0.main_body)
            else:
                raise ValueError("Invalid frame name. Use 'perifocal', 'perifocal_t0', 'bodycentric', 'rotating_bodycentric' or a custom Frame object.")
        if not isinstance(frame, Frame):
            raise ValueError("The 'frame' parameter must be an instance of the Frame class")
        self.frame = frame
        
        # Draw central body once
        self.draw_body(orbit0.main_body, r_vec_bc=np.array([0, 0, 0]), t_clock_bc=orbit0.t0_clock)

        # Different colors for each trajectory
        trajectory_colors = [(0.0, 0.0, 1.0, 1),
                             (1.0, 0.8, 0.0, 1),
                             (0.5, 0.0, 0.8, 1),
                             (1.0, 0.0, 0.0, 1),
                             (0.5, 0.5, 0.5, 1)]
        
        # For each trajectory
        for traj_idx, trajectory in enumerate(trajectories):
            for orbit_number in trajectory.orbits:
                orbit_norm = orbit_number/len(trajectory.orbits)
                orbit_color = self.color_gradient(orbit_norm, colors=[trajectory_colors[traj_idx], self.lighten(trajectory_colors[traj_idx], 0.7)])
                orbit = trajectory.orbits[orbit_number]['orbit']
                if orbits:
                    self.draw_orbit(orbit, t_clock_bc=orbit.t0_clock, groundtrack=groundtrack, 
                                  h_arrow=False, e_arrow=True, n_arrow=True, style="dashed", 
                                  color=orbit_color, 
                                  label=f"Trajectory {traj_idx+1} - Orbit {orbit_number}")
                if points:
                    self.draw_points(orbit)
                self.draw_trajectory(orbit_number, traj_idx=traj_idx, time_step=time_step, velocities=velocities, 
                                  groundtrack=groundtrack, color=orbit_color, positions=positions, v_scale=v_scale)
                if maneuvers:
                    for maneuver in trajectory.orbits[orbit_number]['maneuvers']:
                        self.draw_maneuver(orbit, maneuver, v_scale=v_scale)
        
        self.ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        self.ax.set_xlabel('X (km)')
        self.ax.set_ylabel('Y (km)')
        if self.plot3d:
            self.ax.set_zlabel('Z (km)') # type: ignore
        self.ax.axis('equal')
        self.ax.grid(True)
        plt.title(f'Frame: {self.frame.name}')
        plt.tight_layout()
        plt.get_current_fig_manager().window.showMaximized() # type: ignore
        plt.show()

    def plot_groundtracks(self, trajectories, time_step=300, frame="rotating_bodycentric"):
        """
        Plots the groundtrack for multiple trajectories

        :param trajectories: List of Trajectory objects or single Trajectory object
        :param time_step: Time step for continuous line
        :param frame: Reference frame to use ('perifocal', 'perifocal_t0', 'bodycentric', 'rotating_bodycentric' or Frame object)
        """
        # Convert single trajectory to list for uniform handling
        if isinstance(trajectories, Trajectory):
            trajectories = [trajectories]
        
        # Validate trajectories
        if not all(isinstance(traj, Trajectory) for traj in trajectories):
            raise ValueError("All trajectories must be instances of the Trajectory class")

        # Use first orbit of first trajectory to define frame
        orbit0 = trajectories[0].orbits[0]['orbit']
        if isinstance(frame, str):
            if frame == "perifocal":
                frame = Frame(name="perifocal",
                              Omega0_bc_frame=orbit0.Omega0,
                              Omega_dot_bc_frame=orbit0.Omega_dot,
                              omega0_bc_frame=orbit0.omega0,
                              omega_dot_bc_frame=orbit0.omega_dot,
                              i0_bc_frame=orbit0.i,
                              t_clock_bc_frame=orbit0.t0_clock)
            elif frame == "perifocal_t0":
                frame = Frame(name="perifocal_t0",
                              Omega0_bc_frame=orbit0.Omega0,
                              Omega_dot_bc_frame=0,
                              omega0_bc_frame=orbit0.omega0,
                              omega_dot_bc_frame=0,
                              i0_bc_frame=orbit0.i,
                              t_clock_bc_frame=orbit0.t0_clock)
            elif frame == "bodycentric":
                frame = Frame.bodycentric()
            elif frame == "rotating_bodycentric":
                frame = Frame.rotating_bodycentric(orbit0.main_body)
            else:
                raise ValueError("Invalid frame name. Please use 'perifocal', 'perifocal_t0', 'bodycentric', 'rotating_bodycentric' or a custom Frame object.")
        if not isinstance(frame, Frame):
            raise ValueError("The 'frame' parameter must be an instance of the Frame class")
        
        fig_gt = plt.figure(figsize=(16, 9))
        ax_gt = fig_gt.add_subplot()

        annotations = []  # List to store all annotations        
        
        # Different colors for each trajectory
        trajectory_colors = [(0.0, 0.0, 1.0, 1),
                             (1.0, 0.8, 0.0, 1),
                             (0.5, 0.0, 0.8, 1),
                             (1.0, 0.0, 0.0, 1),
                             (0.5, 0.5, 0.5, 1)]
        
        # Plot each trajectory
        for traj_idx, trajectory in enumerate(trajectories):
            for orbit_number in trajectory.orbits:
                orbit = trajectory.orbits[orbit_number]['orbit']
                t1_clock = trajectory.get_trajectory_position(orbit_number, position_index="last")['t_clock']
                t0_clock = trajectory.get_trajectory_position(orbit_number, position_index="first")['t_clock']
                
                # Plot trajectory points
                for position in trajectory.orbits[orbit_number]['trajectory_positions']:
                    r_vec = orbit.state_vectors_at_t_clock(position['t_clock'], frame=frame)[0]
                    ra, dec = frame.convert_cartesian_to_ra_dec(r_vec)
                    t_clock = position['t_clock']
                    t_norm = (t_clock - t0_clock) / (t1_clock - t0_clock)
                    point_color = self.color_gradient(t_norm, colors=[trajectory_colors[traj_idx], (0,0,0,0.5)])
                    point_plot = ax_gt.plot(ra, dec, 'o', color=point_color,
                                            picker=True, pickradius=5,
                                            markersize=5, zorder=4,
                                            label=f'Trajectory {traj_idx+1} - Orbit {orbit_number}')[0]
                        
                    # Create annotation
                    annotation = ax_gt.annotate(f't = {t_clock:.2f}s',
                        xy=(ra, dec), xytext=(10, 10),
                        textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                        arrowprops=dict(arrowstyle='->'),
                        visible=False,
                        zorder=10)
                    
                    annotations.append((point_plot, annotation))
                
                # Plot continuous line
                samples = int((t1_clock - t0_clock)/time_step)
                times = np.linspace(t0_clock, t1_clock, samples)
                for t in times:
                    r_vec = orbit.state_vectors_at_t_clock(t, frame=frame)[0]
                    ra, dec = frame.convert_cartesian_to_ra_dec(r_vec)
                    t_norm = (t - t0_clock) / (t1_clock - t0_clock)
                    point_color = self.color_gradient(t_norm, colors=[trajectory_colors[traj_idx], (0,0,0,0.5)])
                    ax_gt.plot(ra, dec, '.', color=point_color,
                            alpha=0.6, markersize=3, zorder=3)

        # Hover function that manages all annotations
        def hover(event):
            if event.inaxes == plt.gca():
                for point_plot, annotation in annotations:
                    cont, _ = point_plot.contains(event)
                    annotation.set_visible(cont)
                plt.draw()
        
        fig_gt.canvas.mpl_connect('motion_notify_event', hover)

        ax_gt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        ax_gt.set_xlabel('Right Ascension (°)')
        ax_gt.set_ylabel('Declination (°)')
        ax_gt.set_xlim(-180, 180)
        ax_gt.set_ylim(-90, 90)
        ax_gt.grid(True)
        ax_gt.axhline(y=0, color='k', linestyle='-', linewidth=1)
        ax_gt.axvline(x=0, color='k', linestyle='-', linewidth=1)

        # Add text for coordinates
        coord_text = ax_gt.text(0.02, 0.98, '', transform=ax_gt.transAxes, 
                            bbox=dict(facecolor='white', alpha=0.7),
                            verticalalignment='top')

        # Show coordinates on hover
        def update_coords(event):
            if event.inaxes == ax_gt:
                ra = event.xdata
                dec = event.ydata
                coord_text.set_text(f'RA: {ra:.2f}°\nDec: {dec:.2f}°')
                plt.draw()

        plt.gcf().canvas.mpl_connect('motion_notify_event', update_coords)

        plt.title('Groundtrack')
        plt.tight_layout()
        plt.get_current_fig_manager().window.showMaximized() # type: ignore
        plt.show()
