from astronautica.orbit import Orbit
from astronautica.frame import Frame
import matplotlib.pyplot as plt
import numpy as np

class Plotter:
    def __init__(self, frame="bodycentric", plot3d=True):
        """
        Initializes the Plotter class

        :param frame: Frame of reference to plot the orbit
        :param plot3d: Boolean to plot in 3D
        """
        if not isinstance(frame, Frame):
            if frame == "perifocal":
                frame = Frame.perifocal()
            elif frame == "perifocal_t0":
                frame = Frame.perifocal_t0()
            elif frame == "bodycentric":
                frame = Frame.bodycentric()
            else:
                raise ValueError("Invalid frame name. Please use 'perifocal', 'perifocal_t0', 'bodycentric' or a custom Frame object.")
        self.frame = frame
        self.plot3d = plot3d
    
        fig = plt.figure(figsize=(16, 9))
        if plot3d:
            ax = fig.add_subplot(projection='3d')
        else:
            ax = fig.add_subplot()
        self.ax = ax

    @staticmethod
    def color_gradient(value, colors=[(1,0.5,0), (0,0,1)]):
        """
        Returns an interpolated color based on a number between 0 and 1.

        :param colors: List of RGB colors to interpolate (each color as tuple of 3 values 0-1) 
        :param value: Value between 0 and 1 to interpolate the color
        :return: Interpolated RGB color tuple
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

    def draw_body(self, body, r_vec_bc, t_clock_bc, label="Central Body"):        
        # Body surface
        u = np.linspace(0, 2 * np.pi, 36)
        v = np.linspace(0, np.pi, 36)
        sphere_vec_bc = np.array([r_vec_bc[0] + body.radius * np.outer(np.cos(u), np.sin(v)),
                               r_vec_bc[1] + body.radius * np.outer(np.sin(u), np.sin(v)),
                               r_vec_bc[2] + body.radius * np.outer(np.ones(np.size(u)), np.cos(v))])
        sphere_vec = self.frame.transform_bc_to_frame(sphere_vec_bc, t_clock_bc)
        if label == "Central Body":
            color = 'red'
            # Equator line
            phi = np.linspace(0, 2 * np.pi, 36)
            equator_vec_bc = np.array([r_vec_bc[0] + body.radius * np.cos(phi),
                                       r_vec_bc[1] + body.radius * np.sin(phi),
                                       r_vec_bc[2] * np.ones(np.size(phi))])
            equator_vec = self.frame.transform_bc_to_frame(equator_vec_bc, t_clock_bc)
            #plot equator line
            if self.plot3d:
                self.ax.plot(equator_vec[0], equator_vec[1], equator_vec[2], color='gray', label='Equator')
            else:
                self.ax.plot(equator_vec[0], equator_vec[1], color='gray', label='Equator')
        else:
            color = 'blue'
        
        #plot body surface
        if self.plot3d:
            self.ax.plot_surface(sphere_vec[0], sphere_vec[1], sphere_vec[2], color=color, alpha=0.4, label=label)   # type: ignore
        else:
            center_vec = self.frame.transform_bc_to_frame(np.array([0, 0, 0]), t_clock_bc)
            circle = plt.Circle((center_vec[0], center_vec[1]), body.radius, color=color, alpha=0.4, label=label, zorder=2) # type: ignore
            self.ax.add_patch(circle)

    def draw_orbit(self, orbit, t_clock_bc=0):
        scale = orbit.rp  # Scale for auxiliary vectors
        h_vec = self.frame.transform_bc_to_frame(orbit.h_vec_bc, t_clock_bc)
        e_vec = self.frame.transform_bc_to_frame(orbit.e_vec_bc, t_clock_bc)
        n_vec = self.frame.transform_bc_to_frame(orbit.n_vec_bc, t_clock_bc)
        
        ### Scaling vectors for plotting ###
        if orbit.h == 0:
            h_vec = [0, 0, 0]
        else:
            h_vec = h_vec/orbit.h * scale
        if orbit.n == 0:
            n_vec = [0, 0, 0]
        else:
            n_vec = n_vec/orbit.n * scale
        if orbit.e == 0:
            e_vec = [0, 0, 0]
        else:
            e_vec = e_vec/orbit.e * scale

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
        
        r_vecs_bc = np.array([orbit.state_vectors_at_theta(theta, "bodycentric")[0] for theta in theta_list]).T
        r_vecs = self.frame.transform_bc_to_frame(r_vecs_bc, t_clock_bc)
        
        # Plot orbit
        center_vec = self.frame.transform_bc_to_frame(np.array([0, 0, 0]), t_clock_bc)
        if self.plot3d:
            self.ax.plot(r_vecs[0], r_vecs[1], r_vecs[2], 'b-', label='Orbit', linewidth=2)
            self.ax.quiver(center_vec[0], center_vec[1], center_vec[2], h_vec[0], h_vec[1], h_vec[2], color=(0,0,1), label='h - Angular Momentum')
            self.ax.quiver(center_vec[0], center_vec[1], center_vec[2], e_vec[0], e_vec[1], e_vec[2], color=(1,0,0), label='e - Eccentricity')
            self.ax.quiver(center_vec[0], center_vec[1], center_vec[2], n_vec[0], n_vec[1], n_vec[2], color=(0.8,0.8,0), label='n - Ascending Node')
        else:
            self.ax.plot(r_vecs[0], r_vecs[1], 'b-', label='Orbit', linewidth=2, zorder=1)
            self.ax.quiver(center_vec[0], center_vec[1], h_vec[0], h_vec[1], color='g', label='h', scale_units='xy', scale=1, width=0.003)
            self.ax.quiver(center_vec[0], center_vec[1], e_vec[0], e_vec[1], color='r', label='e', scale_units='xy', scale=1, width=0.003)
            self.ax.quiver(center_vec[0], center_vec[1], n_vec[0], n_vec[1], color='y', label='n', scale_units='xy', scale=1, width=0.003)

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
                    self.ax.plot(r_vec[0], r_vec[1], r_vec[2], 'o', color=color, label=label)
                else:
                    self.ax.plot(r_vec[0], r_vec[1], 'o', color=color, label=label, zorder=2)
    
    def draw_positions(self, orbit, velocities=True, groundtrack=False):
        if orbit.orbital_positions:
            for i, position in enumerate(orbit.orbital_positions):
                label = position.get('name', f'Position {i+1}')
                r_vec_bc, v_vec_bc = orbit.state_vectors_at_theta(position['theta'], "bodycentric")
                r_vec = self.frame.transform_bc_to_frame(r_vec_bc, position['t_clock'])
                v_vec = self.frame.transform_bc_to_frame(v_vec_bc, position['t_clock'])
                label += f' (t={position['t_clock']:.2f}s)'
                if position['t_clock'] < orbit.t0_clock:
                    t_min = min(position['t_clock'] for position in orbit.orbital_positions if position['t_clock'] is not None)
                    t_norm = (position['t_clock'] - t_min) / (orbit.t0_clock - t_min)
                    color = self.color_gradient(t_norm, [(0.5, 0, 0), (1, 0.5, 0)])
                else:
                    t_max = max(position['t_clock'] for position in orbit.orbital_positions if position['t_clock'] is not None)
                    t_norm = (position['t_clock'] - orbit.t0_clock) / (t_max - orbit.t0_clock)
                    color = self.color_gradient(t_norm, [(1, 0.5, 0), (0, 0.6, 0.6)])
                # Plot a marker at the position
                if self.plot3d:
                    self.ax.plot(r_vec[0], r_vec[1], r_vec[2], 'o', color=color, label=label)
                else:
                    self.ax.plot(r_vec[0], r_vec[1], 'o', color=color, label=label, zorder=2)
                
                # Plot velocity vector
                if velocities:
                    v_vec = v_vec * orbit.get_p()/20  # Velocity vector scale
                    if groundtrack:
                        v_r = np.dot(v_vec, r_vec/np.linalg.norm(r_vec))
                        v_vec = v_vec - v_r * r_vec/np.linalg.norm(r_vec)
                    if self.plot3d:
                        self.ax.quiver(r_vec[0], r_vec[1], r_vec[2], v_vec[0], v_vec[1], v_vec[2], color=color)
                    else:
                        self.ax.quiver(r_vec[0], r_vec[1], v_vec[0], v_vec[1], color=color, scale_units='xy', scale=1, width=0.003, zorder=5)

    """
    def draw_trajectory():
            if self.trajectory_points:
                annotations = []  # List to store all annotations
                # Trajectory points for annotations
                for point in self.trajectory_points:
                    r_vec = self.r_vec_at_theta(point['theta'], frame=frame)
                    if point['t_clock'] is not None:
                        t_clock = point['t_clock']
                        t_norm = (t_clock - self.t0_clock) / (self.t1_clock - self.t0_clock)
                        color = self.color_gradient(t_norm)
                        if groundtrack:
                            r_vec = r_vec / np.linalg.norm(r_vec) * self.body1radius
                        if plot3d:
                            point_plot = ax.plot(r_vec[0], r_vec[1], r_vec[2], 'o', color=color,
                                                picker=True, pickradius=5)[0]
                        else:
                            point_plot = plt.plot(r_vec[0], r_vec[1], 'o', color=color,
                                                picker=True, pickradius=5, zorder=4)[0]
                            # Create annotation
                            annotation = plt.annotate(f't = {t_clock:.2f}s',
                                xy=(r_vec[0], r_vec[1]), xytext=(10, 10),
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
                
                plt.gcf().canvas.mpl_connect('motion_notify_event', hover)
                
                # Generate high resolution points for the continuous line
                times = np.linspace(self.t0_clock, self.t1_clock, int((self.t1_clock - self.t0_clock)/60))
                r_vecs = []
                colors = []

                for t in times:
                    r_vec = self.state_vectors_at_t_clock(t, frame=frame)[0]
                    if groundtrack:
                        r_vec = r_vec / np.linalg.norm(r_vec) * self.body1radius
                    t_norm = (t - self.t0_clock) / (self.t1_clock - self.t0_clock)
                    r_vecs.append(r_vec)
                    colors.append(self.color_gradient(t_norm))

                # Convert to numpy array for better performance
                r_vecs = np.array(r_vecs)
                for i in range(len(r_vecs)-1):
                    if plot3d:
                        ax.plot(r_vecs[i:i+2,0], r_vecs[i:i+2,1], r_vecs[i:i+2,2], 
                                '-', color=colors[i], alpha=0.6, linewidth=2)
                    else:
                        plt.plot(r_vecs[i:i+2,0], r_vecs[i:i+2,1], 
                                    '-', color=colors[i], alpha=0.6, linewidth=2, zorder=3)
    """

    def plot_orbit(self, orbit, points=True, positions=True, velocities=True, groundtrack=False):
        """
        Plots the orbit
        """
        if not isinstance(orbit, Orbit):
            raise ValueError("O parâmetro 'orbit' deve ser uma instância da classe Orbit")
        
        if self.frame.name == "perifocal":
            self.frame.Omega0_bc_frame = orbit.Omega0
            self.frame.Omega_dot_bc_frame = orbit.Omega_dot 
            self.frame.omega0_bc_frame = orbit.omega0
            self.frame.omega_dot_bc_frame = orbit.omega_dot
            self.frame.i0_bc_frame = orbit.i
            self.frame.t_clock_bc_frame = orbit.t0_clock

        if self.frame.name == "perifocal_t0":
            self.frame.Omega0_bc_frame = orbit.Omega0
            self.frame.Omega_dot_bc_frame = 0
            self.frame.omega0_bc_frame = orbit.omega0
            self.frame.omega_dot_bc_frame = 0
            self.frame.i0_bc_frame = orbit.i
            self.frame.t_clock_bc_frame = orbit.t0_clock

        self.draw_body(orbit.main_body, r_vec_bc=np.array([0, 0, 0]), t_clock_bc=orbit.t0_clock, label="Central Body")
        self.draw_orbit(orbit, t_clock_bc=orbit.t0_clock)
        self.draw_points(orbit)
        self.draw_positions(orbit, velocities=velocities, groundtrack=groundtrack)

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

