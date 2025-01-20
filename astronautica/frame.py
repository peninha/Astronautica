from .body import Body
import numpy as np

class Frame:
    """
    Class representing a reference frame relative to the bodycentric frame.
    """
    def __init__(self, name, r0_vec_bc_frame=np.zeros(3), v_vec_bc_frame=np.zeros(3),
                 Omega0_bc_frame=0.0, Omega_dot_bc_frame=0.0,
                 omega0_bc_frame=0.0, omega_dot_bc_frame=0.0,
                 i0_bc_frame=0.0, i_dot_bc_frame=0.0,
                 t_clock_bc_frame=0.0):
        """
        Initializes a custom reference frame.
        
        Note:
            For standard frames, use the class methods:
            - Frame.bodycentric()
            - Frame.perifocal()
            - Frame.perifocal_t0()

        Parameters:
            name (str): The name of the frame.
            r0_vec_bc_frame (array): The origin position vector of the frame in the bodycentric frame at t0_clock_bc.
            v_vec_bc_frame (array): The origin velocity vector of the frame in the bodycentric. v_vec is considered constant.
            Omega0_bc_frame (float): The rotation rate of the frame in the bodycentric frame at t0_clock_bc.
            Omega_dot_bc_frame (float): The rotation rate derivative of the frame in the bodycentric frame. Omega_dot is considered constant.
            omega0_bc_frame (float): The rotation rate of the frame in the bodycentric frame at t0_clock_bc.
            omega_dot_bc_frame (float): The rotation rate derivative of the frame in the bodycentric frame. omega_dot is considered constant.
            i0_bc_frame (float): The inclination of the frame in the bodycentric frame at t0_clock_bc.
            i_dot_bc_frame (float): The inclination derivative of the frame in the bodycentric frame. i_dot is considered constant.
            t_clock_bc_frame (float): The time of the frame in the bodycentric frame.
        """
        self.name = name
        self.r0_vec_bc_frame = r0_vec_bc_frame
        self.v_vec_bc_frame = v_vec_bc_frame
        self.Omega0_bc_frame = Omega0_bc_frame
        self.Omega_dot_bc_frame = Omega_dot_bc_frame
        self.omega0_bc_frame = omega0_bc_frame
        self.omega_dot_bc_frame = omega_dot_bc_frame
        self.i0_bc_frame = i0_bc_frame
        self.i_dot_bc_frame = i_dot_bc_frame
        self.t_clock_bc_frame = t_clock_bc_frame

    @classmethod
    def bodycentric(cls):
        """
        Creates a bodycentric frame with zeros, because the bodycentric frame is the reference frame.
        """
        return cls(
            name="bodycentric",
            r0_vec_bc_frame=np.zeros(3),
            v_vec_bc_frame=np.zeros(3),
            Omega0_bc_frame=0,
            Omega_dot_bc_frame=0,
            omega0_bc_frame=0,
            omega_dot_bc_frame=0,
            i0_bc_frame=0,
            i_dot_bc_frame=0,
            t_clock_bc_frame=0
        )

    @classmethod
    def rotating_bodycentric(cls, main_body: Body, Omega0=0.0, t0_clock=0.0):
        """
        Creates a rotating bodycentric frame, from a Body object.
        """
        if not isinstance(main_body, Body):
            raise ValueError("main_body must be an instance of Body")
        
        return cls(
            name="rotating_bodycentric",
            r0_vec_bc_frame=np.zeros(3),
            v_vec_bc_frame=np.zeros(3),
            Omega0_bc_frame=Omega0,
            Omega_dot_bc_frame=main_body.rotation_speed,
            omega0_bc_frame=0,
            omega_dot_bc_frame=0,
            i0_bc_frame=0,
            i_dot_bc_frame=0,
            t_clock_bc_frame=t0_clock
        )

    ########### FRAMES TRANSFORMATION ###########
    @staticmethod
    def Q_from_Euler_angles(alpha, beta=0, gamma=0, pattern="classical"):
        """
        Calculates the direction cosine matrix from the Euler angles.

        """
        alpha = np.radians(alpha)
        beta = np.radians(beta)
        gamma = np.radians(gamma)
        if pattern == "1" or pattern == "x":
            Q = np.array([[1, 0, 0],
                          [0, np.cos(alpha), np.sin(alpha)],
                          [0, -np.sin(alpha), np.cos(alpha)]])
        elif pattern == "2" or pattern == "y":
            Q = np.array([[np.cos(alpha), 0, -np.sin(alpha)],
                          [0, 1, 0],
                          [np.sin(alpha), 0, np.cos(alpha)]])
        elif pattern == "3" or pattern == "z":
            Q = np.array([[np.cos(alpha), np.sin(alpha), 0],
                          [-np.sin(alpha), np.cos(alpha), 0],
                          [0, 0, 1]])
        elif pattern == "313" or pattern == "classical":
            Q = np.array([[np.cos(alpha)*np.cos(gamma) - np.sin(alpha)*np.cos(beta)*np.sin(gamma), np.sin(alpha)*np.cos(gamma) + np.cos(alpha)*np.cos(beta)*np.sin(gamma), np.sin(beta)*np.sin(gamma)],
                          [-np.cos(alpha)*np.sin(gamma) - np.sin(alpha)*np.cos(beta)*np.cos(gamma), -np.sin(alpha)*np.sin(gamma) + np.cos(alpha)*np.cos(beta)*np.cos(gamma), np.sin(beta)*np.cos(gamma)],
                          [np.sin(alpha)*np.sin(beta), -np.cos(alpha)*np.sin(beta), np.cos(beta)]])
        return Q

    def get_Q_bc_frame(self, t_clock_bc):
        """
        Returns the direction cosine matrix from the perifocal frame to the bodycentric frame at time t_clock.
        """
        return self.Q_from_Euler_angles(self.Omega0_bc_frame + self.Omega_dot_bc_frame * t_clock_bc,
                                        self.i0_bc_frame + self.i_dot_bc_frame * t_clock_bc,
                                        self.omega0_bc_frame + self.omega_dot_bc_frame * t_clock_bc,
                                        pattern="classical")

    def transform_bc_to_frame(self, r_vec_bc, t_clock_bc=0):
        """
        Transforms vectors from the bodycentric frame to the actual frame at time t_clock.
        
        Parameters:
            r_vec_bc: Vector in bodycentric frame with shape (3,), (3,n) or (3,n,n)
            t_clock_bc: Time
        Returns:
            Transformed vector in frame - same shape as input
        """
        r_vec_bc = np.array(r_vec_bc)
        is_single_vector = r_vec_bc.shape == (3,)
        
        # Calculate displacement
        displacement = self.r0_vec_bc_frame - self.v_vec_bc_frame * t_clock_bc
        
        Q = self.get_Q_bc_frame(t_clock_bc)
        # Special case for single vector (3,)
        if is_single_vector:
            r_vec_frame_prime = r_vec_bc - displacement
            r_vec_frame = Q @ r_vec_frame_prime  # simple matrix multiplication for single vector
            return r_vec_frame
        
        # For larger arrays (3,n) or (3,n,n)
        r_vec_frame_prime = r_vec_bc - displacement.reshape(3, 1, 1) if r_vec_bc.ndim == 3 else r_vec_bc - displacement.reshape(3, 1)
        r_vec_frame = np.einsum('ij,j...->i...', Q, r_vec_frame_prime)
        
        return r_vec_frame
    
    ########### CONVERSION ###########
    @staticmethod
    def convert_cartesian_to_ra_dec(r_vec_xyz):
        """
        Calculates the right ascension and declination from the Cartesian position vector.
        """
        r = np.linalg.norm(r_vec_xyz)
        ra = np.degrees(np.arctan2(r_vec_xyz[1], r_vec_xyz[0]))
        dec = np.degrees(np.arcsin(r_vec_xyz[2]/r))
        return ra, dec