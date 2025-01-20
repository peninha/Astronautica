from .orbit import Orbit

class Three_body_restricted(Orbit):
    """
    Class to calculate the position and velocity vectors of a body in a three-body restricted problem.
    Frame of reference: noninertial, rotating system, with origin at the center of mass and X axis towards m2.
    Note: Valid only for the case where the two main bodies (m1 and m2) are in circular orbit around each other.
    """

    ########### CONSTRUCTORS ###########
    def __init__(self, body1, body2, _from_classmethod=False):
        """
        Initializes the Three_body_restricted class.

        :param body1: Body object of the primary body
        :param body2: Body object of the secondary body
        """
        self.m1 = m1
        self.m2 = m2
        self.mu1 = self.G * m1
        self.mu2 = self.G * m2
        self.pi1 = m1/(m1 + m2)
        self.pi2 = m2/(m1 + m2)
        self.mu = self.mu1 + self.mu2
        self.r12 = r12
        self.body1_center = np.array([-self.pi2*self.r12, 0, 0])
        self.body2_center = np.array([(1-self.pi2)*self.r12, 0, 0])
        self.body1radius = body1radius
        self.body2radius = body2radius
        self.omega = np.sqrt(self.mu/self.r12**3)

    def lagrange_points(self, add_points=True):
        """
        Calculates the Lagrange points.
        """
        
        #### L4 and L5 ####
        x = self.r12/2 - self.pi2 * self.r12
        y = np.sqrt(3)/2 * self.r12
        z = 0
        L4 = np.array([x, y, z])
        L5 = np.array([x, -y, z])

        #### L1, L2 and L3 ####
        def f(csi, pi2):
            return (1 - pi2) * (csi + pi2)/np.abs(csi + pi2)**3 + pi2 * (csi + pi2 - 1)/np.abs(csi + pi2 -1)**3 - csi
        
        sign = 1 if self.m1 >= self.m2 else -1
        L1_csi0 = (self.pi1 - self.pi2)
        L2_csi0 = (1 + min(self.pi1, self.pi2)) * sign
        L3_csi0 = (1 + min(self.pi1, self.pi2)) * -sign
        L1_csi = fsolve(f, L1_csi0, args=(self.pi2))[0]
        L1 = np.array([L1_csi * self.r12, 0, 0])
        L2_csi = fsolve(f, L2_csi0, args=(self.pi2))[0]
        L2 = np.array([L2_csi * self.r12, 0, 0])
        L3_csi = fsolve(f, L3_csi0, args=(self.pi2))[0]
        L3 = np.array([L3_csi * self.r12, 0, 0])

        #print(L1_csi0, L1_csi)
        #print(L2_csi0, L2_csi)
        #print(L3_csi0, L3_csi)

        if add_points:
            self.add_point_in_orbital_plane(*self.convert_cartesian_to_polar(L1), name="L1")
            self.add_point_in_orbital_plane(*self.convert_cartesian_to_polar(L2), name="L2")
            self.add_point_in_orbital_plane(*self.convert_cartesian_to_polar(L3), name="L3")
            self.add_point_in_orbital_plane(*self.convert_cartesian_to_polar(L4), name="L4")
            self.add_point_in_orbital_plane(*self.convert_cartesian_to_polar(L5), name="L5")
        return L1, L2, L3, L4, L5

    def jacobi_constant(self, r_vec, v):
        """
        Calculates the Jacobi constant.
        """
        r1 = np.linalg.norm(np.array([r_vec[0]+self.pi2*self.r12, r_vec[1], r_vec[2]]))
        r2 = np.linalg.norm(np.array([r_vec[0]-self.pi1*self.r12, r_vec[1], r_vec[2]]))
        C = (v**2 - self.omega**2 * (r_vec[0]**2 + r_vec[1]**2) - 2*self.mu1/r1 - 2*self.mu2/r2)/2
        return C
    
    def v_for_C(self, r_vec, C):
        """
        Calculates the velocity vector from some specific Jacobi constant (energy).
        """
        r1 = np.linalg.norm(np.array([r_vec[0]+self.pi2*self.r12, r_vec[1], r_vec[2]]))
        r2 = np.linalg.norm(np.array([r_vec[0]-self.pi1*self.r12, r_vec[1], r_vec[2]]))
        v = np.sqrt(2*C + self.omega**2 * (r_vec[0]**2 + r_vec[1]**2) + 2*self.mu1/r1 + 2*self.mu2/r2)
        return v

    def trajectory(self, r_vec_0, v_vec_0, t_span, t_eval=None, method='RK45', add_trajectory_points=True):
        """
        Simulates the trajectory of the spacecraft in the rotating frame.
        """
        def f(t, y):
            y1 = y[0] # x
            y2 = y[1] # y
            y3 = y[2] # z

            y4 = y[3] # x'
            y5 = y[4] # y'
            y6 = y[5] # z'
            
            r1 = np.sqrt((y1 + self.pi2*self.r12)**2 + y2**2 + y3**2)
            r2 = np.sqrt((y1 - self.pi1*self.r12)**2 + y2**2 + y3**2)

            y1_dot = y4
            y2_dot = y5
            y3_dot = y6
            y4_dot = 2*self.omega*y5 + self.omega**2*y1 - self.mu1/r1**3 * (y1 + self.pi2*self.r12) - self.mu2/r2**3 * (y1 - self.pi1*self.r12)
            y5_dot = -2*self.omega*y4 + self.omega**2*y2 - self.mu1/r1**3 * y2 - self.mu2/r2**3 * y2
            y6_dot = - self.mu1/r1**3 * y3 - self.mu2/r2**3 * y3
            return np.concatenate(((y1_dot, y2_dot, y3_dot), (y4_dot, y5_dot, y6_dot)))
        
        sol = solve_ivp(f, t_span, np.concatenate((r_vec_0, v_vec_0)), t_eval=t_eval, method=method)
        if add_trajectory_points:
            # Converter cada ponto (x,y) para (r,theta)
            for point in sol.y[:2].T:  # Pegamos apenas x e y, ignorando z
                r, theta = self.convert_cartesian_to_polar(point)
                self.add_trajectory_points(r, theta)
        return sol

