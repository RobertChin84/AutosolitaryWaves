"""
This function will solve a nonlinear diffusion equation using MacCormack method
"""
import numpy as np

class Solver(object):
    """
    pass
    """
    def __init__(self, t, l, n=10, k=100):
        """
        :param dt: temporal step
        :param dx: spatial step
        :param model: the discretised pde
        """
        self.num_of_time_steps = n # number of time steps
        self.num_of_spatial_steps = k # number of spatial steps
        self.l = l
        self.t = t
        self.u = np.zeros([self.num_of_time_steps, self.num_of_spatial_steps])
        self.intialise_u()

    @property
    def dx(self):
        return self.l/float(self.num_of_spatial_steps)

    @property
    def dt(self):
        return self.t/float(self.num_of_time_steps)

    def intialise_u(self):

        x = np.linspace(0, self.l, self.num_of_spatial_steps)
        self.u[0][:] = np.sin(2.0*np.pi*x/self.l)

    def run(self):
        a = 1.0
        for n in range(self.num_of_time_steps):
            # compute predictor step
            u_pred = self.u[n][:-2] - a * (self.dt/self.dx) * (self.u[n][1:-1] - self.u[n][:-2])
            # corrector step
            self.u[n+1][:-2] = 0.5 * (self.u[n][:-2] + u_pred)



if __name__ == '__main__':
    solver = Solver(t=10.0, l=100.0)
    solver.run()
