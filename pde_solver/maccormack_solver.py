"""
This function will solve a nonlinear diffusion equation using MacCormack method
"""
import numpy as np
import matplotlib.pyplot as plt

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
        # set the initial condtions
        self.x = np.linspace(0, self.l, self.num_of_spatial_steps)
        self.u[0][:] = np.sin(2.0*np.pi*self.x/self.l)

    def boundary_conditions(self, n, condition_type='periodic'):
        if condition_type == 'periodic':
             self.u[n+1][-1] = self.u[n+1][0]

    def run(self):
        a = 10.0
        for n in range(self.num_of_time_steps-1):
            # compute predictor step

            u_pred = self.u[n][:-1] - a * (self.dt/self.dx) * (self.u[n][1:] - self.u[n][:-1])
            # corrector step
            for j in range(self.num_of_spatial_steps-1):
                self.u[n+1][j] = 0.5*(self.u[n][j] + u_pred[j]) - a * 0.5*(self.dt/self.dx)*(u_pred[j] - u_pred[j-1])
            self.boundary_conditions(n)
            plt.plot(self.x, self.u[n+1],label="time={}".format(round(n*self.dt,3)))
            plt.legend(loc=0)
            plt.show()



if __name__ == '__main__':
    solver = Solver(t=10.0, l=100.0)
    solver.run()
