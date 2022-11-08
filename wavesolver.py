class Domain(object):
    def __init__(self, nx, nz, nt, dx, dz, dt) -> None:
        self._nx = nx
        xmin = 0
        xmax = nx * dx  # (m)
        zmin = 0
        zmax = nz * dz  # (m)
        tmin = 0
        tmax = nt * dt  # (s)
        self._zmin = zmin
        self._zmax = zmax
        self._tmax = tmax
        self._x = 1

    @property
    def nx(self):
        return self._nx

    @nx.setter
    def nx(self, value):
        self._nx = value

    def mesh(self):
        x = np.linspace(self.xmin, self.xmax, self.nx)
        z = np.linspace(self.zmin, self.zmax, self.nz)
        X, Z = np.meshgrid(x, z, indexing='ij')
        return X, Z


class DifferentialOperator(object):
    def __init__(self, f):
        self._f = f
        self._dfdx = None
        self._dfdz = None

    def diff(self):
        dfdx, dfdz = np.gradient(fz, dx, dz)
        self._dfdx = dfdx
        self._dfdz = dfdz

    @property
    def dfdx(self):
        return self._dfdx

    @property
    def dfdz(self):
        return self._dfdz

    @property
    def laplacian(self):
        d2fdx2, d2fdxdz = np.gradient(self.dfdx, self.dx, self.dz)
        d2fdzdx, d2fdz2 = np.gradient(self.dfdz, self.dx, self.dz)
        return d2fdx2 + d2fdz2


class VelocityProfile(object):
    def __init__(self):
        1


class WaveSolver(object):
    def __init__(self, domain) -> None:
        self._domain = domain
        self._u = np.zeros((nt, nx, nz))
        self._u_xx = np.zeros((nx, nz))
        self._u_zz = np.zeros((nx, nz))
        self._q_1 = np.zeros((nx, nz))
        self._q_2 = np.zeros((nx, nz))

    @property
    def u(self):
        return self._u

    @property
    def q_1(self):
        return self._q_1

    def compute(self):
        C1 = 1 + dt*(sigma_1 + sigma_2)/2
        C2 = sigma_1 * sigma_2 * (dt**2) - 2
        C3 = 1 - dt*(sigma_1 + sigma_2)/2
        C4 = (dt*c)**2
        C5 = 1 + dt*sigma_1/2
        C6 = 1 + dt*sigma_2/2
        C7 = 1 - dt*sigma_1/2
        C8 = 1 - dt*sigma_2/2

        for n in range(1, nt-1):
            u_xx[1:-1, 1:-1] = u[n, :-2, 1:-1] - \
                2*u[n, 1:-1, 1:-1] + u[n, 2:, 1:-1]
            u_zz[1:-1, 1:-1] = u[n, 1:-1, :-2] - \
                2*u[n, 1:-1, 1:-1] + u[n, 1:-1, 2:]

            u[n+1] = (C4*(u_xx/(dx**2) + u_zz/(dz**2) - divergence(q_1*sigma_1, q_2*sigma_2, dx, dz) + sigma_2*dfdx(q_1, dx) + sigma_1*dfdz(q_2, dz) + f[n]) -
                      C2 * u[n] - C3 * u[n-1]) / C1

            q_1 = (dt*dfdx(u[n], dx) + C7*q_1) / C5
            q_2 = (dt*dfdz(u[n], dx) + C8*q_2) / C6

            # Enforce Dirichlet boundary condition
            u[n+1, :, 0] = np.zeros(nx)
            u[n+1, :, nz-1] = np.zeros(nx)
            u[n+1, 0, :] = np.zeros(nz)
            u[n+1, nx-1, :] = np.zeros(nz)

        return u


if __name__ == '__main__':
    domain = Domain(1, 2, 3, 0.1, 0.1, 0.1)
    X, Z = domain.mesh()
    op = DifferentialOperator(f)
    op.diff()
