import <vector>;

export FiniteDifferenceWaveSolvers;

std::vector<double> np_gradient(std::vector<double> f, std::vector<double> dx, std::vector<double> dz)
{
    std::vector<double> dfdx, dfdz;
    dfdx, dfdz = np.gradient(f, dx, dz) return dfdx, dfdz
}

std::vector<double> dfdx(std::vector<double> fx, std::vector<double> dx) : dfdx, _ = np.gradient(fx, dx, dz) return dfdx

                                                                                     def dfdz(fz, dz) : _,
                                                                                 dfdz = np.gradient(fz, dx, dz) return dfdz

                                                                                        def divergence(fx, fz, dx, dz) : dfdx,
                                                                                 _ = np.gradient(fx, dx, dz)
                                                                                         _,
                                                                                 dfdz = np.gradient(fz, dx, dz) return dfdx + dfdz

                                                                                        def laplacian(f, dx, dz) : dfdx,
                                                                                 dfdz = np.gradient(f, dx, dz)
                                                                                            d2fdx2,
                                                                                 d2fdxdz = np.gradient(dfdx, dx, dz)
                                                                                               d2fdzdx,
                                                                                 d2fdz2 = np.gradient(dfdz, dx, dz) return d2fdx2 + d2fdz2

                                                                                          def d2dt2(f, dt) : dfdt,
                                                                                 _, _ = np.gradient(f, dt, dt, dt)
                                                                                            d2fdt2,
                                                                                 _, _ = np.gradient(dfdt, dt, dt, dt) return d2fdt2

                                                                                        def dfdt(f, dt) : dfdt,
                                                                                 _, _ = np.gradient(f, dt, dt, dt) return dfdt

                                                                                        def timer(start, end) : hours,
                                                                                 rem = divmod(end - start, 3600)
                                                                                     minutes,
                                                                                 seconds = divmod(rem, 60)
                                                                                     print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds))
