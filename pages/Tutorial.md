# Tutorial

## Using the executable

If all you want is to run one simulation with custom parameters, you can use the executable that is generated when you build the project. The executable is located in the build folder.

To build the executable you can use the following commands:

```bash
mkdir build
cd build
cmake ..
make WaveSimCExec
```

And then run it with the following command:

```bash
./WaveSimCExec
```

The executable will ask you for the parameters of the simulation. A complete description of the parameters can be found in the the help menu of the executable. To see the help menu you can run the executable with the -h flag:

```bash
./WaveSimCExec -h
```

## Creating a wave simulation and solving it

The source code of this tutorial can be found in src/examples/wave_solver_with_animation.cpp

The first part of creating the wave simulation is to define the constants for the simulation. These constants are the number of grid points in the x and z directions, the number of time steps, the differentiation values, the domain, the source parameters, and the source location.

It is important to pay attention to the fact that the plane is XZ not XY due to geophysical conventions.

```cpp
// Define the constants for the simulation

// Number of x and z grid points
int nx = 100;
int nz = 100;
// Number of time steps
int nt = 1000;

// Differentiation values
double dx = 0.01;
double dz = 0.01;
double dt = 0.001;

// Define the domain
double xmin = 0.0;
double xmax = nx * dx;
double zmin = 0.0;
double zmax = nz * dz;
double tmin = 0.0;
double tmax = nt * dt;

// Define the source parameters
double f_M = 10.0;
double amp = 1e0;
double shift = 0.1;

// Source location
int source_is = 50;
int source_js = 50;
```

The next step is to create the source object and the velocity profile

```cpp
// Create the source
boost::multi_array<double, 3> f = waveSimCore::ricker(source_is, source_js, f_M, amp, shift, tmin, tmax, nt, nx, nz);

// Create the velocity profile
double r = 150.0;
boost::multi_array<double, 2> vel = waveSimCore::get_profile(xmin, xmax, zmin, zmax, nx, nz, r);
```

Then we can proceed to solve the wave equation using the wave solver.

```cpp
// Solve the wave equation
boost::multi_array<double, 3> u = waveSimCore::wave_solver(vel, dt, dx, dz, nt, nx, nz, f);
```

u is the multi_array that contains the result of the simulation. It has the shape (nt, nx, nz). That means that the first index is the time step, the second index is the x grid point, and the third index is the z grid point.

You can access the result of a specific time step by simply using:

```cpp
boost::multi_array<double, 2> u_at_time_20 = u[20];
```

## Obtaining the results

To see the results you have two main options:

- Use the wavePlotter::Plotter class to plot the results and/or create an animation.
- Export all the frames to a series of csv files for processing in another program.

### Plotting the results using the wavePlotter::Plotter class

You should do some light processing to conver the domain from a boost::multi_array to a matplot::vector_2d. This is done using the np::convert_to_matplot function.

After that you can convert the result of each frame to a matplot::vector_2d and plot it using the wavePlotter::Plotter class.

```cpp
// Define the number of different levels for the contour plot
int num_levels = 100;
// Create the levels for the contour plot based on the min and max values of u
double min_u = np::min(u);
double max_u = np::max(u);
std::vector<double> levels = matplot::linspace(min_u, max_u, num_levels);

// Create the x and z axis for the contour plot and convert them to matplot format
boost::multi_array<double, 1> x = np::linspace(xmin, xmax, nx);
boost::multi_array<double, 1> z = np::linspace(zmin, zmax, nz);
const boost::multi_array<double, 1> axis[2] = {x, z};
std::vector<boost::multi_array<double, 2>> XcZ = np::meshgrid(axis, false, np::xy);

matplot::vector_2d Xp = np::convert_to_matplot(XcZ[0]);
matplot::vector_2d Zp = np::convert_to_matplot(XcZ[1]);
```

From this point on you can plot it as you wish, for example, if you want a filled contour plot you can do:

```cpp
matplot::vector_2d Up = np::convert_to_matplot(this->u[frame_index]);
matplot::contourf(this->Xp, this->Zp, Up, this->levels);
matplot::show();
```

Another option is to pass the data to the wavePlotter:Plotter class and use the plot function to render all frames

```cpp
// Create the plotter object and animate the results
wavePlotter::Plotter my_plotter(u, Xp, Zp, num_levels, nt);

// If you want to render a specific frame, use this:
// my_plotter.renderFrame(int frame_index);

// Renders the entire animation from start_frame to end_frame
int start_frame = 20;
int end_frame = nt - 1;
int fps = 30;
my_plotter.animate("example-wave.mp4", start_frame, end_frame, fps);
```

The animation will be saved in . and the frames will be saved to ./output

### Exporting the results to csv files

You can export the results to a series of individual csv files using the wavePlotter::Plotter class.

```cpp
// Create the plotter object and animate the results
wavePlotter::Plotter my_plotter(u, Xp, Zp, num_levels, nt);

// If you want to export a specific frame, use this:
// my_plotter.exportFrame(int frame_index);

// Exports the entire simulation from start_frame to end_frame
int start_frame = 20;
int end_frame = nt - 1;
my_plotter.exportAllFrames(start_frame, end_frame);
```

The frames will be saved to ./output

Each frame is saved as a csv file with the following format:

```csv
data at 0_0, data at 0_1, data at 0_2, ..., data at 0_nz
data at 1_0, data at 1_1, data at 1_2, ..., data at 1_nz
data at 2_0, data at 2_1, data at 2_2, ..., data at 2_nz
.
.
.
```

If you want to import the data into a python program you can use the following code:

```python
frames = []
for i in range(999):
    filename = "output/frame_" + f'{i:08}' +".csv"
    frames.append(pd.read_csv(filename))
```

Please refer to src/examples/PythonLoadingExample.ipynb for a complete example on how to load the data into a python program and render it using matplotlib (of course you can use any other library to render the data).
