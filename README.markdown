Library name: libmultigrid
Author: Jesse Robertson, Australian National University
Date: 17 October 2011
URL: http://github.com/jess-robertson/multigrid
Email: my name with domain anu.edu.au

README:

This is a multigrid solver for solving elliptic PDEs using finite differences on a rectangular grid. It uses red-black updating with a user-specified smoother and fully weighted restriction/bilinear interpolation for solution transfer between grids. These are implemented as methods of a solver class (mgrid::LinearMultigrid) so that you don't have to specify too much to get the solver running.

Unfortunately there isn't much in the way of documentation other than this README but you can have a look at the two examples and get the general idea of what's going on.

I'm releasing this under the the Community Research and Academic Programming License (CRAPL: see http://matt.might.net/articles/crapl/ for details). Hopefully things'll work as advertised but if you have any problems then drop me a line and we'll see if we can fix them.



Installing
----------

You'll need to install cmake, and the boost and blitz libraries to compile the multigrid library. On a Mac the easiest thing to do is use Homebrew (), so that this can be done in one hit with the command `brew install cmake boost blitz`.

Once you've installed these dependencies, you can go to the root folder (where this README is located) and run 

	cmake . && make install 

...this should put the library and headers under /usr/local. Feel free to modify the install directory in the CMakeLists.txt file if you want it to go somewhere else.

When building your own script, you can also use the CMake files that are included in the two example directories. These should find the multigrid library and link it into your own executable.



Specifying your own elliptic PDEs
---------------------------------

You can use this library by subclassing from the mgrid::LinearMultigrid class. You need to specify a differential operator, and a smoothing operator which will be specific to your given PDE. These functions are given a multigrid level, and i and j indexes which they can use to calculate the smoothing step. To make specifying these things easier, the solution attribute also implements finite differences of the specified order in the form:

	solution[level].dxx(i, j)
	solution[level].dzz(i, j)

...which would give you a finite difference approximation to the second derivative in the x and z directions respectively. These finite differences also take the location of the point into account, so if you're near a boundary they will automatically use forward or backward differences as required.

For example, to solve the Poisson equation for a solution u_{h} on a grid with spacing h = (hx, hz), with a source term f_{h}, we have the following differential operator:

L(u_{h}) = d^2u/dx^2 + d^2u/dz^2

and the following Gauss-Seidel update step as a smoother:

S(u_{h}) = u_{h}(i+1, j)((u_{h}(i+1, j) + u_{h}(i-1, j))*/(hx^2)
        + (u_{h}(i, j+1) + u_{h}(i, j-1))/(hz^2) 
        - f_{h}(i, j))/(2*(hx + hz))

These are implemented in the Poisson example as

	inline double differential_operator(mgrid::Level level, int i, int j) {
	    return solution[level].dxx(i, j) + solution[level].dzz(i, j);
	};
	inline void Poisson::relaxation_updater(mgrid::Level level, int i, int j) { 
	    const double hx = solution[level].spacing(0);
	    const double hz = solution[level].spacing(1);
	    const double xxfactor = 1/(hx*hx);
	    const double zzfactor = 1/(hz*hz);
	    solution[level](i, j) = 
	        ((solution[level](i+1, j) + solution[level](i-1, j))*xxfactor
	        + (solution[level](i, j+1) + solution[level](i, j-1))*zzfactor 
	        - source[level](i, j))/(2*(xxfactor + zzfactor));
	};

As you can see, you can get the spacing size from the current level of the solution. The solution and source arrays are stored as attributes containing vectors of grids, arranged from fine to coarse, and are called `solution` and `source` respectively. You shouldn't need to bother with the order too much, just pass the supplied level through.

The source term is automatically moved between grids using the same restriction/prolongation operators as used for the solutions. To make it wasier to specify initial conditions and a source term, the mgrid::LinearMultigrid class also contains finestLevel and coarsestLevel attributes, with the index of the finest and coarsest grid level respectively.



Specifying boundary conditions
------------------------------

Boundary conditions can be set in the initialisation routine of your solver class which has been subclassed from mgrid::LinearMultigrid. You set them for the solution attribute.

Here's how it's done for the Poisson example, with two homogeneous Dirichlet and two Neumann conditions

	// Set boundary conditions for velocity array 
	solution.boundaryConditions.set(mgrid::leftBoundary,   mgrid::zeroNeumannCondition);
	solution.boundaryConditions.set(mgrid::rightBoundary,  mgrid::zeroDirichletCondition);
	solution.boundaryConditions.set(mgrid::topBoundary,    mgrid::zeroNeumannCondition);    
	solution.boundaryConditions.set(mgrid::bottomBoundary, mgrid::zeroDirichletCondition); 

If your conditions are not homogeneous, you can specify a value using the BoundaryPoint class. For example, the conditions above were specified as:

	const mgrid::BoundaryPoint zeroDirichletCondition = { mgrid::dirichlet, 0.0 };

A one-dimensional blitz array of BoundaryPoints makes an instance of mgrid::Boundary - you can also pass one of these to mgrid::BoundaryConditions::set to set the boundary conditions.



Specifying solver settings
--------------------------

Most of the settings are fairly self-explanatory:

	struct Settings { 
	    double aspectRatio; 			# Aspect ratio of the grids
	    int numberOfGrids; 				# Total number of grids
	    int minimumResolution;			# Resolution on the coarsest grid
	    double residualTolerance;		# Stopping tolerance for solver
	    int maximumIterations;			# Maximum allowable iterations for the solver
	    CycleType mgCycleType;			# Either mgrid::wCycle or mgrid::vCycle
	    unsigned long preMGRelaxIter;	# Number of relaxation iterations on way down
	    unsigned long postMGRelaxIter;	# Number of relaxation iterations on way back up
	};

...although a couple need a bit more explanation. The library will work out all the grid sizes based on your minimum grid size and the total number of grids by simply doubling the resolution in both the x and z direction with each grid. 

You specify the minimum grid size as the minimum resolution on the smallest side of the grid - the library will adjust the x and z spacings to have as close to the same resolution in both directions as possible using the aspect ratio setting.



Running the solver
------------------

You can run the multigrid solver using the mgrid::LinearMultigrid::multigrid method. There's an optional mgrid::LinearMultigrid::solve method that you can do more complicated stuff with. For example the viscoplastic channel flow example requires a linear elliptic PDE to be solved at each step, and the source term updated from the last solution. The solve method deals with this recalculation of the source term and then calls the multigrid method.
