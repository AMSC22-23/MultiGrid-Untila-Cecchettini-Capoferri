# Implementation of Multigrid method

## Method description

The Multigrid method is an iterative numerical technique used to solve discretized systems of equations, particularly those arising from partial differential equations (PDEs). It's a powerful method for accelerating the convergence of iterative solvers, especially for large-scale problems.


### Compilation

* Run "make" to compile the entire codebase. This command will build the executable from the source files.
* Execute "make remake" if the code is already compiled, and you wish to recompile it.

### Parallel Compilation

For compiling a parallel version of the program, use "make parallel". This command will compile the code with specific flags, OpenMP flag, to enable parallel processing.

### Execution

* After successful compilation, execute the program by running ./Multigrid.
* To explore available options and input parameters, utilize the --help command. This will provide comprehensive information on the inputs that can be provided to the program.

### Execution Parameters:

When executing ./Multigrid, you can provide the following parameters:
* '-a' (or --alpha):
    Specifies differential constant.
* '-w' (or --width):
    Specifies the width of the squared domain.
* '-n' (or --intervals):
    Sets the discretization number used to discretize the domain.
* '-ml' (or --multigrid):
    Specifies levels of multigrid
* 'test':
    Specifies test for the exe

### Example
```
# Example 1: Compile the code if not already compiled
~$ make

# Example 2: Execute the program
~$ ./Multigrid -a 10.0 -w 5.0 -n 100 -ml 5 test 2

```


### Test

To conduct a test using an Ubuntu machine, navigate to the "/test/test.ipynb" directory. Modify the desired parameters, including the width and length of the squared domain, discretization number, level of multigrid coarseness and the number of test functions. Choose from the following available test functions:

* Test 1:

$$f(x, y) = -5.0 \cdot e^{x} \cdot e^{-2.0 \cdot y}$$
$$g(x, y) = e^{x} \cdot e^{-2.0 \cdot y}$$

* Test 2:

$$ f(x, y) = \sin(k \cdot \sqrt{x^2 + y^2})$$
$$g(x, y) = -k \left( \frac{\cos(kr)}{r} - k\sin(kr) \right)$$


* Test 3:

$$f(x,y) = 1$$
$$g(x,y) = 0$$

Upon completion, refer to the "result" file to examine the matrix solution. Additionally, a 3D graphic representation of the solution is available for visualization. 
Note that f is the function defined within the domain, while g specifies the function on the boundary.

### Example of initialization in ipynb file

```
finestGridN = 25 # number for discretization
levels = 5 # multigrid level
alpha = 1.
width = 10. # dimension of squared domain
testCase = 3 # number of available test

```