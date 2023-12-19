# Implementation of Multigrid method

## Method description

The Multigrid method is an iterative numerical technique used to solve discretized systems of equations, particularly those arising from partial differential equations (PDEs). It's a powerful method for accelerating the convergence of iterative solvers, especially for large-scale problems.

TODO: Adding other things

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
* '-a' (or --length):
    Specifies the length of the squared domain.
* '-w' (or --width):
    Specifies the width of the squared domain.
* '-n' (or --intervals):
    Sets the discretization number used to discretize the domain.

### Example
```
# Example 1: Compile the code if not already compiled
~$ make

# Example 2: Execute the program
~$ ./Multigrid -a 10.0 -w 5.0 -n 100

```