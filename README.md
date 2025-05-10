# SREES - State Estimation for Electrical Systems

This is a cross-platform application built for the course **Structures and Schemes of Electrical Power Systems**. The project is designed to compare the Gauss–Seidel and Newton-Raphson Methods for estimating the state of the power grid.

## Overview

The application provides a graphical user interface (GUI) built using GTK, allowing users to perform state estimation of power grids using the two aforementioned methods. The project is cross-platform and can be compiled to run on various operating systems. 

Currently, the only provided binary is for Linux, named `App`.

## Features

- Compare Gauss–Seidel and Newton-Raphson Methods for power grid state estimation
- Graphical user interface (GTK-based)
- Cross-platform (currently only Linux binary available)
- Designed for educational and research purposes in power grid analysis

## Skills

This project involves the following technical skills:

- **C++ Programming**: Core language for application development.
- **GTK**: Used for creating the graphical user interface (GUI).
- **State Estimation Algorithms**: Implemented Gauss-Seidel and Newton-Raphson Methods for power grid analysis.
- **Cross-Platform Development**: Designed to work across multiple operating systems (currently Linux binary).
- **Software Engineering**: Knowledge of object-oriented design, modular programming, and code optimization.

## Project Structure

```

SREES-main/
├── App                     # Executable for Linux (currently provided)
├── App.cpp                 # Main application source code
├── App.h                   # Header file for main application
├── GSSolver.cpp            # Gauss-Seidel solver implementation
├── GSSolver.h              # Header for Gauss-Seidel solver
├── Line.cpp                # Line class implementation
├── Line.h                  # Header for Line class
├── NRSolver.cpp            # Newton-Raphson solver implementation
├── NRSolver.h              # Header for Newton-Raphson solver
├── Node.cpp                # Node class implementation
├── Node.h                  # Header for Node class
├── README.md               # This readme file
└── main.cpp                # Main entry point for the application

````

## Requirements

- **Operating System**: Linux (for provided binary)
- **GTK**: GTK+3 or later for GUI support
- **Compiler**: g++ or compatible C++ compiler

## Build Instructions

To build the application from source on Linux:

1. Clone the repository:
   ```bash
   git clone <repository_url>
   ````
2. Navigate to the project directory:
   ```bash
   cd SREES-main
   ```
3. Compile the code:
   ```bash
   g++ -o App main.cpp App.cpp GSSolver.cpp NRSolver.cpp Line.cpp Node.cpp `pkg-config --cflags --libs gtk+-3.0`
   ```
4. Run the application:
   ```bash
   ./App
   ```

## License

This project is for educational and demonstration purposes.

## Contributing

Feel free to fork this repository, open issues, and submit pull requests. Contributions are welcome!
