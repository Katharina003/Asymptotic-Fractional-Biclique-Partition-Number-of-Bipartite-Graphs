# Asmyptotic Fractional Biclique Partition Number of Bipartite Graphs

This is the supplementary material for the master's thesis "Asymptotic Fractional Biclique Partition Number of Bipartite Graphs" submitted by Lea Katharina Herrmann.

## Contents
You find all certificates for $k\le 7$ in resultK1.txt to resultK7.txt.
You can certify the files for yourself using the python script in tools.
The output of the python script lies in certificates.txt.

All other files are part of the code of my program.

## Building

The project requires CMake 3.16 or newer, a C++17 compiler, Gurobi with its C++ headers and libraries, and the Boost headers. A valid Gurobi license is required to run the program. `GUROBI_HOME` must point to the Gurobi platform directory containing `include` and `lib`; `BOOST_ROOT` must point to the directory containing `boost/`.

Run the following commands from the project directory, replacing the dependency paths with your installation paths.

### Linux (CMake)

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release \
    -DGUROBI_HOME="/opt/gurobi1300/linux64" \
    -DBOOST_ROOT="/path/to/boost"
cmake --build build/release --parallel
```

The executable is `build/release/kronecker`. To run it from the project directory type:

```sh
./build/release/kronecker
```

To run the verifier for all seven certificates from the project directory, type (Python 3 required):

```sh
python3 tools/verify_primal.py resultK1.txt resultK2.txt resultK3.txt resultK4.txt resultK5.txt resultK6.txt resultK7.txt
```

To verify a single certificate, pass only its filename. The results are printed to the terminal and written to `certificates.txt`, replacing its previous contents.
