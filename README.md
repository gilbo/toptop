This code was written for and is being made available to accompany the SIGGRAPH 2013 paper "Putting Holes in Holey Geometry: Topology Change for Arbitrary Surfaces"

For updates to this code, please check my website www.gilbertbernstein.com and/or my github www.github.com/gilbo .  This code will be located at www.github.com/gilbo/toptop for the forseeable future.

Gilbert Bernstein can be contacted at gilbert@gilbertbernstein.com


DEPENDENCIES
------------

This code depends on C++11 features.  It has been tested on Clang v4.0/LLVM 3.1.  A version which can be compiled on Visual Studio may be made available in the future.  Please inform me if you have prepared such a version yourself, and I can add it to the repository.

This code depends on the following external libraries:

QT (4.8.X and above)
GMP

If you're on Mac, both of these can be easily installed using the 'Homebrew' package manager.


INSTALLATION
------------

Once you have the necessary dependencies, edit 'makeConstants' and set the variable to the location of your GMP installation.

type 'make' in the root directory to build.


GIVE IT A SPIN
--------------

Run everything from the project root.  Try

./bin/modeler

to try out the interactive modeling tool.

look in ./resources/demos/README to learn how to run the code non-interactively and inspect parity fields


OMG -- WARNING
--------------

This research prototype (!!!) probably has bugs.  Please read the CRAPL ( http://matt.might.net/articles/crapl/ ) before complaining!

