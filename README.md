## Synopsis

`iB4e-GTFold-parametrizer` determines the sensitivity of RNA folding to multibranch-loop-related parameters.
It relies on an implementation of Zucker's dynamic-programming algorithm for RNA secondary structure prediction, using adapted code from the `gtmfe` program from the [GTFold][gtfold] project.

## Installation

This project is under active development, so we recommend downloading it using Git.
To do this, run the following in your terminal:

```
git clone https://github.com/AMS-MRC-disc-math-bio/iB4e-GTfold-parametrizer.git
```

This will download the code and extract it into a directory called `iB4e-GTfold-parameterizer`.

### Dependencies
This project depends on the [CGAL][cgal] computational geometry library, several libraries from the [Boost][boost] project, and the [GMP][gmp] arbitrary-precision arithmetic library.
To install the required dependencies on a Debian system, run

    sudo apt-get install libgmp-dev libboost-filesystem-dev libboost-program-options-dev libboost-python-dev libcgal-dev

To install the required dependencies on OSX using Homebrew, run

    brew install boost --with-python
    brew install gmp
    brew install cgal

If you have previously installed `boost` without the `boost_python` component, the first command should be called with `reinstall` instead of `install`.

This project also depends on the [NLTemplate] string templating library.
It has been included here under the terms of the MIT license.


### Building the project code
Next, you will need to build our custom version of `GTFold` and the parametrizer.
To do so, simply run `make` from the `iB4e-GTfold-parametrizer` directory.
If you have multiple cores or processors, you can build in parallel by running

    nice make -j

instead.
The resulting binaries can be found in the `bin/` subdirectory.

## Updating

If you used Git to download your copy of this software, you can update it easily.
Just run `git pull` in a terminal from anywhere inside the repository to fetch the latest version.

## Usage

Given a FASTA file representing an RNA sequence, the program will produce a Sage file representing the regions of ℝ<sup>3</sup> in which given structures are optimal.

To run the calculation on the sequence in `test.fasta`, type

    bin/parametrizer test.fasta

The result will be a file `test.polytope.sage` containing the required Sage commands and a directory called `test` containing structure files representing the MFE structures for each set of parameters.

We also supply a program which can be used to find an MFE structure for a single set of parameters.
To use it on the sequence in `test.fasta` with parameters `A`, `B`, `C`, and `D`, type

    bin/pmfe-param test.fasta -a A -b B -c C -d D

The result will be printed to your terminal.

For more information about either program, run it with the `-h` option.

## Python interface

The project also includes a Python testing interface called `pyparam`.
To use it, just run the following from your favorite Python shell or Sage:

    import pyparam

You can now construct `pyparam.ParameterVector` and `pyparam.ScoreVector` objects (representing folding parameters and structure scores respectively), as well as call the `pyparam.get_mfe_score()` method to run `pmfe` (our modified `gtmfe`) on a sequence of your choice.
See the files `pyparam/example.*.py` for examples of how to use these methods.

## iB4e

This project includes "BBPolytope.h", a headers-only implementation of Huggins' `iB4e` algorithm.
It can be found in the `iB4e` subdirectory.
Note that it requires a replacement for one header file in the CGAL library; this is a small but necessary modification to support the algorithm.

## License

The source for iB4e-GTfold-parametrizer is released under the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

[macports]: //www.macports.org/
[openmp]: http://openmp.org/
[opemmp-dl]: http://openmp.org/wp/openmp-compilers/
[gmp]: //gmplib.org/
[gmp-dl]: //gmplib.org/#DOWNLOAD
[sage]: //sagemath.org
[cgal]: //www.cgal.org
[boost]: //www.boost.org
[boost-getstarted]: //www.boost.org/doc/libs/1_57_0/more/getting_started/unix-variants.html
[cmake]: //www.cmake.org/download/
[NLTemplate]: //github.com/catnapgames/NLTemplate
[gtfold]: //gtfold.sourceforge.net/
