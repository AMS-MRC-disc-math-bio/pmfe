## Warning: old repo

**This version of pmfe is unmaintained. See https://github.com/gtDMMB/pmfe for active development.**

## Synopsis

`pmfe` determines the sensitivity of RNA folding to multibranch-loop-related parameters.
It relies on an implementation of Zucker's dynamic-programming algorithm for RNA secondary structure prediction, using adapted code from the `gtmfe` program from the [GTFold][gtfold] project.

[![DOI](https://zenodo.org/badge/12672/AMS-MRC-disc-math-bio/pmfe.svg)](https://zenodo.org/badge/latestdoi/12672/AMS-MRC-disc-math-bio/pmfe)

## Docker container

This project is available as a [Docker][docker] container.
If you wish to *run* `pmfe` but do not need to modify and recompile the code, we recommend this approach, as it makes it easy to ensure you have all the dependencies.

You can download the latest image by running the following in your shell:

    docker pull agdphd/pmfe

This command can be used in the future to update your image to the latest version.

You can then enter the image by running the following in your shell:

    docker run -it agdphd/pmfe

This will give you a shell in the `pmfe` directory with the binaries built and ready to go.

WARNING: Each time you call `docker run […]`, you get a clean new copy of the image!
Any data you produced in previous runs will *not* be present.
Thus, be sure to copy out anything important before exiting the image.

## Getting started with code

If you want to modify the source and compile `pmfe` yourself, the source is available here.
This project is under active development, so we recommend downloading it using Git.
To do this, run the following in your terminal:

```
git clone https://github.com/AMS-MRC-disc-math-bio/pmfe.git
```

This will download the code and extract it into a directory called `pmfe`.

### Dependencies
This project depends on the [CGAL][cgal] computational geometry library, several libraries from the [Boost][boost] project, and the [GMP][gmp] arbitrary-precision arithmetic library.
To install the required dependencies on a Debian system, run

    sudo apt-get install libgmp-dev libboost-filesystem-dev libboost-program-options-dev libboost-log-dev libcgal-dev

To install the required dependencies on OSX using Homebrew, run

    sudo brew install boost cgal gmp

This project also depends on the [Catch][catch] unit testing library.
It has been included here under the terms of the Boost Software License.


### Building the project code
Next, you will need to build our custom version of `GTFold` and the parametrizer.
To do so, simply run `make` from the `pmfe` directory.
If you have multiple cores or processors, you can build in parallel by running

    nice make -j

instead.
The resulting binaries can be found in the base directory with the `pmfe-` prefix.

## Updating

If you used Git to download your copy of this software, you can update it easily.
Just run `git pull` in a terminal from anywhere inside the repository to fetch the latest version.
Afterwards, run `make` again to build the new version.

## Usage
This project builds programs that perform a variety of tasks.

### `pmfe-findmfe`
Given a FASTA file representing an RNA sequence and (optionally) some modified values for the Turner99 multibranch loop parameters, the `pmfe-findmfe` program will generate a secondary structure which minimizes free energy.
For example, to use it on the sequence in `test_seq/tRNA/c.diphtheriae_tRNA.fasta` with parameters `A`, `B`, `C`, and `D`, type

    pmfe-findmfe test_seq/tRNA/c.diphtheriae_tRNA.fasta -a A -b B -c C -d D

The result will be printed to your terminal.

### `pmfe-subopt`
Given a FASTA file representing an RNA sequence, an energy gap δ, and (optionally) some modified values for the Turner99 multibranch loop parameters, the `pmfe-subopt` program will generate all secondary structures with energy within δ of the minimum.
To use it on the sequence in `test_seq/tRNA/c.diphtheriae_tRNA.fasta` with parameters `A`, `B`, `C`, and `D` and energy gap δ, type

    pmfe-subopt test_seq/tRNA/c.diphtheriae_tRNA.fasta -a A -b B -c C -d D --delta δ

The result will be saved in `test_seq/tRNA/c.diphtheriae_tRNA.rnasubopt`.

### `pmfe-parametrizer`
Given a FASTA file representing an RNA sequence, the `pmfe-parametrizer` program will generate the polytope which is the convex hull of all secondary structures on that sequence in ℚ⁴.
To use it on the sequence in `test_seq/tRNA/c.diphtheriae_tRNA.fasta`, type

    pmfe-parametrizer test_seq/tRNA/c.diphtheriae_tRNA.fasta

The result will be saved in `test_seq/tRNA/c.diphtheriae_tRNA.rnapoly`, which can be read directly or used with `rna_poly.py` to produce a Sage polytope for further investigation.

### `pmfe-tests`
The `pmfe-tests` program runs a suite of unit tests.

## iB4e

This project includes "BBPolytope.h", a headers-only implementation of Huggins' `iB4e` algorithm.
It can be found in the `iB4e` subdirectory.
Note that it requires a replacement for one header file in the CGAL library; this is a small but necessary modification to support the algorithm.
In the future, this will be updated to use the new Triangulations library in CGAL, and the modification will no longer be necessary.

## License

The source for pmfe is released under the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

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
[gtfold]: //gtfold.sourceforge.net/
[docker]: //docker.io/
[catch]: //github.com/philsquared/Catch
