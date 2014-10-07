## Synopsis

iB4e-GTFold-parametrizer determines the sensitivity of RNA folding to multibranch-loop-related parameters.

## Installation

This project is under active development, so we recommend downloading it using Git.
To do this, run the following in your terminal:

```
git clone https://github.com/AMS-MRC-disc-math-bio/iB4e-GTfold-parametrizer.git
```

This will download the code and extract it into a directory called `iB4e-GTfold-parameterizer`.

### Dependencies
The project depends on various libraries and functionality that are implemented in [Sage][sage].
If you do not already have Sage installed on your computer, you must download and install it before you can proceed.

Once you have installed Sage, you must make it available.
To see if this is the case, run `sage` from a terminal.
If you see a version string like `Sage Version 6.4.beta3, Release Date: 2014-09-10`, you're good to go!
If you see a `command not found` error, you need to configure your path.
The easiest way to do this is to link it into the `/usr/local/bin` directory, using the following command:
```
ln -s /path/to/sage/directory/sage /usr/local/bin
```

If you install Sage using the binary distribution (typical), you must also install the Sage version of the compiler suite GCC.
To do this, run the following from a terminal.
```
sage -i gcc
```
If you installed Sage from source or Git, you already have the Sage GCC and no action is required.


### Building the project code
Next, you will need to build our custom version of `GTFold`.
To do so, simply run `make` from the `iB4e-GTfold-parametrizer` directory.

## Updating

If you used Git to download your copy of this software, you can update it easily.
Just run `git pull` in a terminal from anywhere inside the repository to fetch the latest version.

## Usage

Given a FASTA file representing an RNA sequence, the program will produce a Sage file representing the regions of ℝ<sup>3</sup> in which given structures are optimal.

To run the calculation on the sequence in `testseq.fasta`, type

    ./parametrizer.py -s testseq.fasta

The result will be a file `testseq.polytope.sage` containing the required Sage commands and a directory called `testseq` containing structure files representing the MFE structures for each set of parameters.

## Testing

This project includes unit tests in the files `test_GTrunner.py` and `test_RNAscorer.py`.
Simply run each of these scripts from your shell to test the functionality of the associated program.
If any of the tests fail, please be sure your build of the C code in `gtmfe/` is up to date; if the tests continue to fail, please contact the author with information about your computer architecture, operating system, and compiler.

## License

The source for iB4e-GTfold-parametrizer is released under the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

[macports]: //www.macports.org/
[openmp]: http://openmp.org/
[opemmp-dl]: http://openmp.org/wp/openmp-compilers/
[gmp]: //gmplib.org/
[gmp-dl]: //gmplib.org/#DOWNLOAD
[sage]: //sagemath.org
