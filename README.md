## Synopsis

iB4e-GTFold-parametrizer determines the sensitivity of RNA folding to multibranch-loop-related parameters.

## Installation

This project is under active development, so we recommend downloading it using Git.
To do this, run the following in your terminal:

```
git clone https://github.com/AMS-MRC-disc-math-bio/iB4e-GTfold-parametrizer.git
```

This will download the code and extract it into a directory called `iB4e-GTfold-parameterizer`.

Next, you will need to build our custom versions of iB4e and GTFold.
For iB4e, run the following from the `iB4e-GTfold-parameterizer` directory:

```
cd iB4e-GTfold-parametrizer/iB4e
make
```

For GTFold, you will need to decide whether to link against the [OpenMP][openmp] library.
Using OpenMP will allow GTFold to take advantage of multiple cores or processors, speeding it up significantly.
However, it is not supported natively by the cgit ompiler included in OSX.

To link against OpenMP, run the following from the `iB4e-GTfold-parameterizer` directory:
```
cd ../gtfold/gtfold-mfe
python setup-with-openmp.py build_ext --inplace
```

To build without linking against OpenMP, run the following from the `iB4e-GTfold-parameterizer` directory:
```
cd ../gtfold/gtfold-mfe
python setup-without-openmp.py build_ext --inplace
```

## Updating

If you used Git to download your copy of this software, you can update it easily.
Just run `git pull` in a terminal from anywhere inside the repository to fetch the latest version.
Be sure to run `make` in the `iB4e` directory if the iB4e code has been updated!

## Usage

Given a FASTA file representing an RNA sequence, the program will produce a Sage file representing the regions of ℝ<sup>3</sup> in which given structures are optimal.

To run the calculation on the sequence in `testseq.fasta`, type

    ./iB4erunner.py -s testseq.fasta

The result will be a file `testseq.sage` containing the required Sage commands.

## License

The source for iB4e-GTfold-parametrizer is released under the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
later version.

The code in the `iB4e` folder is from the iB4e project by Peter Huggins, used and distributed here under the terms of the GNU General Public License.

[macports]: https://www.macports.org/
