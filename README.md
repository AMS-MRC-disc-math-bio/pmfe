## Synopsis

iB4e-GTFold-parametrizer determines the sensitivity of RNA folding to multibranch-loop-related parameters.

## Prerequisites

You will need the `gtmfe` and `RNAScoring` programs from the [GTFold][gtfold] project.
Download and build them according to the instructions provided at the project site.
Be sure to run `make install` so that these are available in your path.

You will also need a working Python environment.

## Installation

You can obtain the code in two ways:

1. By cloning this repository to your computer, or
1. By downloading a ZIP archive of the current state.

Once you have obtained the code, you will need to build our custom version of iB4e.
Run the `compile_iB4e` script in the `iB4e` directory to set it up.

## Usage

Given a FASTA file representing an RNA sequence, the program will produce a Sage file representing the regions of ℝ<sup>3</sup> in which given structures are optimal.

To run the calculation on the sequence in `testseq.fasta`, type

    ./iB4erunner.py -s testseq.fasta

The result will be a file `testseq.sage` containing the required Sage commands.

## Contributors

Let people know how they can dive into the project, include important links to things like issue trackers, irc, twitter accounts if applicable.

## License

The source for iB4e-GTfold-parametrizer is released under the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
later version.

The code in the `iB4e` folder is from the iB4e project by Peter Huggins, used and distributed here under the terms of the GNU General Public License.

[gtfold]: https://github.com/gtfold/gtfold
[macports]: https://www.macports.org/
