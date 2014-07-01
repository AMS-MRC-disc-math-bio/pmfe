## Synopsis

iB4e-GTFold-parametrizer determines the sensitivity of RNA folding to multibranch-loop-related parameters.

## Prerequisites

You will need the `gtmfe` and `RNAScoring` programs from the [GTFold][gtfold] project.
At present, we recommend using our [fork][mrc-gtfold], especially if you are running a Mac.
To do this, run the following in your terminal:

```
git clone https://github.com/AMS-MRC-disc-math-bio/gtfold.git
cd gtfold/gtfold-mre
make
sudo make install
cd ../rna-scoring
make
sudo make install
```

If you do not have Git on your system, you can instead download a ZIP archive of the current state of the repository.
Extract it, then proceed from the second line of the above instructions.

You will also need a working Python environment.

## Installation

This project is under active development, so we recommend downloading it using Git.
To do this, run the following in your terminal:

```
git clone https://github.com/AMS-MRC-disc-math-bio/iB4e-GTfold-parametrizer.git
cd iB4e-GTfold-parametrizer/iB4e
make
```

This will download the code and build our custom version of iB4e.

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

[gtfold]: https://github.com/gtfold/gtfold
[macports]: https://www.macports.org/
[mrc-gtfold]: https://github.com/AMS-MRC-disc-math-bio/gtfold
