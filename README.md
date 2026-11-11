![Fortran](https://img.shields.io/badge/Fortran-%23734F96.svg?style=for-the-badge&logo=fortran&logoColor=white)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/1274358482.svg)](https://doi.org/10.5281/zenodo.24011229)

![EddyFlow Logo](img/app-logo.svg)

# Welcome to EddyFlow

EddyFlow is a powerful open source software application for processing eddy covariance data. It computes fluxes of water vapor (evapotranspiration), carbon dioxide, and other trace gases, and energy with the Eddy Covariance method.

EddyFlow is based on EddyPro, for which [LI-COR Biosciences](http://www.licor.com) has not published upstream updates for several years (No open source code since 2019 and no binaries since 2022). This repository is an actively maintained fork of EddyPro 6.2.2 that continues development while gratefully acknowledging LI-COR, ECO2S, and the original EddyPro authors for creating and releasing EddyPro. EddyPro originates from [ECO<sub>2</sub>S](http://gaia.agraria.unitus.it/eco2s), the Eddy COvariance COmmunity Software project, which was developed as part of the Infrastructure for Measurement of the European Carbon Cycle (IMECC-EU) research project. We gratefully acknowledge the [IMECC](http://imecc.ipsl.jussieu.fr/index.html) consortium, the ECO<sub>2</sub>S development team, the [University of Tuscia](https://www.unitus.it) (Italy) and scientists around the world who assisted with development and testing of the original version of this software.

## EddyFlow repositories

EddyFlow is a cross-platform application (Windows, Linux and macOS), which consists of a set of command line programs and a graphical user interface (GUI), kept in 4 repositories:

- [eddyflow-engine](https://github.com/kebasaa/eddyflow-engine) (current repository): The engine behind it all
- [eddyflow-gui](https://github.com/kebasaa/eddyflow-gui): The Qt6 GUI
- [eddyflow-documentation](https://github.com/kebasaa/eddyflow-documentation): The documentation of the software and its usage
  It provides the [documentation website](https://kebasaa.github.io/eddyflow-documentation/)
- [eddyflow-build-script](https://github.com/kebasaa/eddyflow-build-script): This builds the software binaries for Windows, Linux and macOS

## What the engine does

The engine is two command line programs. `eddyflow_rp` reads the raw data and writes the fluxes of every averaging period; `eddyflow_fcc` applies the flux and cospectral corrections that need the whole run. Beyond the processing options listed in the [GUI README](https://github.com/kebasaa/eddyflow-gui#data-processing-options-in-eddyflow), the engine:

- splits a run across all processor cores (`-j`): the planar-fit and time-lag pre-passes, the production pass and the flux computation of `eddyflow_fcc` run in worker processes, with output identical to a single-process run
- reads inputs from a Google Drive or Dropbox folder shared with "anyone with the link"; the output folder must be local
- runs EddyPro projects directly (see below)
- handles projects with any number of gases, each with its own unit, limits and moisture reference
- handles LI-COR GHG files whose acquisition frequency changes part way, assessing the spectra at each rate
- detects the time lag of slow gases at their own rate, and can let a gas that has no usable lag take that of another gas on the same analyser
- reports random uncertainty by several methods and, beside them, a limit of detection after Langford et al. (2015)
- corrects many errors that EddyPro 6.2.2 contained (see [Corrections to EddyPro 6.2.2](#corrections-to-eddypro-622))

## Installing EddyFlow

You can download EddyFlow from the GUI [releases](https://github.com/kebasaa/eddyflow-gui/releases) page. On Windows it is a portable `.7z` archive that does not require installation; Linux packages (`.deb`, `.rpm`, `.snap`) and a macOS `.dmg` are provided as well.

To build the EddyFlow engine separately, follow these instructions:

### Engine

The engine needs [gfortran](https://gcc.gnu.org/wiki/GFortran) (the GNU Fortran compiler), CMake 3.21 or later and Ninja. From the repository root, run:

    $ cmake --preset release
    $ cmake --build --preset release

This builds `eddyflow_rp` and `eddyflow_fcc` into `bin/<os>/`. The `checks` preset builds the same two with gfortran's subscript and runtime checks compiled in, into the same folder. The test drivers that some static checks run are separate targets, built into `obj/<os>/`, for example `cmake --build --preset release --target pwbref`.

On Windows, put the MinGW folder with `gfortran.exe` on `PATH` first. WinLibs' GCC builds include CMake and Ninja too.

## Running the engines from the command line

Use `eddyflow_rp` for Express mode, Advanced mode step 1, and Metadata Retriever runs. Use `eddyflow_fcc` for Advanced mode step 2 when flux/cospectral correction outputs are required.

From a shell, the commands are:

```bash
eddyflow_rp -c gui -s [linux,win] -e /path/to/eddyflow-env /path/to/project.eddyflow
eddyflow_fcc -c gui -s [linux,win] -e /path/to/eddyflow-env /path/to/project.eddyflow
```

Command line arguments:

- `-c`, `--caller`: use `gui` to mirror GUI behavior; omit it to use the console default.
- `-s`, `--system`: host operating system, one of `win`, `linux`, or `mac`.
- `-e`, `--environment`: EddyFlow environment directory containing folders such as `ini` and `tmp`; if omitted, the engine defaults to `..`.
- `-m`, `--mode`: `embedded` or `desktop` mode.
- `-j`, `--jobs`: how many worker processes a run may be split across. `0`, or leaving it out, uses as many as there are processor cores; `1` runs serially.
- `-v`, `--version`: print the version. `-h`, `--help`: print a short help.
- Final positional argument: path to the EddyFlow project file (`*.eddyflow`).

### Running an EddyPro project

A project written by EddyPro can be given to either engine as it stands:

```bash
eddyflow_rp -c gui -s win /path/to/project.eddypro
```

It is converted into an ordinary EddyFlow project and metadata pair beside the
file it came from - `project_ep_imported.eddyflow` and
`project_ep_imported.metadata` - and that pair is what runs. The EddyPro files
are only read; nothing is written back to them. A project whose first line is
`;EDDYPRO_PROCESSING` is recognised whatever its extension.

The conversion happens once. Run `eddyflow_fcc` afterwards against the same
`.eddypro` path and it reuses the imported pair rather than rebuilding it,
which is what lets it pick up the FLUXNET file `eddyflow_rp` recorded there.
Delete the imported pair to convert again after editing the EddyPro project.

The metadata is taken from the path the project's `proj_file` names, or - when
that path is not on this machine, which is usual for a project copied from
elsewhere - from the `.metadata` file of the project's own name beside it.

Both imported files are complete EddyFlow files, not merely readable ones: they
carry the settings EddyFlow added and EddyPro has no way to state, each at the
value the engine already applies when the key is absent. So the conversion
changes no result, and the GUI has nothing left to migrate when it opens them.

For a full Advanced mode run, run `eddyflow_rp` first. Run `eddyflow_fcc` afterwards only when the project requires flux/cospectral correction outputs.


## Utilities

To successfully run EddyFlow, the program installation folder must contain the following command line utilities under the 'bin' sub-directory:

- 7-zip (Windows only)
- Quazip
- Zlib

### 7-zip

[7-Zip](http://www.7-zip.org/) is a file archiver.

The console application consists of two files:
- 7z.dll
- 7z.exe

License: [LGPL](https://www.7-zip.org/license.txt).

## Using EddyFlow sample data

You can run EddyFlow using sample data files available at the [eddyflow-build-script](https://github.com/kebasaa/eddyflow-build-script) repository.

## Corrections to EddyPro 6.2.2

EddyFlow corrects errors that EddyPro 6.2.2 contained. These are the ones whose effect on results is the largest, as recorded in the [CHANGELOG](CHANGELOG), which has the details of each:

+ Time lag
  - The optimiser recognised a lag pinned at the edge of its search window only by exact equality, which almost never held, so such lags stayed in the data the nominal lag and window are fitted to and pulled the result towards the edge (8.1.8).
  - Sparse relative-humidity classes for water vapour were filled from neighbours still holding the error code, giving lags of about -5000 s and +10000 s that were accepted and then written outside their arrays (8.1.8).

+ Rotation
  - A planar-fit file read back by a later run shifted every sector boundary by the north offset, so wind near a boundary was rotated with the neighbouring sector's matrix (8.1.8).
  - The counterclockwise fill of sectors with too few periods (`pf_fix`) never ran, so such a sector was double-rotated without a message (8.1.8).

+ Quality control
  - The dropout test never raised its hard flag, and its percentages in the output read 0 (8.1.8).
  - The stationarity test added a sub-period's error code into the mean covariance like a value, giving a deviation of millions of percent and quality flag 2 for a period with a missing sub-period (8.1.8).

+ Spectral correction
  - A very short or mistyped analyser path length made the correction factor negative and tiny, so the flux flipped sign and shrank to nothing without a flag (8.1.8).

+ Day and night
  - The potential radiation took clock noon for solar noon. Where biomet has neither global radiation nor PPFD it decides day and night, so the NIGHT flag, and Burba's day or night regression for an LI-7500, were wrong at sites far from their time zone's meridian (8.1.8).

+ Footprint
  - Footprints were written for periods without u\*, z/L or L, from distances computed on the error code (8.1.8).

+ Density corrections
  - With the WPL correction switched off, only the designated hygrometer obeyed the switch; every other gas was still density-corrected (8.1.3).

+ Detrending
  - Running-mean and exponential detrending crashed when the time constant was longer than the averaging period (8.1.9).

+ Output files
  - The metadata file had two layouts depending on whether `eddyflow_fcc` follows `eddyflow_rp`, with columns under the wrong names and different units; EddyPro 7 has the same split (8.1.3).

How EddyFlow compares with EddyUH is in [docs/eddyuh-comparison.md](docs/eddyuh-comparison.md); that comparison was written for version 8.0.0.

## EddyFlow Trademark and Logo Policy

In order to help users who want to cite EddyFlow on posters or publications, we provide the [EddyFlow logo as vector graphic](img/app-logo.svg) for the proper use of the EddyFlow wordmark and logo.

## Want to Know More?

More information is available at the help website [EddyFlow help](https://kebasaa.github.io/eddyflow-documentation/).

See also the [CHANGELOG](CHANGELOG).

## How to Cite

Muller, J.D. (2026). *eddyflow-engine: EddyFlow flux processing engine*. DOI: https://doi.org/10.5281/zenodo.24011229

## License

This software is distributed under the [GNU General Public License (GPL) v3.0](LICENSE). Any modification of the code in this repository may only be released under the same license, and with attribution of authorship of the original code (i.e., citation above).
