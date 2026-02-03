# 🌌 Laniakea
Laniakea is a Ultra-Light Dark Matter Dynamics Pseudo-Spectral solver written in C/C++.

![Soliton merging](Gallery/merge8.gif "Soliton merging") ![Gravitating solitons](Gallery/gravitating_solitons.gif "Gravitating solitons")

## 🚀 Motivation 
This project is based on the original work by Dr. Faber Edwards et al. ([arXiv:1807.04037](https://arxiv.org/abs/1807.04037)) and aims to develop a **faster and more scalable simulation framework** for Ultra-Light Dark Matter (ULDM) models.

Compared to the existing [PyUltraLight](https://github.com/auckland-cosmo/PyUltraLight) codebase, my implementation is designed for **high-performance computing**, enabling large-scale simulations on modern hardware.
At the current state, the algorithm is **faster** than the original implementation and provides **mostly identical** results. You can plot the results using the provided python scripts.

## 🛠️ Installation
To get started, clone the repository using the command and compile the project using make **AFTER** you install FFTW and FFmpeg: 
```
  git clone git@github.com:Ciriphys/Laniakea.git
  make
```

This project uses [FFmpeg](https://www.ffmpeg.org/) and [FFTW](https://www.fftw.org/), and currently supports macOS and Linux (tested on Ubuntu).
The codebase provides implementation for two different projects: SolitonSim and Laniakea. You'll need to run them both in this order, 
as SolitonSim provides the soliton data that Laniakea uses to compute the full simulation. All the simulations I've conducted are obtained
by running SolitonSim with the following arguments:
```
  ./run-sltn soliton.dat 900000 9.0
```

Laniakea uses a native configuration loader to avoid recompilation after each change.
You can find some examples in the `configs` folder. To run Laniakea use the following command:
```
  ./run-lnk <path-to-config-file>
```
at this point, a simulation should start and all information should appear on terminal.

### ⚙️ FFTW Installation Notes
On **both systems**, you can download and compile FFTW manually from the [official website](https://www.fftw.org/fftw3_doc/Installation-on-Unix.html).
The library should install automatically in `/usr/local/lib` unless specified otherwise. If you do **make sure to modify** the makefile accordingly.
Here is a little guide to install FFTW correctly:
```
./configure --enable-shared --enable-openmp --enable-thread
make && make install
```

### ⚙️ FFmpeg Installation Notes
FFmpeg is launched with a pipe on a separate thread on start by default. Here is a guide to install it on different operating systems:

#### macOS
On **macOS**, you can obtain FFmpeg via [Homebrew](https://brew.sh/):
```
brew install ffmpeg
```
Last tested version: macOS Tahoe 26.0.

#### Linux (Ubuntu)
On **Ubuntu**, you can obtain FFmpeg via apt:
```
sudo apt install ffmpeg
```
Last tested version: Ubuntu 25.10.

## 📋 Roadmap
Here is a planned roadmap of implemented and upcoming features:

- [x] Pseudo-spectral Schrödinger-Poisson and RKIV methods
- [x] Fully parallelized CPU computation using OpenMP.
- [x] Native video rendering using FFmpeg.
- [x] Native system loader using configuration files.
- [ ] Codebase refactoring.
- [ ] GPU computation (I do not own a modern GPU right now!).

## 📜 License
This project is licensed under the [GNU General Public License v3](LICENSE).

Portions of this project are based on code released under the BSD 3-Clause License. The original repository can be found [here](https://github.com/auckland-cosmo/PyUltraLight). 
See [`LICENSE.BSD`](LICENSE.BSD) for full terms.

## 🤝 Contributing
By contributing to this project (via pull request or other means), you agree to the terms of the [Contributor License Agreement](CLA.md).

Your contributions are welcome!

