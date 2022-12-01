![WaveSimC](WaveSimCLogo.png)

# COMSW4995 Final Project: WaveSimC

This is the repository for our final project.

# Building

## Install the boost library

It is important to install the boost library before building the project. The boost library is used for data structures and algorithms. The boost library can be installed using the following command on ubuntu:

```bash
sudo apt-get install libboost-all-dev
```

For Mac:

```bash
brew install boost
```

## Build the project

```bash
mkdir build
cd build
cmake ..
make Main
```

## Running

```bash
./Main
```

## Building the documentation

Docs building script:

```bash
./compileDocs.sh
```

```bash
doxygen dconfig
cd documentation/latex
pdflatex refman.tex
cp refman.pdf ../WaveSimC-0.8-doc.pdf
```
