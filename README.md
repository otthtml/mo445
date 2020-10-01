# MO445
This is a repo for MO445, a course given by professor Alexandre Falcão of Unicamp.


# How to run

## get all dependencies
sudo apt-get install libatlas-base-dev liblapack-dev libblas-dev \
sudo apt-get install make \
sudo apt-get install make-guile \
sudo apt-get install build-essential \
sudo apt-get install libz-dev \
sudo apt-get install libpng-dev \

## download and extract the examples (not required if cloning this repo)
curl https://www.ic.unicamp.br/~afalcao/mo445/libmo445.tar.bz2 --output libmo445.tar.bz2 \
tar -xf libmo445.tar.bz2

## navigate to the examples folder
cd libmo445/examples

## build adjacency.c (this will create a file in libmo445/bin)
make adjacency

## run the program (just a test, without any args)
../bin/adjacency 
