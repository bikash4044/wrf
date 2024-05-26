export angel=$(pwd)
mkdir wrf
sudo apt-get update & sudo apt-get upgrade -y
sudo apt-get install gcc cpp csh perl git gfortran
cd wrf
mkdir Build_WRF
cd Build_WRF
mkdir test
mkdir LIBRARIES
cd LIBRARIES


pushd $(dirname $(which ncl))
cd ..
echo "export NCARG_ROOT=$(pwd)" >> ~/.bashrc
echo "export PATH=$NCARG_ROOT/bin:$PATH" >> ~/.bashrc

echo "export DIR=$(pwd)" >> ~/.bashrc
echo "export CC=gcc" >> ~/.bashrc
echo "export CXX=g++" >> ~/.bashrc
echo "export FC=gfortran" >> ~/.bashrc
echo "export FCFLAGS=-m64" >> ~/.bashrc
echo "export F77=gfortran" >> ~/.bashrc
echo "export FFLAGS=-m64" >> ~/.bashrc
echo "export JASPERLIB=$DIR/grib2/lib" >> ~/.bashrc
echo "export JASPERINC=$DIR/grib2/include" >> ~/.bashrc
echo "export LDFLAGS=-L$DIR/grib2/lib" >> ~/.bashrc
echo "export CPPFLAGS=-I$DIR/grib2/include" >> ~/.bashrc
source ~/.bashrc

popd
wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/netcdf-c-4.7.2.tar.gz
tar -xvf netcdf-c-4.7.2.tar.gz
cd netcdf-c-4.7.2
./configure --prefix=$DIR/netcdf --disable-dap --disable-netcdf-4 --disable-shared
make
make install
cd ..

echo "export PATH=$DIR/netcdf/bin:$PATH" >> ~/.bashrc
echo "export NETCDF=$DIR/netcdf" >> ~/.bashrc
echo "export CPPFLAGS=-I$DIR/netcdf/include" >> ~/.bashrc
echo "export LDFLAGS=-$DIR/netcdf/lib" >> ~/.bashrc

source ~/.bashrc
wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/netcdf-fortran-4.5.2.tar.gz
tar -xvf netcdf-fortran-4.5.2.tar.gz
cd netcdf-fortran-4.5.2
./configure --prefix=$DIR/netcdf --disable-dap --disable-netcdf-4 --disable-shared
make
make install
cd ..

# wget https://www.mpich.org/static/downloads/4.1.2/mpich-4.1.2.tar.gz
# tar -xvf mpich-4.1.2.tar.gz
# cd mpich-4.1.2
# ./configure --prefix=$DIR/mpich
# make
# make install
# echo "export PATH=$DIR/mpich/bin:$PATH" >> ~/.bashrc
# source ~/.bashrc
# cd ..


#zlib
wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/zlib-1.2.11.tar.gz
tar -xvf zlib-1.2.11.tar.gz
cd zlib-1.2.11
./configure --prefix=$DIR/grib2
make
make install
cd ..

wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/libpng-1.2.50.tar.gz
tar -xvf libpng-1.2.50.tar.gz
cd libpng-1.2.50
./configure --prefix=$DIR/grib2
make
make install
cd ..

wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/jasper-1.900.1.tar.gz
tar -xvf jasper-1.900.1.tar.gz
cd jasper-1.900.1
./configure --prefix=$DIR/grib2
make 
make install
cd ..

cd ../test
wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/Fortran_C_NETCDF_MPI_tests.tar
tar -xvf Fortran_C_NETCDF_MPI_tests.tar




cd $angel
unset angel
rm $(basename "$0")
