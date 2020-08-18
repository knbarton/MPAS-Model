# Python setup

from distutils.core import setup, Extension
from numpy import get_include

numpy_inc = get_include()		#  NumPy include path.
objs = "sht_init.o sht_kernels_a.o sht_kernels_s.o sht_fly.o sht_omp.o"
shtns_o = objs.split()			# transform to list of objects
libdir = "/usr/local"
if len(libdir) == 0:
	libdir = []
else:
	libdir = [libdir+"/lib"]
cargs = "-fopenmp"
largs = "-O3 -m64"
libs = "-lfftw3_omp -lfftw3 -lm -L/usr/projects/climate/SHARED_CLIMATE/software/grizzly/pio/1.7.2/gcc-5.3.0/openmpi-1.10.5/netcdf-4.4.1-parallel-netcdf-1.5.0/lib -lpio -L/usr/projects/climate/SHARED_CLIMATE/software/grizzly/netcdf/4.4.1/gcc-5.3.0/lib -lnetcdff -lnetcdf -L/usr/projects/climate/SHARED_CLIMATE/software/grizzly/parallel-netcdf/1.5.0/gcc-5.3.0/openmpi-1.10.5/lib -lpnetcdf"
libslist = libs.replace('-l','').split()	# transform to list of libraries

shtns_module = Extension('_shtns', sources=['shtns_numpy_wrap.c'],
	extra_objects=shtns_o, depends=shtns_o,
	extra_compile_args=cargs.split(),
	extra_link_args=largs.split(),
	library_dirs=libdir,
	libraries=libslist,
	include_dirs=[numpy_inc])

setup(name='SHTns',
	version='3.4.1',
	description='High performance Spherical Harmonic Transform',
	license='CeCILL',
	author='Nathanael Schaeffer',
	author_email='nschaeff@ujf-grenoble.fr',
	url='https://bitbucket.org/nschaeff/shtns',
	ext_modules=[shtns_module],
	py_modules=["shtns"],
	requires=["numpy"],
	)
