
```
# See gHMM download instructions:

http://ghmm.sourceforge.net/download.html

# Check most recent version of code from repository:

[pcahan1@localhost ~]$ svn co https://ghmm.svn.sourceforge.net/svnroot/ghmm/trunk/ghmm
Checked out revision 2209.

[pcahan1@localhost ~]$ mkdir ghmm8
[pcahan1@localhost ~]$ mv ghmm ghmm8/

# Set shell environment variables manually as below or modify your .shell profile file

[pcahan1@localhost ghmm8]$ PYTHONPATH=/home/pcahan1/ghmm8/lib/python2.3/site-packages
[pcahan1@localhost ghmm8]$ export PYTHONPATH
[pcahan1@localhost ghmm8]$ LD_LIBRARY_PATH=/home/pcahan1/ghmm8/lib
[pcahan1@localhost ghmm8]$ export LD_LIBRARY_PATH
[pcahan1@localhost ghmm8]$ PATH=${PATH}/home/pcahan1/ghmm8/bin
[pcahan1@localhost ghmm8]$ export PATH

# See install instructions:
http://ghmm.sourceforge.net/installation.html

# Replace model.c with:

http://wuhmm.googlecode.com/files/model.c

# Compile and install ghmm
[pcahan1@localhost ghmm8]$ cd ghmm/
[pcahan1@localhost ghmm]$ ./autogen.sh
[pcahan1@localhost ghmm]$ ./configure --prefix=/home/pcahan1/ghmm8
[pcahan1@localhost ghmm]$ make
[pcahan1@localhost ghmm]$ make install
###lots of output including:
PATH="$PATH:/sbin" ldconfig -n /home/pcahan1/ghmm8/lib
----------------------------------------------------------------------
Libraries have been installed in:
   /home/pcahan1/ghmm8/lib

If you ever happen to want to link against installed libraries
in a given directory, LIBDIR, you must either use libtool, and
specify the full pathname of the library, or use the `-LLIBDIR'
flag during linking and do at least one of the following:
   - add LIBDIR to the `LD_LIBRARY_PATH' environment variable
     during execution
   - add LIBDIR to the `LD_RUN_PATH' environment variable
     during linking
   - use the `-Wl,--rpath -Wl,LIBDIR' linker flag
   - have your system administrator add LIBDIR to `/etc/ld.so.conf'

See any operating system documentation about shared libraries for
more information, such as the ld(1) and ld.so(8) manual pages.
----------------------------------------------------------------------

# Install python wrappers

[pcahan1@localhost ghmm]$ cd ghmmwrapper/
[pcahan1@localhost ghmmwrapper]$ swig -python -nodefault ghmmwrapper.i
[pcahan1@localhost ghmmwrapper]$ python setup.py build --debug
[pcahan1@localhost ghmmwrapper]$ python setup.py install --prefix=/home/pcahan1/ghmm8
```