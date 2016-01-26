Project to show how to use Gaussian Processes when fitting data.  
Slide Deck [link](https://www.dropbox.com/s/fnw23qs5zdr3fnp/UdeMPizzaTalk.pdf?dl=0)

Using [homebrew](http://brew.sh) is probably the least painful way to get the FORTRAN software working on your Mac.

Install **gfortran**  
`brew install gfortran`  

Install **pgplot**  
`brew install pgplot`  
after installing pgplot add the following to your .bash_profile  
`export PGPLOT_DIR=/usr/local/Cellar/pgplot/5.2.2/share/`  
`export PGPLOT_DEV=/xserve`  

Install **cfitsio**  
`brew install cfitsio`


Install **XQuartz** to get X11 libraries 
Grab the install file from <http://www.xquartz.org>

At this point you should be able to run the configure script and make.  
`./configure`  
`make all`  
The configure script will bomb out if a package is missing.  If successful then the *bin* directory should contain all the binaries. 

To run the program use the following command in the bin directory.  

`./gptest datatest.dat`  

Press `enter` when asked to select graphics output device, the code should 
now run and present a series of plots similar to the presentation.