# LibGeodesicOPPGC
## Implementation of:
# Optimal Point-to-Point Geodesic Path Generation on Point Clouds
### Authors: Dr. Alexander Agathos(1) & Prof. Philip Azariadis(2)
### (1)University of the Aegean, (2)University of West Attica
<p align="center">
<img src="./images/bunny.jpg" alt="results">    
</p>
<p>The sofware requires the following packages:</p>
<ul>
    <li>CGAL</li>
    <li>Boost</li>
    <li>oneAPI MKL</li>
    <li>oneAPI Thread Building Blocks</li>
</ul>
<p>Installation on Windows with Microsoft Visual Studio (Community) 2022 is described in the following.</p>
<p>Please consider using the <a href="https://cmake.org/download/">cmake-gui tool</a>.</p>

<p><b><i>Boost installation</i></b></p>
<p>The most easy way to install boost is to visit <a href="https://sourceforge.net/projects/boost/">this</a> link. Please choose 
to download a file like boost_1_79_0-msvc-14.1-64.exe. This file installs the BOOST library in 
a directory that the user selects, let it be D:\Dev\boost_1_79_0.</p>
<p>The user then needs to set the environmental variables: <br />
BOOST_ROOT pointing to D:\Dev\boost_1_79_0<br />
BOOST_INCLUDEDIR pointing to D:\Dev\boost_1_79_0<br />
BOOST_LIBRARYDIR pointing to D:\Dev\boost_1_79_0\lib64-msvc-14.1<br/>
Now CMake will be able to find Boost in the system.
</p>
<p><b><i>CGAL installation</i></b></p>
<p>Suppose that <a href="https://github.com/CGAL/cgal/releases/download/v5.4/CGAL-5.4.zip">this</a> version
of CGAL is downloaded and unzipped to directory: D:\Dev\CGAL-5.4 then the following environmetal variable needs
to be set:<br />
CGAL_DIR pointing to D:\Dev\CGAL-5.4 <br/>
Now CMAKE can find CGAL.
</p>
<p><b><i>Intel oneAPI MKL & Thread Building Blocks</i></b></p>
<p>The most easy way to install both APIs is to install
<a href="https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html">Intel's oneAPI base toolkit software</a>. Once 
installed CMAKE will find it and build upon it.</p>
<p>The software has a test application testGeodesic which is called by</p>
<p>testGeodesic pointcloud.xyz</p> 
<p>The  argument is the point cloud to be smoothed which is a text file of points of the form "x y z nx ny nz" for each line (x, y, z the point coordinates, nx, ny, nz the normal coordinates).</p><br />

<p align="center">
<img src="./images/Logo_En.jpg" alt="logo" style="width: 85%;">    
</p>

