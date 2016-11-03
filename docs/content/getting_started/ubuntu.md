# Ubuntu

{!content/getting_started/minimum_requirements.md!}

---
## Pre-Reqs
* Install the following using apt-get

    <pre>
    sudo apt-get install build-essential \
    gfortran \
    tcl \
    git \
    m4 \
    freeglut3 \
    doxygen \
    libblas-dev \
    liblapack-dev \
    libx11-dev \
    libnuma-dev \
    libcurl4-gnutls-dev \
    zlib1g-dev \
    libhwloc-dev
    </pre>

* Download one of our redistributables according to your version of Ubuntu

    * Ubuntu 16.04: !MOOSEPACKAGE arch=ubuntu16.04 return=link!
    * Ubuntu 14.04: !MOOSEPACKAGE arch=ubuntu14.04 return=link!
    * Ubuntu 12.04: !MOOSEPACKAGE arch=ubuntu12.04 return=link!

{!content/getting_started/install_redistributable_deb.md!}
{!content/getting_started/clone_moose.md!}
{!content/getting_started/build_libmesh.md!}
{!content/getting_started/conclusion.md!}