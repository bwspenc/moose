# Fedora

{!content/getting_started/minimum_requirements.md!}

---
## Pre-Reqs
* Install the following using dnf

    <pre>
    sudo -E dnf install gcc \
    gcc-c++ \
    gcc-gfortran \
    tcl \
    tk \
    findutils \
    make \
    freeglut-devel \
    libXt-devel \
    libX11-devel \
    m4 \
    blas-devel \
    lapack-devel \
    git
    </pre>

* Download one our redistributables according to your version of Fedora

    * Fedora Core 23: !MOOSEPACKAGE arch=fedora23 return=link!

{!content/getting_started/install_redistributable_rpm.md!}
{!content/getting_started/clone_moose.md!}
{!content/getting_started/build_libmesh.md!}
{!content/getting_started/conclusion.md!}