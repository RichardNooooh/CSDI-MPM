# Installing GMSH on TACC

Last checked on May 29, 2023 with Lonestar6, but it should work for most HPC environments. Not entirely sure if it's needed if you are using the Python package's `gmsh.jl`.

1. Download the latest source code tarball https://gmsh.info/src/gmsh-4.X.X-source.tgz
2. Transfer that file to TACC (`scp`) or download it on TACC directly with `wget`
3. Untar with -xvf
4. Inside that directory, run the following commands:
```
    mkdir build
    cd build
    cmake -DENABLE_BUILD_DYNAMIC=1 --DCMAKE_INSTALL_PREFIX=/your/install/path ..
    make
    make install
```

Note that the `--DCMAKE_INSTALL_PREFIX` is necessary since `/usr/local/bin` is not editable on TACC (or most HPC environments). Using either the `$HOME` or `$WORK` directories will work.

## Use in code

Simply modify `src/Util.jl` to include the `gmsh.jl` located in the installed GMSH directory. As an example, if it's in the `$HOME` directory, you can use `include("$HOME/gmsh/lib64/gmsh.jl")`.
