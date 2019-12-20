# TEMPLATE

Submodules:
- https://github.com/Libavius/cmake-fortran-template.git
- https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git

## HOW TO USE

Clone this repository and initialize submodules:

```bash
git clone https://github.com/kazuki-iwanaga/my-template.git
# git submodule add https://github.com/Libavius/cmake-fortran-template.git
git submodule update --init --recursive
```

Run a docker container and enter in it:

```bash
docker-compose up -d myenv
docker-compose exec myenv /bin/sh
# myenv    : f90 and c++
# tex      : latex
# notebook : jupyter notebook
```

Make directory `build` and compile your sources in there:

```bash
mkdir build
cd build
cmake ..
make
```

See https://github.com/Libavius/cmake-fortran-template.git for more details.
