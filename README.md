# SONO
An ad-hoc attempt of a modern implementation of the classical SONO algorithm. (Experimental)

#Installation

Please clone with the following command to make sure the submodules (dependencies) are pulled correctly:

  git clone --recurse-submodules https://github.com/HenrikSchumacher/SONO

Alternatively, you can clone as usual and then run the following to connect all submodules to their repos:

  git submodule update --init --recursive

The library is head-only. You can find a minimal working example in the subdirectory `SONO_MWE`. See the script `SONO_MWE/compile.sh` for compiling and linking instructions.
