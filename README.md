# Zephyr
*Zephyr* is a secure and non-interactive two-party
inference for transformers. This code is built upon the code of NEXUS scheme,https://github.com/zju-abclab/NEXUS. The differences located 
in the src/matrix_mul.cpp,src/ckks_evaluator.h and src/main.cpp files.

<br/>

## Installing SEAL
HAWK uses a modified (bootstrapping-friendly) version of [Microsoft SEAL version 4.1](https://github.com/microsoft/SEAL/tree/4.1.2).

Please install the SEAL library `thirdparty/SEAL-4.1-bs` included as part of this repository (follow the instructions in the above link) before compiling HAWK.

<br/>

## Compiling HAWK
Once Microsoft SEAL 4.1 is installed, to build HAWK simply run:

```bash
mkdir build
cd build
cmake ..
make
```

This should produce two binary executables inside the `build` directory:
- `bin/main`
- `bin/bootstrapping`

<br/>




