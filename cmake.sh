#!/bin/bash

cmake -B build\
  -D CMAKE_BUILD_TYPE=Release\
  -D ENABLE_HIP=ON\
  -D HIP_DIR=/opt/rocm/lib/cmake/hip\
  -D CMAKE_HIP_ARCHITECTURES=gfx90a
