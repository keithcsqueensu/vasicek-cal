// SPDX-License-Identifier: Apache-2.0
// Must NOT compile: a host-only function reached from device code. Proves the nvcc checks in
// tests/cuda/device_compile_check.cu would catch a device-incompatible callee (K-1, D-100).
#include "core/precision.hpp"

inline double host_only(double x) { return x + 1.0; }

VCAL_HD double calls_host(double x) { return host_only(x); }

__global__ void kernel(double* out) { out[0] = calls_host(out[0]); }
