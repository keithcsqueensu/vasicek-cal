// SPDX-License-Identifier: Apache-2.0
// C++20-only construct. Must NOT compile under the C++17 header check (D-058).
#pragma once

template <class T>
concept Anything = true;
