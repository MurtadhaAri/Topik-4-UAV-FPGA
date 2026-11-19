cmake_minimum_required(VERSION 3.10)
project(System4_FixedPoint C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
add_compile_options(-Wall -Wextra -Wpedantic)

include_directories(include include/core include/model include/analysis third_party/icd)

file(GLOB S4_SRC src/core/*.c src/model/*.c src/analysis/*.c)
add_library(s4core STATIC ${S4_SRC})
target_link_libraries(s4core m)

add_executable(run_analysis tools/run_analysis.c)
target_link_libraries(run_analysis s4core)

add_executable(system4_proof tests/proof/system4_proof.c)
target_link_libraries(system4_proof s4core)

enable_testing()
foreach(t test_fixed_math test_divider test_cordic test_adapter test_quant_stress test_models)
  add_executable(${t} tests/unit/${t}.c)
  target_include_directories(${t} PRIVATE tests/unit)
  target_link_libraries(${t} s4core)
  add_test(NAME ${t} COMMAND ${t})
endforeach()
add_test(NAME system4_proof COMMAND system4_proof ${CMAKE_SOURCE_DIR}/results)
