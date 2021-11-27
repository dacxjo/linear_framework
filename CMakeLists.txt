cmake_minimum_required(VERSION 3.20)
project(BlandonTorrezDavid_prac1 C)

set(CMAKE_C_STANDARD 99)

add_executable(BlandonTorrezDavid_prac1 main.c funs_linalg.c funs_linalg.h )
add_executable(trisup main_trisup.c funs_linalg.c funs_linalg.h)
add_executable(triinf main_triinf.c funs_linalg.c funs_linalg.h)
add_executable(gauss main_gauss.c funs_linalg.c funs_linalg.h)
add_executable(gaussLU main_gaussLU.c funs_linalg.c funs_linalg.h)
add_executable(factLU main_factLU.c funs_linalg.c funs_linalg.h)
add_executable(debug main_debug.c funs_linalg.c funs_linalg.h)
add_executable(gausspiv main_gausspiv.c funs_linalg.c funs_linalg.h)
add_executable(frank main_frank.c funs_linalg.c funs_linalg.h)
