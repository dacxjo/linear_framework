cmake_minimum_required(VERSION 3.20)
project(BlandonTorrezDavid_prac1 C)

set(CMAKE_C_STANDARD 11)

add_executable(BlandonTorrezDavid_prac1 main.c funs_linalg.c funs_linalg.h)
add_executable(trisup main_trisup.c funs_linalg.c funs_linalg.h)
