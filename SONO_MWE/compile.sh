 clang++                                        \
    -Wall                                       \
    -Wextra                                     \
    -std=c++20                                  \
    -O3                                         \
    -ffast-math                                 \
    -march=native                               \
    -mtune=native                               \
    -flto                                       \
    -o SONO_MWE                                 \
    main.cpp                                    \
