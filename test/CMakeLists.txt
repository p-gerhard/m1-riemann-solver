cmake_minimum_required(VERSION 3.5)

# Build the example
add_executable(solve_example solve_example.c)
target_compile_options(solve_example PRIVATE $<$<C_COMPILER_ID:GNU>:-Wall -O3>)
target_link_libraries(solve_example PRIVATE m1_riemann m)