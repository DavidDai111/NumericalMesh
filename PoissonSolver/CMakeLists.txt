# 设置 CMake 的最低版本要求
cmake_minimum_required(VERSION 3.10)

# 设置项目名称和版本
project(PoissonSolver VERSION 1.0)

# 设置使用的 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 添加一个可执行文件
add_executable(PoissonSolver src/main.cpp src/Mesh.cpp)

# add_compile_options(-Wall -Wextra -pedantic -Werror)
