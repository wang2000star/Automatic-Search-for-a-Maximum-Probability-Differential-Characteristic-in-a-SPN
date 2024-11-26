cmake_minimum_required(VERSION 3.10)

# 项目信息
project(Automatic-Search-for-a-Maximum-Probability-Differential-Characteristic-in-a-SPN)

# 设置编译选项
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")

# 添加可执行文件
add_executable(main main.c fenzhidingjie.c)

# 链接 math 库
target_link_libraries(main m)