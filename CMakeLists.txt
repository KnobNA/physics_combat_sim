cmake_minimum_required(VERSION 3.15)
project(RobotCombatSim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. Find ROS 2 (Must be installed on the host OS, e.g., Humble or Galactic)
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)

# 2. Add Project Chrono (Enable the ROS and FEA modules)
set(ENABLE_MODULE_ROS ON CACHE BOOL "Enable Chrono ROS" FORCE)
set(ENABLE_MODULE_FEA ON CACHE BOOL "Enable Chrono FEA" FORCE)
add_subdirectory(external/chrono)

# 3. Add Magnum Engine (Enable SDL2 application and OpenGL)
set(WITH_SDL2APPLICATION ON CACHE BOOL "" FORCE)
set(WITH_OPENGLTESTER ON CACHE BOOL "" FORCE)
add_subdirectory(external/corrade)
add_subdirectory(external/magnum)

# 4. Add Assimp for CAD importing
add_subdirectory(external/assimp)
add_subdirectory(external/magnum-plugins)

# 5. Define your simulation executable
add_executable(combat_sim 
    src/main.cpp 
    # Add other src files here later
)

# 6. Link everything together
target_link_libraries(combat_sim PRIVATE
    Chrono::Engine
    Chrono::ROS
    Magnum::Magnum
    Magnum::Sdl2Application
    Magnum::GL
    rclcpp::rclcpp
)