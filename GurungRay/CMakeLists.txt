project(GurungRay LANGUAGES CXX CUDA)

################################################################################
# Ensure C++ is modern enough on CUDA compiler
################################################################################
set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS_RELEASE} -std=c++11")

################################################################################
# Library Location
################################################################################
set(GLPK_LIB "" CACHE PATH "Path to GLPK lib folder")
set(GLPK_INC "" CACHE PATH "Path to GLPK include/src directory")
set(BOOST_LIB "" CACHE PATH "Path to Boost lib folder")
set(BOOST_INC "" CACHE PATH "Path to Boost root directory")

################################################################################
# Source groups
################################################################################
set(Headers
    dataStructure.h
    parseBenchmark.h
    simplex.cuh
)
source_group("Headers" FILES ${Headers})

set(Sources
    glpk_lp_solver.cpp
    parseBenchmark.cpp
    GPU_Tester.cu
    simplex.cu
)
source_group("Sources" FILES ${Sources})

set(ALL_FILES  ${Headers} ${Sources})

################################################################################
# Target
################################################################################
add_executable(${PROJECT_NAME} ${ALL_FILES})

if(NOT ("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64"))
    message(WARNING "${CMAKE_VS_PLATFORM_NAME} arch is not supported!")
endif()

################################################################################
# Target name
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    string(CONCAT TARGET_NAME
          $<$<CONFIG:Debug>:${PROJECT_NAME}>
        $<$<CONFIG:Release>:${PROJECT_NAME}>
        $<$<NOT:$<OR:$<CONFIG:Debug>,$<CONFIG:Release>>>:${PROJECT_NAME}>
    )
endif()
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${TARGET_NAME})

################################################################################
# Output directory
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    string(CONCAT OUT_DIR
          $<$<CONFIG:Debug>:${CMAKE_SOURCE_DIR}/${CMAKE_VS_PLATFORM_NAME}/$<CONFIG>/>
        $<$<CONFIG:Release>:${CMAKE_SOURCE_DIR}/${CMAKE_VS_PLATFORM_NAME}/$<CONFIG>/>
        $<$<NOT:$<OR:$<CONFIG:Debug>,$<CONFIG:Release>>>:${CMAKE_SOURCE_DIR}/${CMAKE_VS_PLATFORM_NAME}/$<CONFIG>>
    )
endif()
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${OUT_DIR})

################################################################################
# Include directories
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_include_directories(${PROJECT_NAME} PUBLIC
          $<$<CONFIG:Debug>:${CMAKE_SOURCE_DIR}/include;${BOOST_INC};${GLPK_INC}>
        $<$<CONFIG:Release>:${CMAKE_SOURCE_DIR}/include;${BOOST_INC};${GLPK_INC}>
    )
endif()

################################################################################
# Compile definitions
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_compile_definitions(${PROJECT_NAME} PRIVATE
          $<$<CONFIG:Debug>:WIN64;_DEBUG;_CONSOLE;BOOST_UBLAS_USE_INDEXED_ITERATOR;_MBCS>
        $<$<CONFIG:Release>:WIN64;NDEBUG;_CONSOLE;BOOST_UBLAS_USE_INDEXED_ITERATOR;_MBCS>
    )
endif()

################################################################################
# Compile and link options
################################################################################
if(MSVC)
    if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
        target_compile_options(${PROJECT_NAME} PRIVATE
              $<$<CONFIG:Debug>:/Od;/MDd;/W3;/EHsc;/openmp>
            $<$<CONFIG:Release>:/GL;/O2;/Oi;/MD;/Gy;/W3;/openmp;/EHsc;/openmp>
        )
    endif()
    if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
        set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS_DEBUG "/DEBUG /INCREMENTAL")
        set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS_RELEASE "/LTCG /DEBUG /INCREMENTAL:NO")
    endif()
endif()

################################################################################
# Post build events
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    if(MSVC)
        add_custom_command_if(
            TARGET ${PROJECT_NAME}
            POST_BUILD
            COMMANDS
            COMMAND   $<CONFIG:Debug> echo copy "$ENV{CudaToolkitBinDir}\\cudart*.dll" "$<SHELL_PATH:${OUT_DIR}>"
            COMMAND   $<CONFIG:Debug> copy "$ENV{CudaToolkitBinDir}\\cudart*.dll" "$<SHELL_PATH:${OUT_DIR}>"
            COMMAND $<CONFIG:Release> echo copy "$ENV{CudaToolkitBinDir}\\cudart*.dll" "$<SHELL_PATH:${OUT_DIR}>"
            COMMAND $<CONFIG:Release> copy "$ENV{CudaToolkitBinDir}\\cudart*.dll" "$<SHELL_PATH:${OUT_DIR}>"
        )
    endif()
endif()

################################################################################
# Dependencies
################################################################################
# Link with other additional libraries.
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_link_libraries(${PROJECT_NAME}
          $<$<CONFIG:Debug>:boost_timer$<SEMICOLON>boost_system$<SEMICOLON>boost_chrono$<SEMICOLON>glpk>
        $<$<CONFIG:Release>:boost_timer$<SEMICOLON>boost_system$<SEMICOLON>boost_chrono$<SEMICOLON>glpk>
    )
endif()
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_link_directories(${PROJECT_NAME}
        CONDITION   $<CONFIG:Debug> "${BOOST_LIB}" "${GLPK_LIB}"
        CONDITION $<CONFIG:Release> "${BOOST_LIB}" "${GLPK_LIB}"
    )
endif()

