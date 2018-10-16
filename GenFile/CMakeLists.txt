project(GLPK CXX)

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
set(Header_Files
    ../GurungRay/parseBenchmark.h
)
source_group("Header Files" FILES ${Header_Files})

set(Source_Files
    ../GurungRay/glpk_lp_solver.cpp
    ../GurungRay/parseBenchmark.cpp
    GPU_Tester.cpp
)
source_group("Source Files" FILES ${Source_Files})

set(ALL_FILES  ${Source_Files} ${Header_Files})

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
          $<$<CONFIG:Debug>:${CMAKE_SOURCE_DIR}/include;${BOOST_INC};${GLPK_INC};${CMAKE_SOURCE_DIR}/GurungRay>
        $<$<CONFIG:Release>:${CMAKE_SOURCE_DIR}/include;${BOOST_INC};${GLPK_INC};${CMAKE_SOURCE_DIR}/GurungRay>
    )
endif()

################################################################################
# Compile definitions
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_compile_definitions(${PROJECT_NAME} PRIVATE
          $<$<CONFIG:Debug>:WIN64;_CRT_SECURE_NO_WARNINGS;_DEBUG;_CONSOLE;UNICODE;_UNICODE>
        $<$<CONFIG:Release>:WIN64;_CRT_SECURE_NO_WARNINGS;NDEBUG;_CONSOLE;UNICODE;_UNICODE>
    )
endif()

################################################################################
# Compile and link options
################################################################################
if(MSVC)
    if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
        target_compile_options(${PROJECT_NAME} PRIVATE
              $<$<CONFIG:Debug>:/Od;/MDd;/W3;/EHsc;/Y->
            $<$<CONFIG:Release>:/GL;/Od;/Oi;/MD;/Gy;/W3;/EHsc;/Y->
        )
    endif()
    if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
        set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS_DEBUG "/DEBUG /INCREMENTAL")
        set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS_RELEASE "/LTCG /DEBUG /INCREMENTAL:NO")
    endif()
endif()

################################################################################
# Dependencies
################################################################################
# Link with other additional libraries.
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_link_libraries(${PROJECT_NAME}
        $<$<CONFIG:Release>:glpk>
    )
endif()
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_link_directories(${PROJECT_NAME}
        CONDITION   $<CONFIG:Debug> "${GLPK_LIB}"
        CONDITION $<CONFIG:Release> "${GLPK_LIB}"
    )
endif()

################################################################################
# OpenMP
################################################################################
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

