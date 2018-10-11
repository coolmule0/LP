project(GurungRay LANGUAGES CXX CUDA)

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
          $<$<CONFIG:Debug>:${CMAKE_SOURCE_DIR}/include;$ENV{BOOST_DIR};$ENV{GLPK_DIR}/src>
        $<$<CONFIG:Release>:${CMAKE_SOURCE_DIR}/include;$ENV{BOOST_DIR};$ENV{GLPK_DIR}/src>
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

################################################################################
# Dependencies
################################################################################
# Link with other additional libraries.
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_link_libraries(${PROJECT_NAME}
          $<$<CONFIG:Debug>:cudart$<SEMICOLON>kernel32$<SEMICOLON>user32$<SEMICOLON>gdi32$<SEMICOLON>winspool$<SEMICOLON>comdlg32$<SEMICOLON>advapi32$<SEMICOLON>shell32$<SEMICOLON>ole32$<SEMICOLON>oleaut32$<SEMICOLON>uuid$<SEMICOLON>odbc32$<SEMICOLON>odbccp32$<SEMICOLON>libboost_timer-vc140-mt-gd-x64-1_66$<SEMICOLON>libboost_system-vc140-mt-gd-x64-1_66$<SEMICOLON>libboost_chrono-vc140-mt-gd-x64-1_66$<SEMICOLON>glpk_4_64>
        $<$<CONFIG:Release>:cudart$<SEMICOLON>kernel32$<SEMICOLON>user32$<SEMICOLON>gdi32$<SEMICOLON>winspool$<SEMICOLON>comdlg32$<SEMICOLON>advapi32$<SEMICOLON>shell32$<SEMICOLON>ole32$<SEMICOLON>oleaut32$<SEMICOLON>uuid$<SEMICOLON>odbc32$<SEMICOLON>odbccp32$<SEMICOLON>libboost_timer-vc140-mt-x64-1_66$<SEMICOLON>libboost_system-vc140-mt-x64-1_66$<SEMICOLON>libboost_chrono-vc140-mt-x64-1_66$<SEMICOLON>glpk_4_64>
    )
endif()
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_link_directories(${PROJECT_NAME}
        CONDITION   $<CONFIG:Debug> "$ENV{BOOST_DIR}/stage/lib" "$ENV{GLPK_DIR}/w64"
        CONDITION $<CONFIG:Release> "$ENV{BOOST_DIR}/stage/lib" "$ENV{GLPK_DIR}/w64"
    )
endif()

