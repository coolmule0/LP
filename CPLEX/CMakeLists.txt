project(CPLEX CXX)

################################################################################
# Source groups
################################################################################
set(Header_Files
    ../LP/FileIO.h
)
source_group("Header Files" FILES ${Header_Files})

set(Source_Files
    ../LP/FileIO.cpp
    cplex.cpp
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
          $<$<CONFIG:Debug>:${OUT_DIR}${TARGET_NAME}$ENV{TargetExt}>
        $<$<CONFIG:Release>:${OUT_DIR}${TARGET_NAME}$ENV{TargetExt}>
        $<$<NOT:$<OR:$<CONFIG:Debug>,$<CONFIG:Release>>>:${PROJECT_NAME}>
    )
endif()
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${TARGET_NAME})

################################################################################
# Output directory
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    string(CONCAT OUT_DIR
          $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/>
        $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/>
        $<$<NOT:$<OR:$<CONFIG:Debug>,$<CONFIG:Release>>>:${CMAKE_SOURCE_DIR}/${CMAKE_VS_PLATFORM_NAME}/$<CONFIG>>
    )
endif()
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${OUT_DIR})

################################################################################
# Include directories
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_include_directories(${PROJECT_NAME} PUBLIC
          $<$<CONFIG:Debug>:$ENV{CPLEX_STUDIO_DIR1271}/cplex/include;$ENV{CPLEX_STUDIO_DIR1271}/concert/include;${CMAKE_SOURCE_DIR}/include;$ENV{CUDA_PATH}/include>
        $<$<CONFIG:Release>:$ENV{CPLEX_STUDIO_DIR1271}/cplex/include;$ENV{CPLEX_STUDIO_DIR1271}/concert/include;${CMAKE_SOURCE_DIR}/include;$ENV{CUDA_PATH}/include>
    )
endif()

################################################################################
# Compile definitions
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_compile_definitions(${PROJECT_NAME} PRIVATE
          $<$<CONFIG:Debug>:IL_STD;_DEBUG;_CONSOLE;UNICODE;_UNICODE>
        $<$<CONFIG:Release>:IL_STD;NDEBUG;_CONSOLE;UNICODE;_UNICODE>
    )
endif()

################################################################################
# Compile and link options
################################################################################
if(MSVC)
    if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
        target_compile_options(${PROJECT_NAME} PRIVATE
              $<$<CONFIG:Debug>:/Od;/MDd;/W3;/EHsc;/Y->
            $<$<CONFIG:Release>:/GL;/O2;/Oi;/MD;/Gy;/W3;/EHsc;/Y->
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
          $<$<CONFIG:Debug>:cudart$<SEMICOLON>$(CPLEX_STUDIO_DIR1271)\cplex\lib\x64_windows_vs2015\stat_mdd\cplex1271$<SEMICOLON>$(CPLEX_STUDIO_DIR1271)\cplex\lib\x64_windows_vs2015\stat_mdd\ilocplex$<SEMICOLON>$(CPLEX_STUDIO_DIR1271)\concert\lib\x64_windows_vs2015\stat_mdd\concert>
        $<$<CONFIG:Release>:cudart$<SEMICOLON>$(CPLEX_STUDIO_DIR1271)\cplex\lib\x64_windows_vs2015\stat_mda\cplex1271$<SEMICOLON>$(CPLEX_STUDIO_DIR1271)\cplex\lib\x64_windows_vs2015\stat_mda\ilocplex$<SEMICOLON>$(CPLEX_STUDIO_DIR1271)\concert\lib\x64_windows_vs2015\stat_mda\concert>
    )
endif()
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_link_directories(${PROJECT_NAME}
        CONDITION   $<CONFIG:Debug> "$ENV{CUDA_PATH}/lib/x64"
        CONDITION $<CONFIG:Release> "$ENV{CUDA_PATH}/lib/x64"
    )
endif()

