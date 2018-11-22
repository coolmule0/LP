project(CPLEX CXX)

################################################################################
# Library Location
################################################################################
set(CPLEX_DIR "" CACHE PATH "Path to CPLEX root directory, contains cplex and concert as sub directories")
set(CPLEX_CPLEX_LIB "" CACHE PATH "Path to CPLEX lib directory, contains library files within")
set(CPLEX_CONCERT_LIB "" CACHE PATH "Path to CPLEX lib directory, contains library files within")
set(CPLEX_VERSION "" CACHE bool "Version number of CPLEX, e.g. 1271")

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
          $<$<CONFIG:Debug>:${CPLEX_DIR}/cplex/include;${CPLEX_DIR}/concert/include;${CMAKE_SOURCE_DIR}/include;$ENV{CUDA_PATH}/include>
        $<$<CONFIG:Release>:${CPLEX_DIR}/cplex/include;${CPLEX_DIR}/concert/include;${CMAKE_SOURCE_DIR}/include;$ENV{CUDA_PATH}/include>
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
  if(MSVC)
    target_link_libraries(${PROJECT_NAME}
          $<$<CONFIG:Debug>:${CPLEX_CPLEX_LIB}/cplex.lib$<SEMICOLON>${CPLEX_CPLEX_LIB}/ilocplex.lib$<SEMICOLON>${CPLEX_CONCERT_LIB}/concert.lib>
        $<$<CONFIG:Release>:${CPLEX_CPLEX_LIB}/cplex.lib$<SEMICOLON>${CPLEX_CPLEX_LIB}/ilocplex.lib$<SEMICOLON>${CPLEX_CONCERT_LIB}/concert.lib>
    )
    else()
      target_link_libraries(${PROJECT_NAME}
          $<$<CONFIG:Debug>:${CPLEX_CPLEX_LIB}/cplex${CPLEX_VERSION}$<SEMICOLON>${CPLEX_CPLEX_LIB}/ilocplex$<SEMICOLON>${CPLEX_CONCERT_LIB}/concert>
        $<$<CONFIG:Release>:${CPLEX_CPLEX_LIB}/cplex${CPLEX_VERSION}$<SEMICOLON>${CPLEX_CPLEX_LIB}/ilocplex$<SEMICOLON>${CPLEX_CONCERT_LIB}/concert>
    )
  endif()
endif()
#if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
#    target_link_directories(${PROJECT_NAME}
#        CONDITION   $<CONFIG:Debug> "${CPLEX_CPLEX_LIB}" "${CPLEX_CONCERT_LIB}"
#        CONDITION $<CONFIG:Release> "${CPLEX_CPLEX_LIB}" "${CPLEX_CONCERT_LIB}"
#    )
#endif()

