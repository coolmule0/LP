project(LP CXX)

################################################################################
# Source groups
################################################################################
set(Headers
    Auxilary.h
    FileIO.h
)
source_group("Headers" FILES ${Headers})

set(Sources
    FileIO.cpp
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
          $<$<CONFIG:Debug>:${CMAKE_SOURCE_DIR}\/include>
        $<$<CONFIG:Release>:${CMAKE_SOURCE_DIR}\/include>
    )
endif()

################################################################################
# Compile definitions
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_compile_definitions(${PROJECT_NAME} PRIVATE
          $<$<CONFIG:Debug>:WIN64;_DEBUG;_CONSOLE;_MBCS>
        $<$<CONFIG:Release>:WIN64;NDEBUG;_CONSOLE;_MBCS>
    )
endif()

################################################################################
# Compile and link options
################################################################################
if(MSVC)
    if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
        target_compile_options(${PROJECT_NAME} PRIVATE
              $<$<CONFIG:Debug>:/Od;/MDd;/W3;/EHsc>
            $<$<CONFIG:Release>:/GL;/O2;/Oi;/MD;/Gy;/W3;/EHsc>
        )
    endif()
    if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
        set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS_DEBUG "/DEBUG /INCREMENTAL")
        set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS_RELEASE "/LTCG /DEBUG")
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
          $<$<CONFIG:Debug>:cudart$<SEMICOLON>kernel32$<SEMICOLON>user32$<SEMICOLON>gdi32$<SEMICOLON>winspool$<SEMICOLON>comdlg32$<SEMICOLON>advapi32$<SEMICOLON>shell32$<SEMICOLON>ole32$<SEMICOLON>oleaut32$<SEMICOLON>uuid$<SEMICOLON>odbc32$<SEMICOLON>odbccp32>
        $<$<CONFIG:Release>:cudart$<SEMICOLON>kernel32$<SEMICOLON>user32$<SEMICOLON>gdi32$<SEMICOLON>winspool$<SEMICOLON>comdlg32$<SEMICOLON>advapi32$<SEMICOLON>shell32$<SEMICOLON>ole32$<SEMICOLON>oleaut32$<SEMICOLON>uuid$<SEMICOLON>odbc32$<SEMICOLON>odbccp32>
    )
endif()
