
set(MESHSRC ${PDA_SOURCE}/meshing_scripts)
set(DECOMPSRC ${CMAKE_CURRENT_SOURCE_DIR}/../meshing_scripts)

# optional flags for limiting/expanding build
set(TESTWENO3 TRUE)
option(PARTESTS "" OFF)
add_compile_definitions(SCHWARZ_SAVE_TEMPDIR)

# include demoapps headers and Schwarz routines
include_directories(
  ${PDA_SOURCE}/tpls/eigen3
  ${PDA_SOURCE}/include
  ${PRESSIO_SOURCE}/include
  ${PRESSIO_OPS_INCLUDE_DIR}
  ${PRESSIO_LOG_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_compile_definitions(PRESSIO_ENABLE_TPL_EIGEN PRESSIO_ENABLE_DEBUG_PRINT PRESSIO_ENABLE_LOGGING)

link_libraries(stdc++fs pthread)

# ---------------------------------------------------------
# 2d problems
# ---------------------------------------------------------

add_subdirectory(eigen_2d_burgers_outflow_implicit)
add_subdirectory(eigen_2d_burgers_outflow_implicit_schwarz)

add_subdirectory(eigen_2d_euler_riemann_implicit)
add_subdirectory(eigen_2d_euler_riemann_implicit_schwarz)

add_subdirectory(eigen_2d_swe_slip_wall_implicit)
add_subdirectory(eigen_2d_swe_slip_wall_implicit_schwarz)
add_subdirectory(eigen_2d_swe_slip_wall_implicit_roms)
add_subdirectory(eigen_2d_swe_slip_wall_implicit_roms_schwarz)
add_subdirectory(eigen_2d_swe_slip_wall_implicit_hproms)
add_subdirectory(eigen_2d_swe_slip_wall_implicit_hproms_schwarz)
add_subdirectory(eigen_2d_swe_slip_wall_implicit_mixed_schwarz)
add_subdirectory(eigen_2d_swe_slip_wall_implicit_nonoverlap_dd)
add_subdirectory(eigen_2d_swe_slip_wall_implicit_hproms_gpod)
add_subdirectory(eigen_2d_swe_slip_wall_implicit_hproms_gpod_schwarz)
add_subdirectory(eigen_2d_swe_slip_wall_implicit_roms_schwarz_icFile)

# misc subdirectories
if(PARTESTS)
  add_subdirectory(parallel)
endif()
