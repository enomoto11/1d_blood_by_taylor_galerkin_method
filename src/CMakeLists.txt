###################################################################################
#
# 1D flow solver
#
# Copyright (c) 2023 Biomechanics Lab.,
#                    Department of Mechanical Science and Bioengineering,
#                    Graduate School of Engineering Science,
#                    Osaka University.
# All rights reserved.
#
###################################################################################

if(enable_GLOG)
  set(compile_options_tmp ${base_compile_options} " -D_GLOG")
elseif()
  set(compile_options_tmp ${base_compile_options})
endif()

string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

set(main_target "TaylorGalerkin")

add_executable(${main_target} ${main_src})
target_sources(${main_target}
  PRIVATE
    main.cpp
    flow1D.cpp
    preprocess.cpp
)
target_compile_options(${main_target} PUBLIC ${compile_options_tmp})

target_link_libraries(${main_target} PRIVATE
  ${GLOG_LIBRARY}
  iomp5
)

target_include_directories(${main_target} PRIVATE
)

install(TARGETS ${main_target} DESTINATION ${CMAKE_INSTALL_PREFIX})