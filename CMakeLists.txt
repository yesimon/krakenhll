cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(krakenuniq CXX C)
find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
find_package(ZLIB REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wfatal-errors -pipe -std=c++11")

add_library(kraken
  src/krakendb.cpp
  src/quickfile.cpp
  src/krakenutil.cpp
  src/seqreader.cpp
  src/uid_mapping.cpp
  src/gzstream/gzstream.C
  src/hyperloglogplus.cpp)

target_link_libraries(kraken
  ${ZLIB_LIBRARIES})

set(exes
  classify
  db_shrink
  db_sort
  set_lcas
  count_unique
  query_taxdb
  build_taxdb
  dump_taxdb
  read_uid_mapping)
foreach(exe ${exes})
  add_executable(${exe} src/${exe}.cpp)
  target_link_libraries(${exe} kraken)
endforeach(exe)

add_executable(classifyExact src/classify.cpp)
target_link_libraries(classifyExact kraken)
target_compile_options(classifyExact PUBLIC -DEXACT_COUNTING)

install(TARGETS ${exes} DESTINATION libexec)

file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/bin ${CMAKE_INSTALL_PREFIX}/libexec)
file(GLOB scripts scripts/*)
foreach(s ${scripts})
  get_filename_component(basename ${s} NAME)
  file(COPY ${s} DESTINATION scripts)
  if(IS_DIRECTORY scripts/${basename})
    install(DIRECTORY scripts/${basename} DESTINATION libexec)
  else()
    execute_process(COMMAND cat ${CMAKE_SOURCE_DIR}/VERSION OUTPUT_VARIABLE VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND sed -i "s!#####=KRAKEN_DIR=#####!${CMAKE_INSTALL_PREFIX}/libexec!g" scripts/${basename})
    execute_process(COMMAND sed -i "s/#####=VERSION=#####/${VERSION}/g" scripts/${basename})
    install(FILES ${CMAKE_BINARY_DIR}/scripts/${basename} DESTINATION libexec PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE
      GROUP_READ GROUP_EXECUTE
      WORLD_READ WORLD_EXECUTE)
    if(${basename} MATCHES ^krakenuniq.*)
      install(CODE "execute_process( \
        COMMAND ${CMAKE_COMMAND} -E create_symlink \
        ${CMAKE_INSTALL_PREFIX}/libexec/${basename} \
        ${CMAKE_INSTALL_PREFIX}/bin/${basename})")
    endif()
  endif()
endforeach(s)

set(test_exes dump_db_kmers
  grade_classification)
foreach(exe ${test_exes})
  add_executable(${exe} EXCLUDE_FROM_ALL src/${exe}.cpp)
  target_link_libraries(${exe} kraken)
endforeach(exe)

add_custom_target(tests)
add_dependencies(tests grade_classification dump_db_kmers)
