foreach(example quickstart custom_excitation phonon_benchmark)
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} PRIVATE cmpmc::cmpmc cmpmc_build_options)
endforeach()
