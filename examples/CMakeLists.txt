foreach(example phonon_benchmark custom_excitation)
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} PRIVATE cmpmc::cmpmc cmpmc_build_options)
endforeach()
