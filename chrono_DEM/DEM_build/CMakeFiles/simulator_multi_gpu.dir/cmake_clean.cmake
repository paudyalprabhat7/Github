file(REMOVE_RECURSE
  "libsimulator_multi_gpu.a"
  "libsimulator_multi_gpu.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CUDA CXX)
  include(CMakeFiles/simulator_multi_gpu.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
