set(ITKIOBioRad_SRC
itkBioRadImageIOFactory.cxx
itkBioRadImageIO.cxx
)

add_library(ITKIOBioRad ${ITKIOBioRad_SRC})
target_link_libraries(ITKIOBioRad  ${ITKIOImageBase_LIBRARIES})
itk_module_target(ITKIOBioRad)
