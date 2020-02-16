cmake_minimum_required(VERSION 3.5)

project(icyFlow2 LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -qopenmp")



find_package(VTK COMPONENTS

  vtkCommonColor
  vtkCommonCore
  vtkFiltersSources
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  vtkViewsQt
  QUIET
)

if (NOT VTK_FOUND)
  message("Skipping RenderWindowUISingleInheritance: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")

if(vtkGUISupportQt_LOADED OR TARGET VTK::GUISupportQt)
  if (VTK_VERSION VERSION_LESS "8.90.0")
    # Instruct CMake to run moc automatically when needed.
    set(CMAKE_AUTOMOC ON)
  else()
     # Instruct CMake to run moc automatically when needed.
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)
  endif()
  # We have ui files, this will bring in the macro: qt5_wrap_ui
  find_package(Qt5Widgets REQUIRED QUIET)
else()
  message(STATUS "GuiSupportQt not found.")
  return()
endif()

include(${VTK_USE_FILE})

set(qtpb_dir /home/s2/Qt/5.12.6/Src/qttools/src/shared/qtpropertybrowser)
SET(qtpropertyeditor_SRCS
${qtpb_dir}/qtpropertybrowser.cpp
${qtpb_dir}/qtpropertymanager.cpp
${qtpb_dir}/qteditorfactory.cpp
${qtpb_dir}/qtvariantproperty.cpp
${qtpb_dir}/qttreepropertybrowser.cpp
${qtpb_dir}/qtbuttonpropertybrowser.cpp
${qtpb_dir}/qtgroupboxpropertybrowser.cpp
${qtpb_dir}/qtpropertybrowserutils.cpp
${qtpb_dir}/qtpropertybrowser.h
${qtpb_dir}/qtpropertymanager.h
${qtpb_dir}/qteditorfactory.h
${qtpb_dir}/qtvariantproperty.h
${qtpb_dir}/qttreepropertybrowser.h
${qtpb_dir}/qtbuttonpropertybrowser.h
${qtpb_dir}/qtgroupboxpropertybrowser.h
${qtpb_dir}/qtpropertybrowserutils_p.h
${qtpb_dir}/qtpropertybrowser.qrc
)


find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt5Charts)

  add_executable(icyFlow2
    main.cpp
    mainwindow.cpp
    mainwindow.h
    mainwindow.ui
    ${qtpropertyeditor_SRCS}
    pb/objectpropertybrowser.cpp
    pb/objectpropertybrowser.h
    generator/beamparams.cpp
    generator/beamparams.h
    generator/generatortool.cpp
    generator/generatortool.h
    generator/czinsertiontool.cpp
    generator/czinsertiontool.h
    geometry/node.cpp
    geometry/node.h
    geometry/face.h
    geometry/face.cpp
    geometry/element.cpp
    geometry/element.h
    geometry/cz.h
    geometry/cz.cpp
    geometry/mesh.cpp
    geometry/mesh.h
    geometry/meshcollection.cpp
    geometry/meshcollection.h
    simulation/modelprms.h
    simulation/csrdictionary.cpp
    simulation/csrdictionary.h
    simulation/frameinfo.h
    simulation/implicitmodel4.cpp
    simulation/implicitmodel4.h
    simulation/linearsystem.cpp
    simulation/linearsystem.h
    simulation/bvh/bvht.cpp
    simulation/bvh/bvht.h
    simulation/bvh/bvhn.cpp
    simulation/bvh/bvhn.h
    simulation/bvh/kdop24.cpp
    simulation/bvh/kdop24.h
    simulation/point_triangle_distance.cpp
    simulation/GteSymmetricEigensolver3x3.h
    simulation/numbercrunching.cpp
    simulation/numbercrunching.h
    SimpleObjectPool.h
    backgroundworker.cpp
    backgroundworker.h
  )

target_include_directories(icyFlow2 PUBLIC ${qtpb_dir})
target_include_directories(icyFlow2 PUBLIC generator pb)
target_include_directories(icyFlow2 PUBLIC /home/s2/intel/mkl/include)
target_link_directories(icyFlow2 PUBLIC /home/s2/intel/mkl/lib/intel64_lin)

target_link_libraries(icyFlow2 PRIVATE Qt5::Widgets Qt5::Charts)
target_link_libraries(icyFlow2 PRIVATE ${VTK_LIBRARIES})
target_link_libraries(icyFlow2 PRIVATE gmsh)
target_link_libraries(icyFlow2 PRIVATE mkl_rt)
