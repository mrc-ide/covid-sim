set(GEOMETRY_SRC_FILES Vector2.cpp)
set(GEOMETRY_HDR_FILES Size.h Vector2.h BoundingBox.h)
source_group(covidsim\\geometry FILES ${GEOMETRY_SRC_FILES} ${GEOMETRY_HDR_FILES})

add_library(geometrylib ${GEOMETRY_SRC_FILES} ${GEOMETRY_HDR_FILES})
