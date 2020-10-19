boundaryToSpatialPolygon<-function(bnd,tag) {
    SpatialPolygons(list(
         Polygons(list(Polygon(bnd[,c("X","Y")])),tag)))
}


