package kopenvoronoi.util

import kopenvoronoi.geometry.Face
import kopenvoronoi.geometry.Point

class KdPoint(p: Point, face: Face) {
    var p: Point
    var face: Face

    init {
        this.p = p
        this.face = face
    }
}
