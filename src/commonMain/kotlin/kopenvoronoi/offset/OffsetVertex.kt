package kopenvoronoi.offset

import kopenvoronoi.geometry.Edge
import kopenvoronoi.geometry.Face
import kopenvoronoi.geometry.Point

//\brief Line- or arc-vertex of an offset curve.
///
//\todo this duplicates the idea of the Ofs class. Remove this or Ofs!
class OffsetVertex {
    var p // < position (start)
            : Point
    var r // < arc radius (line-vertex is indicated by radius of -1)
            : Double
    var c // < arc center
            : Point?
    var cw // < clockwise (or not)
            : Boolean
    var f // < corresponding face in the vd-graph
            : Face?
    var e // corresponding edge
            : Edge

    constructor(p: Point, r: Double, c: Point?, cw: Boolean, f: Face?, e: Edge) {
        this.p = p
        this.r = r
        this.c = c
        this.cw = cw
        this.f = f
        this.e = e
    }

    constructor(p: Point, e: Edge) {
        this.p = p
        r = -1.0
        c = null
        cw = false
        f = null
        this.e = e
    }
}
