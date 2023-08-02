package kopenvoronoi.site

import kopenvoronoi.geometry.Face
import kopenvoronoi.geometry.Point
import kopenvoronoi.offset.ArcOfs
import kopenvoronoi.offset.Ofs
import kopenvoronoi.vertex.Vertex

//vertex Site
class PointSite : Site {
    private var _p // < position
            : Point
    var v: Vertex? = null // < vertex descriptor of this PointSite

    constructor(p: Point) {
        _p = p
//        face = null
        eq.q = true
        eq.a = -2 * p.x
        eq.b = -2 * p.y
        eq.k = 0.0
        eq.c = p.x * p.x + p.y * p.y
    }

    // ctor
    constructor(p: Point, f: Face) {
        _p = p
        face = f
        eq.q = true
        eq.a = -2 * p.x
        eq.b = -2 * p.y
        eq.k = 0.0
        eq.c = p.x * p.x + p.y * p.y
    }

    // ctor
    constructor(p: Point, f: Face, vert: Vertex?) {
        v = vert
        _p = p
        face = f
        eq.q = true
        eq.a = -2 * p.x
        eq.b = -2 * p.y
        eq.k = 0.0
        eq.c = p.x * p.x + p.y * p.y
    }

    override fun apex_point(p: Point): Point {
        return _p
    }

    override fun offset(p1: Point, p2: Point): Ofs {
        val rad = p1.sub(_p).norm()
        return ArcOfs(p1, p2, _p, rad)
    }

    override fun position(): Point {
        return _p
    }

    override fun x(): Double {
        return _p.x
    }

    override fun y(): Double {
        return _p.y
    }

    override fun r(): Double {
        return 0.0
    }

    override fun k(): Double {
        return 0.0
    }

    override val isPoint: Boolean
        get() = true

    override fun in_region(p: Point): Boolean {
        return true
    }

    override fun in_region_t(p: Point): Double {
        return (-1).toDouble()
    }

    override fun vertex(): Vertex? {
        return v
    }

    override fun toString(): String {
        return "PS(${_p})"
    }
}
