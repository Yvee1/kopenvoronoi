package kopenvoronoi.site

import kopenvoronoi.geometry.Edge
import kopenvoronoi.geometry.Point
import kopenvoronoi.offset.ArcOfs
import kopenvoronoi.offset.Ofs
import kopenvoronoi.util.Numeric
import kotlin.math.abs

//circular arc Site
class ArcSite(startpt: Point, endpt: Point, centr: Point, dir: Boolean) : Site() {
    var _start // < start Point of arc
            : Point
    var _end // < end Point of arc
            : Point
    var _center // < center Point of arc
            : Point
    var _dir // < CW or CCW direction flag
            : Boolean
    var _radius // < radius of arc
            : Double
    var _k // < offset-direction. +1 for enlarging, -1 for shrinking circle
            : Double
    var e: Edge? = null // < edge_descriptor to ::ARCSITE pseudo-edge

    // create arc-site
    init {
        _start = startpt
        _end = endpt
        _center = centr
        _dir = dir
        _radius = _center.sub(_start).norm()
        eq.q = true
        eq.a = -2 * _center.x
        eq.b = -2 * _center.y
        _k = 1.0
        eq.k = -2 * _k * _radius
        eq.c = _center.x * _center.x + _center.y * _center.y - _radius * _radius
    }

    override fun offset(p1: Point, p2: Point): Ofs {
        return ArcOfs(p1, p2, _center, -1.0)
    } // FIXME: radius

    override fun in_region(p: Point): Boolean {
        if (p === _center) {
            return true
        }
        val t = in_region_t(p)
        return t >= 0 && t <= 1
    }

    // \todo fix arc-site in_region_t test!!
    override fun in_region_t(pt: Point): Double {
        var t = in_region_t_raw(pt) // (diangle_pt - diangle_min) / (diangle_max-diangle_min);
        val eps = 1e-7
        if (abs(t) < eps) {
            t = 0.0
        } else if (abs(t - 1.0) < eps) {
            t = 1.0
        }
        return t
    }

    override fun in_region_t_raw(pt: Point): Double {
        // projection onto circle
        val cen_start = _start.sub(_center)
        val cen_end = _end.sub(_center)
        val cen_pt = pt.sub(_center)
        val proj = _center.add(cen_pt.mult(_radius / cen_pt.norm()))
        val diangle_min: Double
        val diangle_max: Double
        if (!_dir) {
            diangle_min = Numeric.diangle(cen_start.x, cen_start.y)
            diangle_max = Numeric.diangle(cen_end.x, cen_end.y)
        } else {
            diangle_max = Numeric.diangle(cen_start.x, cen_start.y)
            diangle_min = Numeric.diangle(cen_end.x, cen_end.y)
        }
        val diangle_pt = Numeric.diangle(cen_pt.x, cen_pt.y)
        return (diangle_pt - diangle_min) / (diangle_max - diangle_min)
    }

    override fun apex_point(p: Point): Point {
        return if (in_region(p)) {
            projection_point(p)
        } else {
            closer_endpoint(p)
        }
    }

    override fun x(): Double {
        return _center.x
    }

    override fun y(): Double {
        return _center.y
    }

    override fun r(): Double {
        return _radius
    }

    override fun k(): Double {
        return _k
    } // ?

    // return start Point of ArcSite
    override fun start(): Point {
        return _start
    }

    // return end Point of ArcSite
    override fun end(): Point {
        return _end
    }

    // return center Point of ArcSite
    fun center(): Point {
        return _center
    }

    override fun position(): Point {
        return center()
    }

    // return radius of ArcSite
    fun radius(): Double {
        return _radius
    }

    // return true for CW ArcSite and false for CCW
    override fun cw(): Boolean {
        return _dir
    }

    override val isArc: Boolean
        get() = true

    // projection of given Point onto the ArcSite
    private fun projection_point(p: Point): Point {
        return if (p === _center) {
            _start
        } else {
            val dir = p.sub(_center)
            dir.normalize()
            _center.add(dir.mult(_radius)) // this point should lie on the arc
        }
    }

    // return the end Point (either _start or _end) that is closest to the given
    // Point
    private fun closer_endpoint(p: Point): Point {
        val d_start = _start.sub(p).norm()
        val d_end = _end.sub(p).norm()
        return if (d_start < d_end) {
            _start
        } else {
            _end
        }
    }
}
