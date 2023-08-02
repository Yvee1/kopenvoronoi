package kopenvoronoi.site

import kopenvoronoi.geometry.Edge
import kopenvoronoi.geometry.Face
import kopenvoronoi.geometry.Point
import kopenvoronoi.offset.LineOfs
import kopenvoronoi.offset.Ofs
import kotlin.math.abs
import kotlin.math.sqrt

//line segment Site
class LineSite : Site {
    var _start // < start Point of LineSite
            : Point
    var _end // < end Point of LineSite
            : Point
    private var _p // < position
            : Point
    var e: Edge? = null // < edge_descriptor to the ::LINESITE pseudo-edge

    constructor(st: Point, en: Point, koff: Double) : this(st, en, koff, null)
    constructor(s: Site) {
        this.eq = s.eqp()
        this.face = s.face
        _start = s.start()
        _end = s.end()
        _p = s.position() // inherit position
    }

    // create line-site between start and end Point.
    constructor(st: Point, en: Point, koff: Double, f: Face?) {
        _start = st
        _end = en
        _p = Point(start().add(end()).mult(0.5))
        if (f != null)
            face = f
        eq.q = false
        eq.a = _end.y - _start.y
        eq.b = _start.x - _end.x
        eq.k = koff // ??
        eq.c = _end.x * _start.y - _start.x * _end.y
        // now normalize
        val d: Double = sqrt(eq.a * eq.a + eq.b * eq.b)
        eq.a /= d
        eq.b /= d
        eq.c /= d
//        assert(abs(eq.a * eq.a + eq.b * eq.b - 1.0) < 1e-5) { " abs( eq.a*eq.a + eq.b*eq.b -1.0 ) < 1e-5" }
    }

    override fun offset(p1: Point, p2: Point): Ofs {
        return LineOfs(p1, p2)
    }

    // closest point on start-end segment to given point.
    // project onto line and return either the projected point
    // or one endpoint of the linesegment
    override fun apex_point(p: Point): Point {
        val s_p = p.sub(_start)
        val s_e = _end.sub(_start)
        val t = s_p.dot(s_e) / s_e.dot(s_e)
        if (t < 0) {
            return _start
        }
        return if (t > 1) {
            _end
        } else {
            _start.add(_end.sub(_start).mult(t))
        }
    }

    override fun in_region(p: Point): Boolean {
        val t = in_region_t(p)
        return t >= 0 && t <= 1
    }

    override fun in_region_t(p: Point): Double {
        val s_p = p.sub(_start)
        val s_e = _end.sub(_start)
        var t = s_p.dot(s_e) / s_e.dot(s_e)
        val eps = 1e-7
        if (abs(t) < eps) {
            t = 0.0
        } else if (abs(t - 1.0) < eps) {
            t = 1.0
        }
        return t
    }

    override fun in_region_t_raw(p: Point): Double {
        val s_p = p.sub(_start)
        val s_e = _end.sub(_start)
        return s_p.dot(s_e) / s_e.dot(s_e)
    }

    override val isLine: Boolean
        get() = true

    override fun position(): Point {
        return _p
    }

    override fun a(): Double {
        return eq.a
    }

    override fun b(): Double {
        return eq.b
    }

    override fun c(): Double {
        return eq.c
    }

    override fun k(): Double {
//        assert(eq.k === 1 || eq.k === -1) { " eq.k==1 || eq.k==-1 " }
        return eq.k
    }

    override fun set_c(p: Point) {
        eq.c = -(eq.a * p.x + eq.b * p.y)
    }

    override fun start(): Point {
        return _start
    }

    override fun end(): Point {
        return _end
    }

    override fun edge(): Edge? {
        return e
    }

    override fun toString(): String {
        return "LS(${_start}>${_end})"
    }
}
