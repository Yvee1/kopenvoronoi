package kopenvoronoi.geometry

import kopenvoronoi.site.Site
import kopenvoronoi.util.Numeric.chop
import kopenvoronoi.util.Numeric.sq
import kopenvoronoi.vertex.Vertex
import kotlin.math.abs
import kotlin.math.sqrt

/*
* bisector formulas
* x = x1 - x2 - x3*t +/- x4 * sqrt( square(x5+x6*t) - square(x7+x8*t) )
* (same formula for y-coordinate)
* line (line/line)
* parabola (circle/line)
* hyperbola (circle/circle)
* ellipse (circle/circle)
*/
class Edge(source: Vertex, target: Vertex) {
    var source: Vertex
    var target: Vertex
    var twin: Edge? = null
    lateinit var next: Edge
    lateinit var face: Face
    var null_face: Face? = null
    var has_null_face: Boolean
    var k = 0.0 // < offset-direction from the adjacent site, either +1 or -1
    var type: EdgeType? = null
    var valid: Boolean
    var x = DoubleArray(8) // < 8-parameter parametrization
    var y = DoubleArray(8) // < 8-parameter parametrization
    var sign = false // < flag to choose either +/- in front of sqrt()
    var inserted_direction = false // < true if ::LINESITE-edge inserted in this direction

    init {
        x[0] = 0.0
        x[1] = 0.0
        x[2] = 0.0
        x[3] = 0.0
        x[4] = 0.0
        x[5] = 0.0
        x[6] = 0.0
        x[7] = 0.0
        y[0] = 0.0
        y[1] = 0.0
        y[2] = 0.0
        y[3] = 0.0
        y[4] = 0.0
        y[5] = 0.0
        y[6] = 0.0
        y[7] = 0.0
        has_null_face = false
        valid = true
        this.source = source
        this.target = target
    }

    fun copyFrom(other: Edge) {
        sign = other.sign
        face = other.face
        k = other.k
        null_face = other.null_face
        has_null_face = other.has_null_face
        type = other.type
        valid = other.valid
        x[0] = other.x[0]
        x[1] = other.x[1]
        x[2] = other.x[2]
        x[3] = other.x[3]
        x[4] = other.x[4]
        x[5] = other.x[5]
        x[6] = other.x[6]
        x[7] = other.x[7]
        y[0] = other.y[0]
        y[1] = other.y[1]
        y[2] = other.y[2]
        y[3] = other.y[3]
        y[4] = other.y[4]
        y[5] = other.y[5]
        y[6] = other.y[6]
        y[7] = other.y[7]
    }

    // \brief return point on edge at given offset-distance t
    ///
    // the eight-parameter formula for a point on the edge is:
    // x = x1 - x2 - x3*t +/- x4 * sqrt( square(x5+x6*t) - square(x7+x8*t) )
    fun point(t: Double): Point {
        val discr1 = chop(sq(x[4] + x[5] * t) - sq(x[6] + x[7] * t), 1e-14)
        val discr2 = chop(sq(y[4] + y[5] * t) - sq(y[6] + y[7] * t), 1e-14)
        return if (discr1 >= 0 && discr2 >= 0) {
            val psig = (if (sign) +1 else -1).toDouble()
            val nsig = (if (sign) -1 else +1).toDouble()
            val xc: Double = x[0] - x[1] - x[2] * t + psig * x[3] * sqrt(discr1)
            val yc: Double = y[0] - y[1] - y[2] * t + nsig * y[3] * sqrt(discr2)
            if (xc != xc) { // test for NaN!
                throw RuntimeException()
            }
            Point(xc, yc)
        } else {
            Point(x[0] - x[1] - x[2] * t, y[0] - y[1] - y[2] * t) // coordinates without sqrt()
        }
    }

    /**
     * Returns the midpoint of the edge source and target positions
     *
     * @return
     */
    fun position(): Point {
        return source.position.add(target.position).mult(0.5) // TODO pre-calculate?
    }

    // dispatch to setter functions based on type of \a s1 and \a s2
    fun set_parameters(s1: Site, s2: Site, sig: Boolean) {
        sign = sig // sqrt() sign for edge-parametrization
        if (s1.isPoint && s2.isPoint) {
            set_pp_parameters(s1, s2)
        } else if (s1.isPoint && s2.isLine) {
            set_pl_parameters(s1, s2)
        } else if (s2.isPoint && s1.isLine) { // LP
            set_pl_parameters(s2, s1)
            sign = !sign
        } else if (s1.isLine && s2.isLine) {
            set_ll_parameters(s2, s1)
        } else if (s1.isPoint && s2.isArc) {
            set_pa_parameters(s1, s2)
        } else if (s2.isPoint && s1.isArc) { // AP
            sign = !sign
            set_pa_parameters(s2, s1)
        } else if (s1.isLine && s2.isArc) {
            set_la_parameters(s1, s2)
        } else if (s2.isLine && s1.isArc) {
            set_la_parameters(s2, s1)
        } else {
            throw RuntimeException("Unexpected combination of sites")
            // AA
        }
    }

    // set edge parameters for PointSite-PointSite edge
    fun set_pp_parameters(s1: Site, s2: Site) {
//        assert(s1.isPoint && s2.isPoint) { " s1.isPoint && s2.isPoint " }
        val d = s1.position().sub(s2.position()).norm()
        val alfa1 = (s2.x() - s1.x()) / d
        val alfa2 = (s2.y() - s1.y()) / d
        val alfa3 = -d / 2
        type = EdgeType.LINE
        x[0] = s1.x()
        x[1] = (alfa1 * alfa3).toDouble() //
        x[2] = 0.0
        x[3] = -alfa2
        x[4] = 0.0
        x[5] = +1.0
        x[6] = alfa3.toDouble()
        x[7] = 0.0
        y[0] = s1.y()
        y[1] = (alfa2 * alfa3).toDouble()
        y[2] = 0.0
        y[3] = -alfa1
        y[4] = 0.0
        y[5] = +1.0
        y[6] = alfa3.toDouble()
        y[7] = 0.0
    }

    // set ::PARABOLA edge parameters (between PointSite and LineSite).
    fun set_pl_parameters(s1: Site, s2: Site) {
//        assert(s1.isPoint && s2.isLine) { " s1.isPoint && s2.isLine " }
        type = EdgeType.PARABOLA
        val alfa3 = s2.a() * s1.x() + s2.b() * s1.y() + s2.c() // signed distance to line

        // figure out kk, i.e. offset-direction for LineSite
        x[0] = s1.x() // xc1
        x[1] = s2.a() * alfa3 // alfa1*alfa3
        x[2] = s2.a() // *kk; // -alfa1 = - a2 * k2?
        x[3] = s2.b() // alfa2 = b2
        x[4] = 0.0 // alfa4 = r1 (PointSite has zero radius)
        x[5] = +1.0 // lambda1 (allways positive offset from PointSite)
        x[6] = alfa3 // alfa3= a2*xc1+b2*yc1+d2?
        x[7] = +1.0 // kk; // -1 = k2 side of line??
        y[0] = s1.y() // yc1
        y[1] = s2.b() * alfa3 // alfa2*alfa3
        y[2] = s2.b() // *kk; // -alfa2 = -b2
        y[3] = s2.a() // alfa1 = a2
        y[4] = 0.0 // alfa4 = r1 (PointSite has zero radius)
        y[5] = +1.0 // lambda1 (allways positive offset from PointSite)
        y[6] = alfa3 // alfa3
        y[7] = +1.0 // kk; // -1 = k2 side of line??
    }

    // set ::SEPARATOR edge parameters
    fun set_sep_parameters(endp: Point, p: Point) {
        type = EdgeType.SEPARATOR
        val dx = p.x - endp.x
        val dy = p.y - endp.y
        val d = p.sub(endp).norm()
        x[0] = endp.x
        x[2] = (-dx / d).toDouble() // negative of normalized direction from endp to p
        y[0] = endp.y
        y[2] = (-dy / d).toDouble()
        x[1] = 0.0
        x[3] = 0.0
        x[4] = 0.0
        x[5] = 0.0
        x[6] = 0.0
        x[7] = 0.0
        y[1] = 0.0
        y[3] = 0.0
        y[4] = 0.0
        y[5] = 0.0
        y[6] = 0.0
        y[7] = 0.0
    }

    // set edge parametrization for LineSite-LineSite edge (parallel case)
    fun set_ll_para_parameters(s1: Site, s2: Site) {
//        assert(s1.isLine && s2.isLine) { " s1.isLine && s2.isLine " }
        type = EdgeType.PARA_LINELINE

        // find a point (x1,y1) on the line s1
        // ax+by+c=0
        var x1 = 0.0
        var y1 = 0.0
        if (abs(s1.a()) > abs(s1.b())) {
            y1 = 0.0
            x1 = -s1.c() / s1.a()
        } else {
            x1 = 0.0
            y1 = -s1.c() / s1.b()
        }

        // find a point (x2,y2) on the line s2
        // ax+by+c=0
        var x2 = 0.0
        var y2 = 0.0
        if (abs(s2.a()) > abs(s2.b())) {
            y2 = 0.0
            x2 = -s2.c() / s2.a()
        } else {
            x2 = 0.0
            y2 = -s2.c() / s2.b()
        }

        // now e.g. the s2 line is given by
        // p = (x2,y2) + t*(-b2, a)
        // and we can find the projection of (x1,y1) onto s2 as
        // p1 = p2 = p0 + t*v
        val p1 = Point(x1, y1)
        val p2 = Point(x2, y2)
        val v = Point(-s2.b(), s2.a())
        val t = p1.sub(p2).dot(v) / v.dot(v)
        val p1_proj = p2.add(v.mult(t))
//        assert(p1.sub(p1_proj).norm() > 0) { " p1.sub(p1_proj).norm() > 0 " }

        // from this point, go a distance d/2 in the direction of the normal
        // to find a point through which the bisector passes
        x1 = x1 + p1_proj.sub(p1).x / 2
        y1 = y1 + p1_proj.sub(p1).y / 2
        // the tangent of the bisector (as well as the two line-sites) is a vector
        // (-b , a)
        x[0] = x1
        x[1] = -s1.b()
        y[0] = y1
        y[1] = s1.a()
        x[2] = 0.0
        x[3] = 0.0
        x[4] = 0.0
        x[5] = 0.0
        x[6] = 0.0
        x[7] = 0.0
        y[2] = 0.0
        y[3] = 0.0
        y[4] = 0.0
        y[5] = 0.0
        y[6] = 0.0
        y[7] = 0.0
    }

    // set edge parametrization for LineSite-LineSite edge
    fun set_ll_parameters(s1: Site, s2: Site) { // Held thesis p96
//        assert(s1.isLine && s2.isLine) { " s1.isLine && s2.isLine " }
        type = EdgeType.LINELINE
        val delta = s1.a() * s2.b() - s1.b() * s2.a()

        // (numerically) parallel line segments - the generic LLL solver
        // is numerically unstable for parallel cases
        if (abs(delta) <= 1e-14) {
            set_ll_para_parameters(s1, s2)
            return
        }
//        assert(delta != 0) { " delta != 0 " }
        val alfa1 = (s1.b() * s2.c() - s2.b() * s1.c()) / delta
        val alfa2 = (s2.a() * s1.c() - s1.a() * s2.c()) / delta
        val alfa3 = -(s2.b() - s1.b()) / delta
        val alfa4 = -(s1.a() - s2.a()) / delta

        // point (alfa1,alfa2) is the intersection point between the line-segments
        // vector (-alfa3,-alfa4) is the direction/tangent of the bisector
        x[0] = alfa1
        x[2] = -alfa3
        y[0] = alfa2
        y[2] = -alfa4
        x[1] = 0.0
        x[3] = 0.0
        x[4] = 0.0
        x[5] = 0.0
        x[6] = 0.0
        x[7] = 0.0
        y[1] = 0.0
        y[3] = 0.0
        y[4] = 0.0
        y[5] = 0.0
        y[6] = 0.0
        y[7] = 0.0
    }

    // set edge parameters when s1 is PointSite and s2 is ArcSite
    fun set_pa_parameters(s1: Site, s2: Site) {
//        assert(s1.isPoint && s2.isArc) { " s1.isPoint && s2.isArc " }
        // std::cout << "set_pa_parameters()\n";
        type = EdgeType.HYPERBOLA // hyperbola or ellipse?
        var lamb2 = 1.0

        // distance between centers
        val d: Double =
            sqrt((s1.x() - s2.x()) * (s1.x() - s2.x()) + (s1.y() - s2.y()) * (s1.y() - s2.y()))
//        assert(d > 0) { " d > 0 " }
        if (d <= s2.r()) {
            lamb2 = -1.0
            sign = !sign
        }
        val alfa1 = (s2.x() - s1.x()) / d
        val alfa2 = (s2.y() - s1.y()) / d
        val alfa3 = (s2.r() * s2.r() - d * d) / (2 * d)
        val alfa4 = lamb2 * s2.r() / d
        x[0] = s1.x()
        x[1] = (alfa1 * alfa3).toDouble()
        x[2] = (alfa1 * alfa4).toDouble()
        x[3] = alfa2
        x[4] = 0.0 // r1; PointSite has zero radius
        x[5] = +1.0 // lamb1; allways outward offset from PointSite
        x[6] = alfa3
        x[7] = alfa4
        y[0] = s1.y()
        y[1] = (alfa2 * alfa3).toDouble()
        y[2] = (alfa2 * alfa4).toDouble()
        y[3] = alfa1
        y[4] = 0.0 // r1; PointSite has zero radius
        y[5] = +1.0 // lamb1; allways outward offset from PointSite
        y[6] = alfa3
        y[7] = alfa4
    }

    // set edge parameters when s1 is ArcSite and s2 is LineSite
    fun set_la_parameters(s1: Site, s2: Site) {
//        assert(s1.isLine && s2.isArc) { " s1.isLine && s2.isArc " }
        type = EdgeType.PARABOLA
        val lamb2: Double
        lamb2 = if (s2.cw()) {
            +1.0
        } else {
            -1.0
        }
        val alfa1 = s1.a() // a2
        val alfa2 = s1.b() // b2
        val alfa3 = s1.a() * s2.x() + s1.b() * s2.y() + s1.c()
        val alfa4 = s2.r()
        val kk = +1.0 // # positive line-offset
        // figure out sign?
        x[0] = s2.x()
        x[1] = (alfa1 * alfa3).toDouble()
        x[2] = alfa1 * kk
        x[3] = alfa2
        x[4] = alfa4
        x[5] = lamb2
        x[6] = alfa3
        x[7] = kk
        y[0] = s2.y()
        y[1] = (alfa2 * alfa3).toDouble()
        y[2] = alfa2 * kk
        y[3] = alfa1
        y[4] = alfa4
        y[5] = lamb2
        y[6] = alfa3
        y[7] = kk
    }

    // \return minumum t-value for this edge
    // this function dispatches to a helper-function based on the Site:s \a s1 and
    // \a s2
    // used only for positioning APEX vertices?
    fun minimum_t(s1: Site, s2: Site): Double {
        return if (s1.isPoint && s2.isPoint) {
            minimum_pp_t(s1, s2)
        } else if (s1.isPoint && s2.isLine) {
            minimum_pl_t(s1, s2)
        } else if (s2.isPoint && s1.isLine) {
            minimum_pl_t(s2, s1)
        } else if (s1.isLine && s2.isLine) {
            0.0
        } else if (s1.isPoint && s2.isArc) {
            minimum_pa_t(s1, s2)
        } else if (s2.isPoint && s1.isArc) {
            minimum_pa_t(s2, s1)
        } else {
            throw RuntimeException("Unexpected site types")
            // todo: AP, AL, AA
        }
    }

    // minimum t-value for LINE edge between PointSite and PointSite
    fun minimum_pp_t(s1: Site, s2: Site): Double {
//        assert(s1.isPoint && s2.isPoint) { " s1.isPoint && s2.isPoint " }
        val p1p2 = s1.position().sub(s2.position()).norm()
//        assert(p1p2 >= 0) { " p1p2 >=0 " }
        return p1p2 / 2 // this splits point-point edges at APEX
    }

    // minimum t-value for ::PARABOLA edge
    fun minimum_pl_t(s1: Site?, s2: Site?): Double {
        val mint = -x[6] / (2.0 * x[7])
//        assert(mint >= 0) { " mint >=0 " }
        return mint
    }

    // minimum t-value for edge between PointSite and ArcSite
    fun minimum_pa_t(s1: Site, s2: Site): Double {
//        assert(s1.isPoint && s2.isArc) { " s1.isPoint && s2.isArc " }
        val p1p2 = s1.position().sub(s2.apex_point(s1.position())).norm() // - s2->r() ;
//        assert(p1p2 >= 0) { " p1p2 >=0 " }
        return p1p2 / 2 // this splits point-point edges at APEX
    }

    override fun toString(): String {
        return "E(${source.position}>${target.position})"
    }
}
