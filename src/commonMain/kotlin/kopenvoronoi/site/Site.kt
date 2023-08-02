package kopenvoronoi.site

import kopenvoronoi.geometry.Edge
import kopenvoronoi.geometry.Face
import kopenvoronoi.geometry.Point
import kopenvoronoi.offset.Ofs
import kopenvoronoi.solver.Eq
import kopenvoronoi.vertex.Vertex

//Base-class for a voronoi-diagram site, or generator.
abstract class Site protected constructor() {
    // the HEFace of this Site
    lateinit var face: Face

    // equation parameters
    var eq: Eq = Eq()

    // return closest point on site to given point p
    abstract fun apex_point(p: Point): Point

    // return offset of site
    abstract fun offset(p1: Point, p2: Point): Ofs

    // position of site for PointSite
    abstract fun position(): Point

    // start point of site (for LineSite and ArcSite)
    open fun start(): Point {
        throw UnsupportedOperationException()
    }

    // end point of site (for LineSite and ArcSite)
    open fun end(): Point {
        throw UnsupportedOperationException()
    }

    // return equation parameters
    fun eqp(): Eq {
        return eq
    }

    // return equation parameters
    fun eqp(kk: Double): Eq {
        val eq2 = Eq(eq)
        eq2.k *= kk
        return eq2
    }

    // true for LineSite
    fun is_linear(): Boolean {
        return isLine
    }

    // true for PointSite and ArcSite
    fun is_quadratic(): Boolean {
        return isPoint
    }

    // x position
    open fun x(): Double {
        throw UnsupportedOperationException()
    }

    // y position
    open fun y(): Double {
        throw UnsupportedOperationException()
    }

    // radius (zero for PointSite)
    open fun r(): Double {
        throw UnsupportedOperationException()
    }

    // offset direction
    open fun k(): Double {
        throw UnsupportedOperationException()
    }

    // LineSite a parameter
    open fun a(): Double {
        throw UnsupportedOperationException()
    }

    // LineSite b parameter
    open fun b(): Double {
        throw UnsupportedOperationException()
    }

    // LineSite c parameter
    open fun c(): Double {
        throw UnsupportedOperationException()
    }

    open fun set_c(p: Point) {
        throw UnsupportedOperationException()
    }

    open val isPoint: Boolean
        // true for PointSite
        get() = false
    open val isLine: Boolean
        // true for LineSite
        get() = false
    open val isArc: Boolean
        // true for ArcSite
        get() = false

    // true for CW oriented ArcSite
    open fun cw(): Boolean {
        return false
    }

    // is given Point in_region ?
    abstract fun in_region(p: Point): Boolean

    // is given Point in region?
    open fun in_region_t(p: Point): Double {
        throw UnsupportedOperationException()
    }

    // in-region t-valye
    open fun in_region_t_raw(p: Point): Double {
        throw UnsupportedOperationException()
    }

    // return edge (if this is a LineSite or ArcSite
    open fun edge(): Edge? {
        throw UnsupportedOperationException()
    }

    // return vertex, if this is a PointSite
    open fun vertex(): Vertex? {
        throw UnsupportedOperationException()
    }
}
