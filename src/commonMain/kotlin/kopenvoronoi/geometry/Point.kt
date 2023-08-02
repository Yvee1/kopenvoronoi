package kopenvoronoi.geometry

import kotlin.math.sqrt

/**
 * A point or vector in 2D with coordinates (x, y)
 */
class Point {
    var x: Double
    var y: Double

    constructor() {
        x = 0.0
        y = 0.0
    }

    constructor(xi: Double, yi: Double) {
        x = xi
        y = yi
    }

    constructor(p: Point) {
        x = p.x
        y = p.y
    }

    fun dot(p: Point): Double {
        return x * p.x + y * p.y
    }

    fun cross(p: Point): Double {
        return x * p.y - y * p.x
    }

    fun norm(): Double {
        return sqrt(x * x + y * y)
    }

    fun norm_sq(): Double {
        return x * x + y * y
    }

    fun normalize() {
        if (norm() != 0.0) {
            multEq(1 / norm())
        }
    }

    // return perpendicular in the xy plane, rotated 90 degree to the left
    fun xy_perp(): Point {
        return Point(-y, x)
        // 2D rotation matrix:
        // cos -sin
        // sin cos
        // for theta = 90
        // 0 -1 ( x )
        // 1 0 ( y ) = ( -y x )
    }

    // is this Point right of line through points \a p1 and \a p2 ?
    fun is_right(p1: Point, p2: Point): Boolean {
        // this is an ugly way of doing a determinant
        // should be prettyfied sometime...
        // \todo FIXME: what if p1==p2 ? (in the XY plane)
        val a1 = p2.x - p1.x
        val a2 = p2.y - p1.y
        val t2 = -a1
        val b1 = x - p1.x
        val b2 = y - p1.y
        val t = a2 * b1 + t2 * b2
        return if (t > 0.0) {
            true
        } else {
            false
        }
    }

    fun addEq(p: Point) {
        x += p.x
        y += p.y
    }

    fun subEq(p: Point) {
        x -= p.x
        y -= p.y
    }

    fun add(p: Point): Point {
        return Point(x + p.x, y + p.y)
    }

    fun sub(p: Point): Point {
        return Point(x - p.x, y - p.y)
    }

    fun distance(p: Point): Double {
        val dx = x - p.x
        val dy = y - p.y
        return sqrt(dx * dx + dy * dy)
    }

    fun multEq(a: Double) {
        x *= a
        y *= a
    }

    fun mult(a: Double): Point {
        return Point(x * a, y * a)
    }

    override fun equals(other: Any?): Boolean {
        return this === other || x == (other as Point?)!!.x && y == other!!.y
    }

    override fun hashCode(): Int {
        TODO()
//        return Double.valueOf(x).hashCode() * 42 + Double.valueOf(y).hashCode()
    }

    override fun toString(): String {
//        return String.format("%.3f,%.3f", x, y)
        return "${x},${y}"
    }
}
