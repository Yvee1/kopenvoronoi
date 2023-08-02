package kopenvoronoi.solver

import kopenvoronoi.geometry.Point
import kopenvoronoi.site.Site
import kopenvoronoi.vertex.Solution
import kotlin.math.abs

//\brief line-line-line Solver (parallel line-segment case)
///
//solves 3x3 system.
class LLLPARASolver : Solver() {
    // parallel linesegment edge case.
    //  a1 x + b1 y + c1 + k1 t = 0
    //  a2 x + b2 y + c2 + k2 t = 0
    //  a3 x + b3 y + c3 + k3 t = 0
    //
    // s1 and s2 are parallel, so they have a PARA_LINELINE edge between them
    //
    // this constrains the solution to lie on a line parallel to s1/s2
    // passing through a point equidistant from s1/s2
    //
    // equation of bisector is:
    // ab x + bb y + cb = 0
    // ab = a1
    // bb = b1
    // cb = (c1+c2)2
    // all points on the bisector have a t value
    // tb = fabs(c1-c2)/2
    //
    // find intersection of bisector and offset of third site
    //  ab x + bb y + cb = 0
    //  a3 x + b3 y + c3 + k3 tb = 0
    //  or
    //  ( ab  bb ) ( x ) = ( -cb )
    //  ( a3  b3 ) ( y ) = ( -c3-k3*tb )
    //
    override fun solve(s1: Site, k1: Double, s2: Site, k2: Double, s3: Site, k3: Double, slns: MutableList<Solution>): Int {
//        assert(s1.isLine && s2.isLine && s3.isLine) { " s1.isLine && s2.isLine && s3.isLine " }
        val bisector = Eq()
        bisector.a = s1.a()
        bisector.b = s1.b()
        var s2c = s2.c()

        // if s1 and s2 have opposite (a,b) normals, flip the sign of s2c
        val n0 = Point(s1.a(), s1.b())
        val n1 = Point(s2.a(), s2.b())
        if (n0.dot(n1) < 0) {
            s2c = -s2c
        }
        bisector.c = (s1.c() + s2c) * 0.5
        val tb = 0.5 * abs(s1.c() - s2c) // bisector offset distance
        val xy = two_by_two_solver(bisector.a, bisector.b, s3.a(), s3.b(), -bisector.c, -s3.c() - k3 * tb)
        return if (xy != null) {
            slns.add(Solution(Point(xy.first, xy.second), tb, k3))
            1
        } else {
            0
        }
    }

    // solve 2z2 system Ax = y by inverting A
    // x = Ainv * y
    // returns false if det(A)==0, i.e. no solution found
    fun two_by_two_solver(a: Double, b: Double, c: Double, d: Double, e: Double, f: Double): Pair<Double, Double>? {
        // [ a b ] [u] = [ e ]
        // [ c d ] [v] = [ f ]
        // matrix inverse is
        // [ d -b ]
        // 1/det * [ -c a ]
        // so
        // [u] [ d -b ] [ e ]
        // [v] = 1/det * [ -c a ] [ f ]
        val det = a * d - c * b
        if (abs(det) < 1e-15) {
            return null
        }
        val u = 1.0 / det * (d * e - b * f)
        val v = 1.0 / det * (-c * e + a * f)
        return Pair(u, v)
    }
}
