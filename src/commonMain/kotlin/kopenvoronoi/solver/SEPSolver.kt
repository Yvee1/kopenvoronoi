package kopenvoronoi.solver

import kopenvoronoi.geometry.Point
import kopenvoronoi.site.Site
import kopenvoronoi.vertex.Solution

// this solver is called when we want to position a vertex on a SEPARATOR edge
// a SEPARATOR edge exists between a LineSite and one of its PointSite end-points
// the input sites are thus s1=LineSite and s2=PointSite  (if arcs are supported in the future then s1=ArcSite is possible)
// s3 can be either a LineSite or a PointSite (of arcs are supported, s3=ArcSite is possible)
//
//  s1 (LineSite) offset eq. is     a1 x + b1 y + c1 + k1 t = 0   (1)
//  s2 (PointSite) offset eq. is    (x-x2)^2 + (y-y2)^2 = t^2     (2)
//
// Two possibilities for s3:
//   s3 (LineSite)   a3 x + b3 y + c3 + k3 t = 0
//   s3 (PointSite)  (x-x3)^2 + (y-y3)^2 = t^2
//
// This configuration constrains the solution to lie on the separator edge.
// The separator is given by
// SEP = p2 + t* sv
// where p2 is the location of s2, and the separator direction sv is
// sv = (-a1,-b1)   if k1=-1
// sv = (a1,b1)   if k1=-1
// thus points on the separator are located at:
//
//  x_sep = x2 + t*sv.x
//  y_sep = y2 + t*sv.y
//
//  This can be inserted into (1) or (2) above, which leads to a linear equation in t.
//
//  Insert into (1):
//    a3 (x2 + t*sv.x) + b3 (y2 + t*sv.y) + c3 + k3 t = 0
//      ==>
//          t = -( a3*x2 + b3*y2 + c3 ) / (sv.x*a3 + sv.y*b3 + k3)
//
//  Insert into (2):
//    (x2 + t*sv.x-x3)^2 + (y2 + t*sv.y-y3)^2 = t^2
//    ==> (using dx= x2-x3 and dy = x2-x3)
//   t^2 (sv.x^2 + sv.y^1 - 1)  + t (2*dx*sv.x + 2*dy*sv.y) + dx^2 + dy^2 = 0
//    ==>  (since sv is a unit-vector sv.x^2 + sv.y^1 - 1 = 0)
//         t = - (dx^2+dy^2) / (2*(dx*sv.x + dy*sv.y))
//
//  FIXME: what happens if we get a divide by zero situation ??
//
//\brief ::SEPARATOR Solver
class SEPSolver : Solver() {
    override fun solve(s1: Site, k1: Double, s2: Site, k2: Double, s3: Site, k3: Double, slns: MutableList<Solution>): Int {
//        assert(s1.isLine && s2.isPoint) { " s1.isLine && s2.isPoint " }
//        assert(s3.isLine) { "s3.isLine" }

        // separator direction
        val sv = Point(-s1.a(), -s1.b())
        val tsln = -(s3.a() * s2.x() + s3.b() * s2.y() + s3.c()) / (sv.x * s3.a() + sv.y * s3.b() + k3)
        val psln = Point(s2.x(), s2.y()).add(sv.mult(tsln))
        slns.add(Solution(psln, tsln, k3))
        return 1
    }
}
