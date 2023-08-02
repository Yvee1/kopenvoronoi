package kopenvoronoi.solver

import assert
import kopenvoronoi.geometry.Point
import kopenvoronoi.site.Site
import kopenvoronoi.util.Numeric.chop
import kopenvoronoi.util.Numeric.determinant
import kopenvoronoi.vertex.Solution
import kotlin.math.abs

//\brief line-line-line Solver
///
//solves 3x3 system.
class LLLSolver : Solver() {
    //  a1 x + b1 y + c1 + k1 t = 0
    //  a2 x + b2 y + c2 + k2 t = 0
    //  a3 x + b3 y + c3 + k3 t = 0
    //
    // or in matrix form
    //
    //  ( a1 b1 k1 ) ( x )    ( c1 )
    //  ( a2 b2 k2 ) ( y ) = -( c2 )          Ax = b
    //  ( a3 b3 k3 ) ( t )    ( c3 )
    //
    //  Cramers rule x_i = det(A_i)/det(A)
    //  where A_i is A with column i replaced by b
    override fun solve(s1: Site, k1: Double, s2: Site, k2: Double, s3: Site, k3: Double, slns: MutableList<Solution>): Int {
        assert(s1.isLine && s2.isLine && s3.isLine) { " s1.isLine && s2.isLine && s3.isLine " }
        val eq: MutableList<Eq> = mutableListOf() // equation-parameters, in quad-precision
        val sites: Array<Site> = arrayOf<Site>(s1, s2, s3)
        val kvals = doubleArrayOf(k1, k2, k3)
        for (i in 0..2) {
            eq.add(sites[i].eqp(kvals[i]))
        }
        var i = 0
        var j = 1
        val k = 2
        val d = chop(
            determinant(
                eq[i].a, eq[i].b, eq[i].k, eq[j].a, eq[j].b, eq[j].k,
                eq[k].a, eq[k].b, eq[k].k
            )
        )
        val det_eps = 1e-6
        if (abs(d) > det_eps) {
            val t = determinant(
                eq[i].a, eq[i].b, -eq[i].c, eq[j].a, eq[j].b, -eq[j].c,
                eq[k].a, eq[k].b, -eq[k].c
            ) / d
            if (t >= 0) {
                val sol_x = determinant(
                    -eq[i].c, eq[i].b, eq[i].k, -eq[j].c, eq[j].b, eq[j].k,
                    -eq[k].c, eq[k].b, eq[k].k
                ) / d
                val sol_y = determinant(
                    eq[i].a, -eq[i].c, eq[i].k, eq[j].a, -eq[j].c, eq[j].k,
                    eq[k].a, -eq[k].c, eq[k].k
                ) / d
                slns.add(Solution(Point(sol_x, sol_y), t, k3)) // kk3 just passes through without any effect!?
                return 1
            }
        } else {
            // Try parallel solver as fallback, if the small determinant is due to nearly
            // parallel edges
            i = 0
            while (i < 3) {
                j = (i + 1) % 3
                val delta = abs(eq[i].a * eq[j].b - eq[j].a * eq[i].b)
                if (delta <= 1e-300) {
                    val para_solver = LLLPARASolver()
                    val paraSolutions: MutableList<Solution> = mutableListOf()
                    para_solver.solve(
                        sites[i], kvals[i], sites[j], kvals[j], sites[(i + 2) % 3], kvals[(i + 2) % 3],
                        paraSolutions
                    )
                    var solution_count = 0
                    for (s in paraSolutions) {
                        // check that solution has proper offset-direction
                        if (s3.end().sub(s3.start()).cross(s.p.sub(s3.start())) * k3 >= 0) {
                            slns.add(s)
                            solution_count++
                        }
                    }
                    return solution_count
                }
                i++
            }
        }
        return 0 // no solution if determinant zero, or t-value negative
    }
}
