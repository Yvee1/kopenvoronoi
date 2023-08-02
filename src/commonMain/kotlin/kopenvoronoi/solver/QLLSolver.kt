package kopenvoronoi.solver

import assert
import kopenvoronoi.geometry.Point
import kopenvoronoi.site.Site
import kopenvoronoi.util.Numeric.chop
import kopenvoronoi.util.Numeric.quadratic_roots
import kopenvoronoi.vertex.Solution

//\brief quadratic-linear-linear Solver
class QLLSolver : Solver() {
    override fun solve(s1: Site, k1: Double, s2: Site, k2: Double, s3: Site, k3: Double, slns: MutableList<Solution>): Int {
        // equation-parameters, in quad-precision
        val quads: MutableList<Eq> = mutableListOf()
        val lins: MutableList<Eq> = mutableListOf()
        val sites: Array<Site> = arrayOf<Site>(s1, s2, s3)
        val kvals = doubleArrayOf(k1, k2, k3)
        for (i in 0..2) {
            val eqn = sites[i].eqp(kvals[i])
            if (sites[i].is_linear()) {
                lins.add(eqn)
            } else {
                quads.add(eqn)
            }
        }
        assert(!quads.isEmpty()) { " !quads.isEmpty() " }
        if (lins.size == 1 || lins.size == 0) {
            assert(quads.size == 3 || quads.size == 2) { " quads.size() == 3 || quads.size() == 2 " }
            for (i in 1 until quads.size) {
                quads[i].subEq(quads[0])
                lins.add(quads[i])
            }
        }
        assert(lins.size == 2) { " lins.size() == 2" }

        // TODO: pick the solution appraoch with the best numerical stability.
        // call all three permutations
        // index shuffling determines if we solve:
        // x and y in terms of t
        // y and t in terms of x
        // t and x in terms of y
        qll_solver(lins, 0, 1, 2, quads[0], k3, slns)
        qll_solver(lins, 2, 0, 1, quads[0], k3, slns)
        qll_solver(lins, 1, 2, 0, quads[0], k3, slns)
        return slns.size
    }

    // \brief qll solver
    // l0 first linear eqn
    // l1 second linear eqn
    // xi,yi,ti indexes to shuffle around
    // xk, yk, kk, rk = params of one ('last') quadratic site (point or arc)
    // solns = output solution triplets (x,y,t) or (u,v,t)
    // returns number of solutions found
    private fun qll_solver(
        lins: List<Eq>,
        xi: Int,
        yi: Int,
        ti: Int,
        quad: Eq,
        k3: Double,
        solns: MutableList<Solution>
    ): Int {
        assert(lins.size == 2) { " lins.size() == 2 " }
        val ai = lins[0].get(xi) // first linear
        val bi = lins[0].get(yi)
        val ki = lins[0].get(ti)
        val ci = lins[0].c
        val aj = lins[1].get(xi) // second linear
        val bj = lins[1].get(yi)
        val kj = lins[1].get(ti)
        val cj = lins[1].c
        val d = chop(ai * bj - aj * bi) // chop! (determinant for 2 linear eqns (?))
        if (d == 0.0) {
            return -1
        }
        // these are the w-equations for qll_solve()
        // (2) u = a1 w + b1
        // (3) v = a2 w + b2
        val a0 = (bi * kj - bj * ki) / d
        val a1 = -(ai * kj - aj * ki) / d
        val b0 = (bi * cj - bj * ci) / d
        val b1 = -(ai * cj - aj * ci) / d
        // based on the 'last' quadratic of (s1,s2,s3)
        val aargs = Array(3) { DoubleArray(2) }
        aargs[0][0] = 1.0
        aargs[0][1] = quad.a
        aargs[1][0] = 1.0
        aargs[1][1] = quad.b
        aargs[2][0] = -1.0
        aargs[2][1] = quad.k
        val isolns = Array(2) { DoubleArray(3) }
        // this solves for w, and returns either 0, 1, or 2 triplets of (u,v,t) in
        // isolns
        // NOTE: indexes of aargs shuffled depending on (xi,yi,ti) !
        val scount = qll_solve(
            aargs[xi][0], aargs[xi][1], aargs[yi][0], aargs[yi][1], aargs[ti][0], aargs[ti][1],
            quad.c,  // xk*xk + yk*yk - rk*rk,
            a0.toDouble(), b0.toDouble(), a1.toDouble(), b1.toDouble(), isolns
        )
        val tsolns = Array(2) { DoubleArray(3) }
        for (i in 0 until scount) {
            tsolns[i][xi] = isolns[i][0] // u x
            tsolns[i][yi] = isolns[i][1] // v y
            tsolns[i][ti] = isolns[i][2] // t t chop!
            solns.add(Solution(Point(tsolns[i][0], tsolns[i][1]), tsolns[i][2], k3))
        }
        // std::cout << " k3="<<kk3<<" qqq_solve found " << scount << " roots\n";
        return scount
    }

    // Solve a system of one quadratic equation, and two linear equations.
    ///
    // (1) a0 u^2 + b0 u + c0 v^2 + d0 v + e0 w^2 + f0 w + g0 = 0
    // (2) u = a1 w + b1
    // (3) v = a2 w + b2
    // solve (1) for w (can have 0, 1, or 2 roots)
    // then substitute into (2) and (3) to find (u, v, t)
    private fun qll_solve(
        a0: Double, b0: Double, c0: Double, d0: Double, e0: Double, f0: Double, g0: Double, a1: Double,
        b1: Double, a2: Double, b2: Double, soln: Array<DoubleArray>
    ): Int {
        // std::cout << "qll_solver()\n";
        // TODO: optimize using abs(a0) == abs(c0) == abs(d0) == 1
        val a = chop(a0 * (a1 * a1) + c0 * (a2 * a2) + e0)
        val b = chop(2 * a0 * a1 * b1 + 2 * a2 * b2 * c0 + a1 * b0 + a2 * d0 + f0)
        val c = a0 * (b1 * b1) + c0 * (b2 * b2) + b0 * b1 + b2 * d0 + g0
        val roots = quadratic_roots(a, b, c) // solves a*w^2 + b*w + c = 0
        return if (roots.isEmpty()) { // No roots, no solutions
            0
        } else {
            for (i in 0 until roots.size) {
                val w: Double = roots.get(i)
                soln[i][0] = a1 * w + b1 // u
                soln[i][1] = a2 * w + b2 // v
                soln[i][2] = w // t
            }
            roots.size
        }
    }
}
