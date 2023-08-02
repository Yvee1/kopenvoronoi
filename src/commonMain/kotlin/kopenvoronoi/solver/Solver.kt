package kopenvoronoi.solver

import kopenvoronoi.site.Site
import kopenvoronoi.vertex.Solution

/**
 * Abstract base-class for voronoi vertex position solvers
 *
 * The input to the solver is three Sites (s1,s2,s3) and three offset-directions
 * (k1,k2,k3). The optput is a vector with one or more Solution.
 */
abstract class Solver {
    // \brief solve for position of VoronoiVertex with given adjacent sites and
    // directions
    ///
    // \param s1 first adjacent Site
    // \param k1 direction from \a s1 to new VoronoiVertex
    // \param s2 second adjacent Site
    // \param k2 direction from \a s2 to new VoronoiVertex
    // \param s3 third adjacent Site
    // \param k3 direction from \a s3 to new VoronoiVertex
    // \param slns Solution vector, will be updated by Solver
    abstract fun solve(s1: Site, k1: Double, s2: Site, k2: Double, s3: Site, k3: Double, slns: MutableList<Solution>): Int

    // used by alt_sep_solver
//    fun set_type(t: Int) {
//        type = t
//    }

    // set the debug mode to \a b
//    fun set_debug(b: Boolean) {
//        debug = b
//    }

    // no warnings/messages to stdout will be written, if silent is set true.
//    fun set_silent(b: Boolean) {
//        silent = b
//    }

    // flag for debug output
    var debug = false

    // separator case type.
    // - type = 0 means l3 / p1 form a separator
    // - type = 1 means l3 / p2 form a separator
    var type = 0
    var silent = false // < suppress all warnings or other stdout output
}
