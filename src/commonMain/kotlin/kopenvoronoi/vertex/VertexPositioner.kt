package kopenvoronoi.vertex

import assert
import kopenvoronoi.HalfEdgeDiagram
import kopenvoronoi.geometry.Edge
import kopenvoronoi.geometry.EdgeType
import kopenvoronoi.geometry.Point
import kopenvoronoi.site.Site
import kopenvoronoi.solver.*
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min
import koptimize.minimize

//Calculates the (x,y) position of a VoronoiVertex in the VoronoiDiagram
class VertexPositioner(gi: HalfEdgeDiagram) {
    // solvers, to which we dispatch, depending on the input sites
    var ppp_solver // < point-point-point solver
            : Solver
    var lll_solver // < line-line-line solver
            : Solver
    var lll_para_solver // < solver
            : Solver
    var qll_solver // < solver
            : Solver
    var sep_solver // < separator solver
            : Solver
    var alt_sep_solver // < alternative separator solver
            : Solver

    // DATA
    var g // < reference to the VD graph.
            : HalfEdgeDiagram
    var t_min = 0.0 // < minimum offset-distance
    var t_max = 0.0 // < maximum offset-distance
    lateinit var edge: Edge  // < the edge on which we position a new vertex
    var errstat: MutableList<Double> = mutableListOf() // < error-statistics
    var silent // < silent mode (outputs no warnings to stdout)
            : Boolean

    // create positioner, set graph.
    init {
        g = gi
        ppp_solver = PPPSolver()
        lll_solver = LLLSolver()
        qll_solver = QLLSolver()
        sep_solver = SEPSolver()
        alt_sep_solver = ALTSEPSolver()
        lll_para_solver = LLLPARASolver()
        silent = false
        errstat.clear()
    }

    // \brief position a new vertex on given HEEdge \a e when inserting the new
    // Site \a s3
    ///
    // calculate the position of a new voronoi-vertex lying on the given edge.
    // The new vertex is equidistant to the two sites that defined the edge
    // and to the new site.
    // the edge e holds information about which face it belongs to.
    // each face holds information about which site created it
    // so the three sites defining the position of the vertex are:
    // - site to the left of HEEdge e
    // - site to the right of HEEdge e
    // - given new Site s
    fun position(e: Edge, s3: Site): Solution {
        edge = e
        val face = e.face
        val twin = e.twin!!
        val twin_face = twin.face
        val src = e.source
        val trg = e.target
        val t_src = src.dist()
        val t_trg = trg.dist()
        t_min = min(t_src, t_trg) // the solution we seek must have t_min<t<t_max
        t_max = max(t_src, t_trg)
        val s1 = face.site
        val s2 = twin_face.site
        val sl: Solution = position(s1, e.k, s2, twin.k, s3)
        assert(solution_on_edge(sl)) { " solution_on_edge(sl) " }
        assert(check_dist(edge, sl, s3)) { " check_dist(edge, sl, s3) " }
        return sl
    }

    // position new vertex
    // find vertex that is equidistant from s1, s2, s3
    // should lie on the k1 side of s1, k2 side of s2
    // we try both k3=-1 and k3=+1 for s3
    fun position(s1: Site, k1: Double, s2: Site, k2: Double, s3: Site): Solution {
        assert(k1 == 1.0 || k1 == -1.0) { " (k1==1) || (k1 == -1) " }
        assert(k2 == 1.0 || k2 == -1.0) { " (k2==1) || (k2 == -1) " }
        if (s3.isLine) {
            // special handling for the case when site and edge endpoints share a common
            // point -
            // simply select one of the edge endpoints as solution
            val e_src = edge.source.position
            if ((s1.isPoint && s1.position().equals(e_src) || s1.isLine && (s1.start().equals(e_src) || s1.end()
                    .equals(e_src)))
                && (s2.isPoint && s2.position().equals(e_src) || s2.isLine && (s2.start().equals(e_src) || s2.end()
                    .equals(e_src)))
                && (s3.isPoint && s3.position().equals(e_src) || s3.isLine && (s3.start().equals(e_src) || s3.end()
                    .equals(e_src)))
            ) {
                val src_se = s3.start()
                val trg_se = s3.end()
                val k: Double
                k = if (edge.target.position.is_right(src_se, trg_se)) {
                    (if (s3.k().toInt() == 1) -1 else 1).toDouble()
                } else {
                    (if (s3.k().toInt() == 1) 1 else -1).toDouble()
                }
                return Solution(edge.source.position, edge.source.dist(), k)
            }
            val e_trg = edge.target.position
            if ((s1.isPoint && s1.position().equals(e_trg) || s1.isLine && (s1.start().equals(e_trg) || s1.end()
                    .equals(e_trg)))
                && (s2.isPoint && s2.position().equals(e_trg) || s2.isLine && (s2.start().equals(e_trg) || s2.end()
                    .equals(e_trg)))
                && (s3.isPoint && s3.position().equals(e_trg) || s3.isLine && (s3.start().equals(e_trg) || s3.end()
                    .equals(e_trg)))
            ) {
                val src_se = s3.start()
                val trg_se = s3.end()
                val k: Double
                k = if (edge.source.position.is_right(src_se, trg_se)) {
                    (if (s3.k().toInt() == 1) -1 else 1).toDouble()
                } else {
                    (if (s3.k().toInt() == 1) 1 else -1).toDouble()
                }
                return Solution(edge.target.position, edge.target.dist(), k)
            }
        }
        val solutions = mutableListOf<Solution>()
        if (s3.isLine && ((s1.isPoint && s2.isLine
                    && (s3.start().equals(s1.position()) || s3.end().equals(s1.position())))
                    || (s2.isPoint && s1.isLine
                    && (s3.start().equals(s2.position()) || s3.end().equals(s2.position()))))
        ) {
            val ptsite: Site = if (s1.isPoint) s1 else s2
            var ed = edge
            if (ed.face.site !== ptsite) {
                ed = edge.twin!!
            }
            assert(ed.source.status === VertexStatus.IN || ed.target.status === VertexStatus.IN) { "edge to be split has no IN vertex" }
            var k: Double
            k = if (ed.source.status === VertexStatus.IN) {
                -1.0
            } else {
                +1.0
            }
            if (s3.start().equals(ptsite.position())) {
                k = -k
            }
            solver_dispatch(s1, k1, s2, k2, s3, k, solutions)
        } else {
            solver_dispatch(s1, k1, s2, k2, s3, +1.0, solutions) // a single k3=+1 call for s3->isPoint
            if (!s3.isPoint) {
                solver_dispatch(s1, k1, s2, k2, s3, -1.0, solutions) // for lineSite or ArcSite we try k3=-1 also
            }
        }
        if (solutions.size == 1 && t_min <= solutions[0].t && t_max >= solutions[0].t
            && s3.in_region(solutions[0].p)
        ) {
            return solutions[0]
        }

        // choose only in_region() solutions
        val acceptable_solutions = mutableListOf<Solution>()
        for (s in solutions) {
            if (s3.in_region(s.p) && s.t >= t_min && s.t <= t_max) {
                acceptable_solutions.add(s)
            }
        }
        if (acceptable_solutions.size == 1) { // if only one solution is found, return that.
            return acceptable_solutions[0]
        } else if (acceptable_solutions.size > 1) {
            // two or more points remain so we must further filter here!
            // filter further using edge_error
            var min_error = 100.0
            var min_solution = Solution(Point(0.0, 0.0), 0.0, 0.0)
            for (s in acceptable_solutions) {
                val err = edge_error(s)
                if (err < min_error) {
                    min_solution = s
                    min_error = err
                }
            }
            return min_solution
        }
        return if (solutions.isEmpty()) {
            desperate_solution(s3)
        } else {
            // choose solution that is best by dist_error
            var leastBad: Solution = solutions[0]
            var leastErr = Double.MAX_VALUE
            for (s in solutions) {
                // punish wrong solutions
                val derr = dist_error(edge, s, s3)
                // punish solutions outside t range
                var terr = max(0.0, max(s.t - t_max, t_min - s.t))
                if (edge.type === EdgeType.PARA_LINELINE) {
                    val s_p = s.p.sub(edge.source.position)
                    val s_e = edge.target.position.sub(edge.source.position)
                    val dist = s_p.dot(s_e) / s_e.dot(s_e)
                    terr = max(0.0, max(dist - 1, -dist))
                }
                val err: Double = derr + terr
                if (err < leastErr) {
                    leastBad = s
                    leastErr = err
                }
            }
            if (edge.type === EdgeType.PARA_LINELINE) {
                return leastBad
            }

            // determine clamp direction
            val t = max(t_min, min(t_max, leastBad.t))
            val p_sln = edge.point(t)

            // find out on which side the solution lies
            var desp_k3 = 0.0
            if (s3.isPoint) {
                desp_k3 = 1.0
            } else if (s3.isLine) {
                val src_se = s3.start()
                val trg_se = s3.end()
                desp_k3 = if (p_sln.is_right(src_se, trg_se)) {
                    (if (s3.k() == 1.0) -1 else 1).toDouble()
                } else {
                    (if (s3.k() == 1.0) 1 else -1).toDouble()
                }
            }
            Solution(p_sln, t, desp_k3)
        }
    }

    // search numerically for a desperate solution along the solution-edge
    fun desperate_solution(s3: Site): Solution {
        val err_functor = VertexError(g, edge, s3)
        val src = edge.source
        val trg = edge.target
        val src_p = src.position
        val trg_p = trg.position

        val t_sln = minimize(t_min, t_max, 1e-14, 1000, err_functor::value).x
        val p_sln = err_functor.edge_point(t_sln) // g[edge].point(t_sln);
        var desp_k3 = 0.0
        if (s3.isPoint) {
            desp_k3 = 1.0
        } else if (s3.isLine) {
            // find out on which side the desperate solution lies
            val src_se = s3.start()
            val trg_se = s3.end()
            desp_k3 = if (p_sln.is_right(src_se, trg_se)) {
                (if (s3.k() == 1.0) -1 else 1).toDouble()
            } else {
                (if (s3.k() == 1.0) 1 else -1).toDouble()
            }
        }
        return Solution(p_sln, t_sln, desp_k3)
    }

    // dispatch to the correct solver based on the sites
    fun solver_dispatch(s1: Site, k1: Double, s2: Site, k2: Double, s3: Site, k3: Double, solns: MutableList<Solution>): Int {
        var s1: Site = s1
        var k1 = k1
        var s2: Site = s2
        var k2 = k2
        if (edge.type === EdgeType.SEPARATOR) {
            // this is a SEPARATOR edge with two LineSites adjacent.
            // find the PointSite that defines the SEPARATOR, so that one LineSite and one
            // PointSite
            // can be submitted to the Solver.
            if (s1.isLine && s2.isLine) {
                // the parallell lineseg case v0 --s1 --> pt -- s2 --> v1
                // find t
                if (edge.has_null_face) {
                    s2 = edge.null_face!!.site
                    assert(s2.isPoint) { " s2.isPoint " }
                    k2 = +1.0
                } else if (edge.twin!!.has_null_face) {
                    s2 = edge.twin!!.null_face!!.site
                    assert(s2.isPoint) { " s2.isPoint " }
                    k2 = +1.0
                }
            } else if (s1.isPoint && s2.isLine) {
                // a normal SEPARATOR edge, defined by a PointSite and a LineSite
                // swap sites, so SEPSolver can assume s1=line s2=point
                val tmp: Site = s1
                val k_tmp = k1
                s1 = s2
                s2 = tmp
                k1 = k2
                k2 = k_tmp
                assert(s1.isLine) { " s1.isLine " }
                assert(s2.isPoint) { " s2.isPoint " }
            }
            assert(s1.isLine && s2.isPoint) { " s1.isLine && s2.isPoint " }
            return sep_solver.solve(s1, k1, s2, k2, s3, k3, solns)
        } else if (edge.type === EdgeType.PARA_LINELINE && s3.isLine) { // an edge betwee parallel LineSites
            // std::cout << " para lineline! \n";
            return lll_para_solver.solve(s1, k1, s2, k2, s3, k3, solns)
        } else if (s1.isLine && s2.isLine && s3.isLine) {
            return lll_solver.solve(s1, k1, s2, k2, s3, k3, solns) // all lines.
        } else if (s1.isPoint && s2.isPoint && s3.isPoint) {
            return ppp_solver.solve(s1, 1.0, s2, 1.0, s3, 1.0, solns) // all points, no need to specify k1,k2,k3, they are
            // all +1
        } else if (s3.isLine && s1.isPoint || s1.isLine && s3.isPoint || s3.isLine && s2.isPoint || s2.isLine && s3.isPoint // bad coverage for this line?
        ) {
            // if s1/s2 form a SEPARATOR-edge, this is dispatched automatically to
            // sep-solver
            // here we detect for a separator case between
            // s1/s3
            // s2/s3
            if (s3.isLine && s1.isPoint) {
                if (detect_sep_case(s3, s1)) {
                    alt_sep_solver.type = 0
                    return alt_sep_solver.solve(s1, k1, s2, k2, s3, k3, solns)
                }
            }
            if (s3.isLine && s2.isPoint) {
                if (detect_sep_case(s3, s2)) {
                    alt_sep_solver.type = 1
                    return alt_sep_solver.solve(s1, k1, s2, k2, s3, k3, solns)
                }
            }
        }

        // if we didn't dispatch to a solver above, we try the general solver
        return qll_solver.solve(s1, k1, s2, k2, s3, k3, solns) // general case solver
    }

    // detect separator-case, so we can dispatch to the correct Solver
    fun detect_sep_case(lsite: Site, psite: Site): Boolean {
        val le = lsite.edge()!!
        val src = le.source
        val trg = le.target
        // now from segment end-points get the null-vertex
        lateinit var src_out: Edge
        for (e in src.out_edges) {
            if (e.type === EdgeType.NULLEDGE) {
                src_out = e
            }
        }
        lateinit var trg_out: Edge
        for (e in trg.out_edges) {
            if (e.type === EdgeType.NULLEDGE) {
                trg_out = e
            }
        }
        var src_null_face = src_out.face
        if (src_null_face.is_null_face == false) {
            // take twin face instead
            val src_out_twin = src_out.twin!!
            src_null_face = src_out_twin.face
        }
        var trg_null_face = trg_out.face
        if (trg_null_face.is_null_face == false) {
            val trg_out_twin = trg_out.twin!!
            trg_null_face = trg_out_twin.face
        }
        assert(src_null_face.is_null_face && trg_null_face.is_null_face) { " src_null_face.is_null_face && trg_null_face.is_null_face " }

        // do we want src_out face??
        // OR src_out_twin face??
        // we want the null-face !
        val src_site = src_null_face.site
        val trg_site = trg_null_face.site
        if (src_site == null || trg_site == null) {
            throw RuntimeException()
        }
        if (!src_site.isPoint || !trg_site.isPoint) {
            throw RuntimeException()
        }
        val src_vertex = src_site.vertex()
        val trg_vertex = trg_site.vertex()
        if (src_vertex == psite.vertex()) {
            return true
        }
        return if (trg_vertex == psite.vertex()) {
            true
        } else false
    }

    // error from solution to corresponding point on the edge
    fun edge_error(sl: Solution): Double {
        val p: Point
        p = if (edge.type === EdgeType.PARA_LINELINE) {
            projection_point(sl)
        } else {
            edge.point(sl.t)
        }
        return p.sub(sl.p).norm()
    }

    // when the edge is not parametrized by t-value as normal edges
    // so we need a projection of sl onto the edge instead
    fun projection_point(sl: Solution): Point {
        assert(edge.type === EdgeType.PARA_LINELINE) { " edge.type == EdgeType.PARA_LINELINE " }
        // edge given by
        // p = p0 + t * (p1-p0) with t in [0,1]
        val p0 = Point(edge.source.position)
        val p1 = Point(edge.target.position)
        val v = p1.sub(p0)
        var t = sl.p.sub(p0).dot(v) / v.dot(v)
        // clamp to [0,1]
        if (t > 1) {
            t = 1.0
        } else if (t < 0) {
            t = 0.0
        }
        return p0.add(v.mult(t))
    }

    // check that the new solution lies on the edge
    fun solution_on_edge(s: Solution): Boolean {
        val err = edge_error(s)
        val limit = 9E-4
        return err < limit
    }

    // new vertices should lie within the far_radius
    fun check_far_circle(s: Solution): Boolean {
        return if (s.p.norm() >= 18 * 1) {
            false
        } else true
    }

    // distance sanity check
    // all vertices should be of degree three, i.e. three adjacent faces/sites
    // distance to the three adjacent sites should be equal
    fun check_dist(e: Edge, sl: Solution, s3: Site): Boolean {
        val face = e.face
        val tw_edge = e.twin!!
        val twin_face = tw_edge.face
        val s1 = face.site
        val s2 = twin_face.site
        val d1 = sl.p.sub(s1.apex_point(sl.p)).norm()
        val d2 = sl.p.sub(s2.apex_point(sl.p)).norm()
        val d3 = sl.p.sub(s3.apex_point(sl.p)).norm()
        val maxd = max(
            max(abs(sl.t - d1), abs(sl.t - d2)),
            abs(sl.t - d3)
        )
        errstat.add(maxd)
        return !(!equal(d1, d2) || !equal(d1, d3) || !equal(d2, d3) || !equal(sl.t, d1) || !equal(sl.t, d2)
                || !equal(sl.t, d3))
    }

    // distance-error
    // new vertices should be equidistant to the three adjacent sites that define
    // the vertex
    // we here calculate the distances d1, d2, d3 from the Solution to the three
    // sites s1, s2, s3
    // and return the max deviation from the solution t-value.
    // this works as a sanity check for the solver.
    // a high error value here is also an indication of numerical instability in the
    // solver
    fun dist_error(e: Edge, sl: Solution, s3: Site): Double {
        val face = e.face
        val tw_edge = e.twin!!
        val twin_face = tw_edge.face
        val s1 = face.site
        val s2 = twin_face.site
        val d1 = sl.p.sub(s1.apex_point(sl.p)).norm()
        val d2 = sl.p.sub(s2.apex_point(sl.p)).norm()
        val d3 = sl.p.sub(s3.apex_point(sl.p)).norm()
        return max(
            max(abs(sl.t - d1), abs(sl.t - d2)),
            abs(sl.t - d3)
        )
    }

    // are \a d1 and \a d2 roughly equal?
    fun equal(d1: Double, d2: Double): Boolean {
        val tol = 1e-3
        if (abs(d1 - d2) < 1e-15) {
            return true
        }
        return if (abs(d1 - d2) > tol * max(d1, d2)) {
            false
        } else true
    }
}
