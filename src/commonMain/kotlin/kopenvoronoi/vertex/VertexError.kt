package kopenvoronoi.vertex

import kopenvoronoi.numerical.UnivariateFunction
import kopenvoronoi.HalfEdgeDiagram
import kopenvoronoi.geometry.Edge
import kopenvoronoi.geometry.EdgeType
import kopenvoronoi.geometry.Point
import kopenvoronoi.site.Site
import kotlin.math.abs

//\brief error functor for edge-based desperate solver
///
//minimize error by searching for a point on the solution-edge
class VertexError(gi: HalfEdgeDiagram, sln_edge: Edge, si3: Site) : UnivariateFunction {
    var g // < vd-graph
            : HalfEdgeDiagram
    var edge // < existing edge on which we have positioned a new vertex
            : Edge
    var s3 // < newly inserted Site
            : Site

    // \param gi vd-graph
    // \param sln_edge solution edge
    // \param si3 newly inserted Site
    init {
        g = gi
        edge = sln_edge
        s3 = si3
    }

    // return the vertex-error t-d3 where
    // t3 is the distance from edge-point(t) to s3, and
    // t is the offset-distance of the solution
    override fun value(t: Double): Double {
        val p: Point = edge_point(t)
        val s3_dist = p.sub(s3.apex_point(p)).norm()
        return abs(t - s3_dist)
    }

    // return a point on the edge at given offset-distance
    // \param t offset-distance ( >= 0 )
    fun edge_point(t: Double): Point {
        val p: Point
        if (edge.type === EdgeType.LINELINE) { // this is a workaround because the LINELINE edge-parameters are wrong? at
            // least in some cases?
            val src = edge.source
            val trg = edge.target
            val src_p = src.position
            val trg_p = trg.position
            val src_t = src.dist()
            val trg_t = trg.dist()
            // edge is src_p -> trg_p
            if (trg_t > src_t) {
                val frac: Double = (t - src_t) / (trg_t - src_t)
                p = src_p.add(trg_p.sub(src_p).mult(frac))
            } else {
                val frac: Double = (t - trg_t) / (src_t - trg_t)
                p = trg_p.add(src_p.sub(trg_p).mult(frac))
            }
        } else {
            p = edge.point(t)
        }
        return p
    }
}
