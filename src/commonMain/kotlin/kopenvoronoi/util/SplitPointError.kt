package kopenvoronoi.util

import kopenvoronoi.numerical.UnivariateFunction
import kopenvoronoi.HalfEdgeDiagram
import kopenvoronoi.geometry.Edge
import kopenvoronoi.geometry.Point

//\brief error-functor to locate ::SPLIT vertices
///
//for passing to numerical boost::toms748 root-finding algorithm
class SplitPointError(gi: HalfEdgeDiagram, split_edge: Edge, pt1: Point, pt2: Point) : UnivariateFunction {
    private val g // < reference to vd-graph
            : HalfEdgeDiagram
    private val edge // < the HEEdge on which we position the new SPLIT vertex
            : Edge
    private val p1 // < first point of the split-line
            : Point
    private val p2 // < second point of the split-line
            : Point

    // \param gi graph
    // \param split_edge the edge on which we want to position a SPLIT vertex
    // \param pt1 first point of split-line
    // \param pt2 second point of split-line
    init {
        g = gi
        edge = split_edge
        p1 = pt1
        p2 = pt2
    }

    // \return signed distance to the pt1-pt2 line from edge-point at given offset
    // \a t
    override fun value(t: Double): Double {
        val p = edge.point(t)
        // line: pt1 + u*(pt2-pt1) = p
        // (p-pt1) dot (pt2-pt1) = u* (pt2-pt1) dot (pt2-pt1)
        val u = p.sub(p1).dot(p2.sub(p1)) / p2.sub(p1).dot(p2.sub(p1))
        val proj = p1.add(p2.sub(p1).mult(u))
        val dist = proj.sub(p).norm()
        val sign: Double
        sign = if (p.is_right(p1, p2)) {
            +1.0
        } else {
            -1.0
        }
        return sign * dist
    }
}
