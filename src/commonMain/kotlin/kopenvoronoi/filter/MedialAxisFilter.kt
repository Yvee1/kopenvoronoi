package kopenvoronoi.filter

import assert
import kopenvoronoi.geometry.Edge
import kopenvoronoi.geometry.EdgeType
import kopenvoronoi.vertex.Vertex
import kopenvoronoi.vertex.VertexType

/**
 * Filter for retaining the medial-axis of a voronoi diagram (approximate).
 *
 * Marks the valid-property true for edges belonging to the medial axis and
 * false for other edges.
 *
 * @author MCarleton
 */
class MedialAxisFilter : Filter {
    /**
     * A dot-product threshold in [0,1] for filtering out edges between nearly
     * parallel LineSite segments
     */
    var _dot_product_threshold: Double

    constructor() {
        _dot_product_threshold = 0.8
    }

    /**
     *
     * @param thr A dot-product threshold in [0,1] for filtering out edges between
     * nearly parallel LineSite segments
     */
    constructor(thr: Double) {
        _dot_product_threshold = thr
    }

    // predicate that decides if an edge is to be included or not.
    override fun apply(e: Edge): Boolean {
        if (e.type === EdgeType.LINESITE || e.type === EdgeType.NULLEDGE) {
            return true // we keep linesites and nulledges
        }
        if (e.type === EdgeType.SEPARATOR) {
            return false // separators are allways removed
        }
        if (both_endpoints_positive(e)) {
            return true
        }

        // this leaves us with edges where one end connects to the polygon (dist==0)
        // and the other end does not.
        // figure out the angle between the adjacent line-segments and decide based on
        // the angle.
        return if (segments_parallel(e)) {
            false
        } else true
        // otherwise we keep the edge
    }

    // return true if this is an internal edge, i.e. both endpoints have a nonzero
    // clearance-disk radius
    private fun both_endpoints_positive(e: Edge): Boolean {
        val src = e.source
        val trg = e.target
        return src.dist() > 0 && trg.dist() > 0
    }

    // return true if the segments that connect to the given Edge are nearly
    // parallel
    private fun segments_parallel(e: Edge): Boolean {
        val endp1: Vertex = find_endpoint(e)
        val endp2: Vertex = find_endpoint(e.twin!!)
        // find the segments
        val e1: Edge = find_segment(endp1)
        var e2: Edge = find_segment(endp2)
        e2 = e2.twin!! // this makes the edges oriented in the same direction
        val dotprod = edge_dotprod(e1, e2)
        return dotprod > _dot_product_threshold
    }

    // \brief calculate the dot-product between unit vectors aligned along edges
    // e1->e2
    ///
    // since e1 and e2 are both line-sites the direction is easy to find
    // FIXME: more code needed here for tangent calculation if we have arc-sites
    private fun edge_dotprod(e1: Edge, e2: Edge): Double {
        val src1 = e1.source
        val trg1 = e1.target
        val src2 = e2.source
        val trg2 = e2.target
        val sp1 = src1.position
        val tp1 = trg1.position
        val sp2 = src2.position
        val tp2 = trg2.position
        val dir1 = tp1.sub(sp1)
        val dir2 = tp2.sub(sp2)
        dir1.normalize()
        dir2.normalize()
        return dir1.dot(dir2)
    }

    // find the LineSite edge that connects to \a v
    fun find_segment(v: Vertex): Edge {
        for (e in v.out_edges) {
            if (e.type === EdgeType.LINESITE) {
                return e
            }
        }
        throw RuntimeException("Failed to find line segment from vertex")
    }

    // find an ::ENDPOINT vertex that connects to Edge e through a ::NULLEDGE at
    // either the source or target of e.
    fun find_endpoint(e: Edge): Vertex {
        val next = e.next
        val prev = g!!.previous_edge(e)
        val endp: Vertex
        if (next.type === EdgeType.NULLEDGE) {
            endp = next.target
            assert(endp.type === VertexType.ENDPOINT) { "endp.type == VertexType.ENDPOINT " }
        } else if (prev.type === EdgeType.NULLEDGE) {
            endp = prev.source
            assert(endp.type === VertexType.ENDPOINT) { "endp.type == VertexType.ENDPOINT " }
        } else {
            throw RuntimeException("Failed to find endpoint")
        }
        return endp
    }
}
