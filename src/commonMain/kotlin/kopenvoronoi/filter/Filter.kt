package kopenvoronoi.filter

import kopenvoronoi.HalfEdgeDiagram
import kopenvoronoi.geometry.Edge

//\brief base-class for voronoi-diagram filters
///
//concrete sub-classes of Filter provide a predicate
//for determining if the edge belongs to the filtered graph.
abstract class Filter {
    protected var g: HalfEdgeDiagram? = null // < vd-graph

    // set graph
    fun set_graph(g: HalfEdgeDiagram?) {
        this.g = g
    }

    // does this edge belong to the filtered graph?
    abstract fun apply(e: Edge): Boolean
}
