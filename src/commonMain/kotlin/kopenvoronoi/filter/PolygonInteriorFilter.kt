package kopenvoronoi.filter

import kopenvoronoi.geometry.Edge
import kopenvoronoi.geometry.EdgeType
import kopenvoronoi.geometry.Face

//\brief Filter for retaining voronoi-diagram inside a polygon
///
//this filter sets the valid-property of edges
//all interior edges are marked valid=true
//all exterior edges are marked valid=false
///
//a polygon/pocket boundary shoud be specified in CW order
//islands within the polygon should be specified in CCW order
class PolygonInteriorFilter // \brief create a polygon interior Filter with given \a side
// \param side set true (false) for polygons inserted in CW (CCW) order and
// islands inserted in CCW (CW) order.
    (private val side: Boolean) : Filter() {
    // determine if an edge is valid or not
    override fun apply(e: Edge): Boolean {
        if (e.type === EdgeType.LINESITE || e.type === EdgeType.NULLEDGE) {
            return true
        }

        // if polygon inserted ccw as (id1->id2), then the linesite should occur on
        // valid faces as id1->id2
        // for islands and the outside the edge is id2->id1
        val f = e.face
        val s = f.site
        if (s.isLine && linesite_ccw(f)) {
            return true
        } else if (s.isPoint) {
            // we need to search for an adjacent linesite.
            // (? can we have a situation where this fails?)
            val linetwin: Edge? = find_adjacent_linesite(f)
            if (linetwin != null) {
                val twin = linetwin.twin
                val twin_face = twin!!.face
                if (linesite_ccw(twin_face)) {
                    return true
                }
            } else {
                return false
            }
        }
        return false
    }

    // on the face f, find the adjacent linesite
    private fun find_adjacent_linesite(f: Face): Edge? {
        var current = f.edge
        val start = current
        do {
            val twin = current.twin
            if (twin != null) {
                val twf = twin.face
                if (twf.site.isLine) {
                    return current
                }
            }
            current = current.next
        } while (current != start)
        return null
    }

    // return true if linesite was inserted in the direction indicated by _side
    private fun linesite_ccw(f: Face): Boolean {
        var current = f.edge
        val start = current
        do {
            if ((if (current.inserted_direction) side && current.type === EdgeType.LINESITE else !side && current.type === EdgeType.LINESITE)) {
                return true
            }
            current = current.next
        } while (current != start)
        return false
    }
}
