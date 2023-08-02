package kopenvoronoi.offset

import assert
import kopenvoronoi.HalfEdgeDiagram
import kopenvoronoi.geometry.Edge
import kopenvoronoi.geometry.Face
import kopenvoronoi.geometry.Point
import kotlin.math.max
import kotlin.math.min

//\brief From a voronoi-diagram, generate offsets.
///
//an offset is allways a closed loop.
//the loop consists of offset-elements from each face that the loop visits.
//each face is associated with a Site, and the offset element from
//- a point-site is a circular arc
//- a line-site is a line
//- an arc is a circular arc
///
//This class produces offsets at the given offset-distance on the entire
//voronoi-diagram. To produce offsets only inside or outside a given geometry,
//use a filter first. The filter sets the valid-property of edges, so that offsets
//are not produced on faces with one or more invalid edge.
class Offset(g: HalfEdgeDiagram) {
    var g // < vd-graph
            : HalfEdgeDiagram
    var remaining_faces: MutableSet<Face> = mutableSetOf()
    var offset_list: MutableList<OffsetLoop>? = null // < list of output offsets

    // \param gi vd-graph
    init {
        this.g = g
    }

    // create offsets at offset distance \a t
    fun offset(t: Double): List<OffsetLoop>? {
        offset_list = mutableListOf()
        set_flags(t)
        var start: Face?
        var c = 0
        while (find_start_face().also { start = it } != null) { // while there are faces that still require offsets
            offset_loop_walk(start!!, t) // start on the face, and do an offset loop
            if (c > 30000) {
                throw AssertionError("c > 30000, hang in offset walk")
            }
            c++
        }
        return offset_list
    }

    // find a suitable start face
    private fun find_start_face(): Face? {
        return if (!remaining_faces.isEmpty()) {
            remaining_faces.iterator().next()
        } else {
            null
        }
    }

    // perform an offset walk at given distance \a t,
    // starting at the given face
    private fun offset_loop_walk(start: Face, t: Double) {
        var out_in_mode = false
        val start_edge: Edge = find_next_offset_edge(start.edge, t, out_in_mode) // the first edge on the start-face
        var current_edge: Edge = start_edge
        val loop = OffsetLoop() // store the output in this loop
        loop.offset_distance = t
        loop.add(OffsetVertex(current_edge.point(t), current_edge))
        do {
            out_in_mode = edge_mode(current_edge, t)
            // find the next edge
            val next_edge: Edge = find_next_offset_edge(current_edge.next, t, out_in_mode)
            val current_face = current_edge.face
            loop.add(offset_element_from_face(current_face, current_edge, next_edge, t))
            remaining_faces.remove(current_face) // although we may revisit current_face (if it is non-convex), it
            // seems safe to mark it "done" here.
            current_edge = next_edge.twin!!
        } while (current_edge !== start_edge)
        offset_list!!.add(loop) // append the created loop to the output
    }

    // return an offset-element corresponding to the current face
    private fun offset_element_from_face(
        current_face: Face,
        current_edge: Edge,
        next_edge: Edge,
        t: Double
    ): OffsetVertex {
        val s = current_face.site
        val o = s.offset(current_edge.point(t), next_edge.point(t)) // ask the Site for offset-geometry here.
        var cw = true
        if (!s.isLine) { // point and arc-sites produce arc-offsets, for which cw must be set.
            cw = find_cw(o.start(), o.center(), o.end()) // figure out cw or ccw arcs?
        }
        // add offset to output
        return OffsetVertex(next_edge.point(t), o.radius(), o.center(), cw, current_face, next_edge)
    }

    // \brief figure out mode (?)
    private fun edge_mode(e: Edge, t: Double): Boolean {
        val src = e.source
        val trg = e.target
        val src_r = src.dist()
        val trg_r = trg.dist()
        return if (src_r < t && t < trg_r) {
            true
        } else if (trg_r < t && t < src_r) {
            false
        } else {
            assert(false) { "failed to determine edge mode" }
            false
        }
    }

    // figure out cw or ccw for an arc
    private fun find_cw(start: Point, center: Point, end: Point): Boolean {
        return center.is_right(start, end) // NOTE: this only works for arcs smaller than a half-circle !
    }

    // \brief starting at e, find the next edge on the face that brackets t
    ///
    // we can be in one of two modes.
    // if mode==false then we are looking for an edge where src_t < t < trg_t
    // if mode==true we are looning for an edge where trg_t < t < src_t
    private fun find_next_offset_edge(e: Edge, t: Double, mode: Boolean): Edge {
        val start: Edge = e
        var current: Edge = start
        var ofs_edge: Edge = e
        do {
            val src = current.source
            val trg = current.target
            val src_r = src.dist()
            val trg_r = trg.dist()
            if (!mode && src_r < t && t < trg_r) {
                ofs_edge = current
                break
            } else if (mode && trg_r < t && t < src_r) {
                ofs_edge = current
                break
            }
            current = current.next
        } while (current !== start)
        return ofs_edge
    }

    // go through all faces and set flag=0 if the face requires an offset.
    private fun set_flags(t: Double) {
        for (f in g.faces) {
            val start = f.edge
            var current = start
            do {
                val src = current.source
                val trg = current.target
                val src_r = src.dist()
                val trg_r = trg.dist()
                if (t_bracket(src_r, trg_r, t)) {
                    remaining_faces.add(f)
                }
                current = current.next
            } while (current != start)
        }

        // again go through faces again, and set flag=1 if any edge on the face is
        // invalid
        // this is required because an upstream filter will set valid=false on some
        // edges,
        // but not all, on a face where we do not want offsets.
        for (f in g.faces) {
            val start = f.edge
            var current = start
            do {
                if (!current.valid) {
                    remaining_faces.remove(f) // don't offset faces with invalid edges
                }
                current = current.next
            } while (current != start)
        }
    }

    // is t in (a,b) ?
    private fun t_bracket(a: Double, b: Double, t: Double): Boolean {
        val min_t: Double = min(a, b)
        val max_t: Double = max(a, b)
        return min_t < t && t < max_t
    }
}
