package kopenvoronoi

import assert
import kopenvoronoi.geometry.Edge
import kopenvoronoi.geometry.Face
import kopenvoronoi.vertex.Vertex

/**
 * bundled BGL properties, see:
 * http:www.boost.org/doc/libs/1_44_0/libs/graph/doc/bundles.html
 *
 * dcel notes from http:www.holmes3d.net/graphics/dcel/
 *
 *
 * vertex (boost::out_edges) -leaving pointer to HalfEdge that has this vertex
 * as origin if many HalfEdges have this vertex as origin, choose one
 * arbitrarily
 *
 *
 * HalfEdge - origin pointer to vertex (boost::source) - face to the left of
 * halfedge - twin pointer to HalfEdge (on the right of this edge) - next
 * pointer to HalfEdge this edge starts from h->twin->origin and ends at next
 * vertex in h->face traveling ccw around boundary (allows face traverse, follow
 * h->next until we arrive back at h)
 *
 *
 * Face - edge pointer to HalfEdge this edge has this Face object as face
 * half-edge can be any one on the boundary of face special "infinite face",
 * face on "outside" of boundary may or may not store edge pointer
 *
 *
 * \brief half-edge diagram, based on the boost graph-library /
 * half_edge_diagram is a half-edge diagram class. Templated on Vertex/Edge/Face
 * property classes which allow attaching information to vertices/edges/faces
 * that is required for a particular algorithm.
 *
 *
 * Inherits from boost::adjacency_list minor additions allow storing
 * face-properties.
 *
 *
 * the hedi namespace contains functions for manipulating HEDIGraphs
 *
 *
 * For a general description of the half-edge data structure see e.g.: -
 * http:www.holmes3d.net/graphics/dcel/ - http:openmesh.org/index.php?id=228
 */
class HalfEdgeDiagram {
    var vertices: ArrayList<Vertex> = ArrayList<Vertex>()
    var edges: ArrayList<Edge> = ArrayList<Edge>()
    var faces: ArrayList<Face> = ArrayList<Face>()
    fun add_vertex(): Vertex {
        val v = Vertex()
        vertices.add(v)
        return v
    }

    // add a vertex with given properties, return vertex descriptor
    fun add_vertex(v: Vertex): Vertex {
        vertices.add(v)
        return v
    }

    // return number of faces in graph
    fun num_faces(): Int {
        return faces.size
    }

    // return number of vertices in graph
    fun num_vertices(): Int {
        return vertices.size
    }

    // return number of edges in graph
    fun num_edges(): Int {
        return edges.size
    }

    // return number of edges on Face f
    fun num_edges(f: Face): Int {
        return face_edges(f).size
    }

    // add an edge between vertices v1-v2
    fun add_edge(v1: Vertex, v2: Vertex): Edge {
        val e = Edge(v1, v2)
        v1.out_edges.add(e)
        v2.in_edges.add(e)
        edges.add(e)
        return e
    }

    // return true if v1-v2 edge exists
    fun has_edge(v1: Vertex, v2: Vertex): Boolean {
        for (e in v1.out_edges) {
            if (e.target === v2) {
                return true
            }
        }
        return false
    }

    // return v1-v2 Edge
    fun edge(v1: Vertex, v2: Vertex): Edge {
        for (e in v1.out_edges) {
            if (e.target === v2) {
                return e
            }
        }
        throw RuntimeException("Edge not found in graph!")
    }

    // clear given vertex. this removes all edges connecting to the vertex.
    fun clear_vertex(v: Vertex) {
        for (e in v.out_edges) {
            e.target.in_edges.remove(e)
            edges.remove(e)
        }
        v.out_edges.clear()
        for (e in v.in_edges) {
            e.source.out_edges.remove(e)
            edges.remove(e)
        }
        v.in_edges.clear()
    }

    // remove given vertex. call clear_vertex() before this!
    fun remove_vertex(v: Vertex?) {
        vertices.remove(v)
    }

    // remove given edge
    fun remove_edge(e: Edge) {
        e.source.out_edges.remove(e)
        e.target.in_edges.remove(e)
        edges.remove(e)
    }

    // delete a vertex. clear and remove.
    fun delete_vertex(v: Vertex) {
        clear_vertex(v)
        remove_vertex(v)
    }

    /**
     * Insert vertex v into the middle of Edge e
     */
    fun add_vertex_in_edge(v: Vertex, e: Edge) {
        // the vertex v is inserted into the middle of edge e
        // edge e and its twin are replaced by four new edges: e1,e2 and their twins
        // te2,te1
        // before: face
        // e
        // previous-> source ------> target -> next
        // tw_next<- tw_trg <----- tw_src <- tw_previous
        // twin
        // twin_face
        //
        // after: face
        // e1 e2
        // previous-> source -> v -> target -> next
        // tw_next<- tw_trg <- v <- tw_src <- tw_previous
        // te2 te1
        // twin_face
        //
        val e_twin = e.twin!!
        val esource = e.source
        val etarget = e.target
        val face = e.face
        val twin_face = e_twin.face
        val previous: Edge = previous_edge(e)
        val twin_previous: Edge = previous_edge(e_twin)
        assert(previous.face === e.face) { " previous.face == e.face " }
        assert(twin_previous.face === e_twin.face) { " twin_previous.face == e_twin.face " }
        val e1: Edge = add_edge(esource, v)
        val te2: Edge = add_edge(v, esource)
        e1.twin = te2
        te2.twin = e1
        val e2: Edge = add_edge(v, etarget)
        val te1: Edge = add_edge(etarget, v)
        e2.twin = te1
        te1.twin = e2

        // next-pointers
        previous.next = e1
        e1.next = e2
        e2.next = e.next
        twin_previous.next = te1
        te1.next = te2
        te2.next = e_twin.next

        // this copies params, face, k, type
        e1.copyFrom(e)
        e2.copyFrom(e)
        te1.copyFrom(e_twin)
        te2.copyFrom(e_twin)

        // update the faces
        face.edge = e1
        twin_face.edge = te1

        // finally, remove the old edge
        remove_edge(e)
        remove_edge(e_twin)
    }



    /**
     * Adds two edges: one from v1 to v2, and one from v2 to v1
     *
     * @return
     */
    fun add_twin_edges(v1: Vertex, v2: Vertex): Pair<Edge, Edge> {
        val e1: Edge = add_edge(v1, v2)
        val e2: Edge = add_edge(v2, v1)
        e1.twin = e2
        e2.twin = e1
        return Pair<Edge, Edge>(e1, e2)
    }

    // make e1 the twin of e2 (and vice versa)
    fun twin_edges(e1: Edge, e2: Edge) {
        assert(e1.target === e2.source) { "e1.target == e2.source" }
        assert(e1.source === e2.target) { "e1.source == e2.target" }
        e1.twin = e2
        e2.twin = e1
    }

    // add a face, with given properties
    fun add_face(): Face {
        val f = Face()
        faces.add(f)
        return f
    }

    // return all vertices adjecent to given vertex
    fun adjacent_vertices(v: Vertex): List<Vertex> {
        val adj: MutableList<Vertex> = mutableListOf()
        for (e in v.out_edges) {
            adj.add(e.target)
        }
        return adj
    }

    // return all vertices of given face
    fun face_vertices(face: Face): List<Vertex> {
        val verts: MutableList<Vertex> = mutableListOf()
        val startedge = face.edge // the edge where we start
        val start_target = startedge.target
        verts.add(start_target)
        var current = startedge.next
        var count = 0
        do {
            val current_target = current.target
            verts.add(current_target)
            assert(current.face === current.next.face) { "current.face == current.next.face" }
            current = current.next
            if (count > 30000) {
                throw AssertionError("count < 30000")
            }
            count++
        } while (current != startedge)
        return verts
    }

    // return edges of face f as a vector
    // NOTE: it is faster to write a do-while loop in client code than to call this
    // function!
    fun face_edges(f: Face): List<Edge> {
        val start_edge = f.edge
        var current_edge = start_edge
        val out: MutableList<Edge> = mutableListOf()
        do {
            out.add(current_edge)
            current_edge = current_edge.next
        } while (current_edge != start_edge)
        return out
    }

    // return the previous edge. traverses all edges in face until previous found.
    fun previous_edge(e: Edge): Edge {
        var previous = e.next
        while (previous.next !== e) {
            previous = previous.next
        }
        return previous
    }

    // return adjacent faces to the given vertex
    fun adjacent_faces(q: Vertex): List<Face> {
        val face_set: MutableSet<Face> = mutableSetOf()
        for (e in q.out_edges) {
            face_set.add(e.face)
        }
        return face_set.toList()
    }

    // remove given v1-v2 edge
    fun remove_edge(v1: Vertex, v2: Vertex) {
        remove_edge(edge(v1, v2))
    }

    // remove given v1-v2 edge and its twin
    fun remove_twin_edges(v1: Vertex, v2: Vertex) {
        assert(has_edge(v1, v2)) { " has_edge(v1,v2) " }
        assert(has_edge(v2, v1)) { " has_edge(v2,v1) " }
        remove_edge(edge(v1, v2))
        remove_edge(edge(v2, v1))
    }

    // remove a degree-two Vertex from the middle of an Edge
    // preserve edge-properties (next, face, k)
    fun remove_deg2_vertex(v: Vertex) {
        // face1 e[1]
        // v1_prev -> v1 -> SPLIT -> v2 -> v2_next
        // v1_next <- v1 <- SPLIT <- v2 <- v2_prev
        // e[0] face2
        //
        // is replaced with a single edge:
        // face1
        // v1_prev -> v1 ----------> v2 -> v2_next
        // v1_next <- v1 <---------- v2 <- v2_prev
        // face2
        val v_edges = v.out_edges
        assert(v_edges.size == 2) { " v_edges.size() == 2" }
        assert(v_edges.get(0).source === v && v_edges.get(1).source === v) { " v_edges.get(0).source == v && v_edges.get(1).source == v " }
        val v1 = v_edges.get(0).target
        val v2 = v_edges.get(1).target
        val v1_next = v_edges.get(0).next
        val v1_prev: Edge = previous_edge(v_edges.get(0).twin!!)
        val v2_next = v_edges.get(1).next
        val v2_prev: Edge = previous_edge(v_edges.get(1).twin!!)
        val face1 = v_edges.get(1).face
        val face2 = v_edges.get(0).face
        val (new1, new2) = add_twin_edges(v1, v2)
        set_next(new1, v2_next)
        set_next(new2, v1_next)
        set_next(v2_prev, new2)
        set_next(v1_prev, new1)
        face1.edge = new1
        face2.edge = new2
        new1.copyFrom(v_edges.get(1))
        new2.copyFrom(v_edges.get(0))
        remove_twin_edges(v, v1)
        remove_twin_edges(v, v2)
        remove_vertex(v)
    }

    // set next-pointer of e1 to e2
    fun set_next(e1: Edge, e2: Edge) {
        assert(e1.target === e2.source) { " e1.target == e2.source " }
        e1.next = e2
    }

    // form a face from the edge-list:
    // e1->e2->...->e1
    // for all edges, set edge.face=f, and edge.k=k
    fun set_next_cycle(list: List<Edge>, f: Face, k: Double) {
        f.edge = list[0]
        for (q in list.indices) {
            val e: Edge = list[q]
            e.face = f
            e.k = k
            if (q == list.size - 1) {
                set_next(e, list[0])
            } else {
                set_next(e, list[q + 1])
            }
        }
    }

    // set next-pointers for the given list (but don't close to form a cycle)
    // also set face and k properties for edge
    fun set_next_chain(list: List<Edge>, f: Face, k: Double) {
        f.edge = list[0]
        for (q in list.indices) {
            val e: Edge = list[q]
            e.face = f
            e.k = k
            if (q != list.size - 1) {
                set_next(e, list[q + 1])
            }
        }
    }

    // set next-pointers for the list
    fun set_next_chain(list: List<Edge>) {
        for (q in list.indices) {
            val e: Edge = list[q]
            if (q != list.size - 1) {
                set_next(e, list[q + 1])
            }
        }
    }

    // on a face, search and return the left/right edge from endp
    fun find_next_prev(f: Face, endp: Vertex): Pair<Edge, Edge> {
        var current = f.edge
        val start_edge = current
        var next_edge: Edge? = null
        var prev_edge: Edge? = null
        do {
            val src = current.source
            val trg = current.target
            if (src == endp) {
                next_edge = current
            }
            if (trg == endp) {
                prev_edge = current
            }
            current = current.next
        } while (current != start_edge)
        assert(next_edge != null) { " next_edge != null " }
        assert(prev_edge != null) { " prev_edge != null " }
        return Pair(next_edge!!, prev_edge!!)
    }
}
