package kopenvoronoi.util

import assert
import kopenvoronoi.HalfEdgeDiagram
import kopenvoronoi.geometry.Edge
import kopenvoronoi.geometry.Face
import kopenvoronoi.vertex.Vertex
import kopenvoronoi.vertex.VertexStatus

/**
 * Provides sanity-checks for the VoronoiDiagram class
 */
class VoronoiDiagramChecker(gi: HalfEdgeDiagram) {
    private val g // < vd-graph
            : HalfEdgeDiagram

    init {
        g = gi
    }

    // overall sanity-check for the diagram, calls other sanity-check functions
    fun is_valid(): Boolean {
        return all_faces_ok() && vertex_degree_ok() && face_count_equals_generator_count()
    }

    // check that number of faces equals the number of generators
    // \todo not implemented!
    fun face_count_equals_generator_count(): Boolean {
        // Euler formula for planar graphs
        // v - e + f = 2
        // in a half-edge diagram all edges occur twice, so:
        // f = 2-v+e
        // int vertex_count = hedi::num_vertices(g);
        /*
		 * int vertex_count = 0; BOOST_FOREACH( HEVertex v, hedi::vertices( g ) ) { if (
		 * g[v].type == NORMAL ) vertex_count++; } int face_count = (vertex_count- 4)/2
		 * + 3; // degree three graph //int face_count = hed.num_faces(); if (face_count
		 * != gen_count) { std::cout << " face_count_equals_generator_count() ERROR:\n";
		 * std::cout << " num_vertices = " << vertex_count << "\n"; std::cout <<
		 * " gen_count = " << gen_count << "\n"; std::cout << " face_count = " <<
		 * face_count << "\n"; } return ( face_count == gen_count );
		 */
        return true
    }

    // check that the diagram is of degree three.
    // however ::SPLIT and ::APEX vertices are of degree 2.
    fun vertex_degree_ok(): Boolean {
        for (v in g.vertices) {
            if (v.degree() !== Vertex.expected_degree.get(v.type)) {
                return false
            }
        }
        return true
    }

    // check that all vertices in the input vector have status ::IN
    fun all_in(q: List<Vertex>): Boolean {
        for (v in q) {
            if (v.status !== VertexStatus.IN) {
                return false
            }
        }
        return true
    }

    // check that no undecided vertices remain in the face
    fun noUndecidedInFace(f: Face): Boolean { // is this true??
        val face_verts = g.face_vertices(f)
        for (v in face_verts) {
            if (v.status === VertexStatus.UNDECIDED) {
                return false
            }
        }
        return true
    }

    // check that for HEFace f the vertices TYPE are connected
    fun faceVerticesConnected(f: Face, Vtype: VertexStatus): Boolean {
        val face_verts = g.face_vertices(f)
        val type_verts: MutableList<Vertex> = mutableListOf()
        for (v in face_verts) {
            if (v.status === Vtype) {
                type_verts.add(v) // build a vector of all Vtype vertices
            }
        }
        assert(!type_verts.isEmpty()) { " !type_verts.isEmpty() " }
        if (type_verts.size == 1) {
            return true
        }

        // check that type_verts are connected
        var currentEdge = f.edge
        val endVertex = currentEdge.source // stop when target here
        val startEdges: MutableList<Edge> = mutableListOf()
        var done = false
        while (!done) {
            val src = currentEdge.source
            val trg = currentEdge.target
            if (src.status !== Vtype) { // seach ?? - Vtype
                if (trg.status === Vtype) { // we have found ?? - Vtype
                    startEdges.add(currentEdge)
                }
            }
            currentEdge = currentEdge.next
            if (trg == endVertex) {
                done = true
            }
        }
        assert(!startEdges.isEmpty()) { " !startEdges.isEmpty() " }
        return if (startEdges.size != 1) {
            false
        } else {
            true
        }
    }

    // check that all faces are ok. calls face_ok()
    fun all_faces_ok(): Boolean {
        for (f in g.faces) {
            if (!face_ok(f)) {
                return false
            }
        }
        return true
    }

    // check that the face is ok
    fun face_ok(f: Face): Boolean {
        var current_edge = f.edge
        val start_edge = current_edge
        val k = current_edge.k
        if (!(k == 1.0 || k == -1.0)) {
            // std::cout << " VoronoiDiagramChecker::face_ok() f=" << f << " ERROR:\n";
            // std::cout << " illegal k-value for edge:";
            // std::cout << g[ g.source(current_edge)].index << " - ";
            // std::cout << g[ g.target(current_edge)].index ;
            // std::cout << " k= " << k << "\n";
            return false
        }
        if (f.site != null) { // guard against null-faces that dont have Site
            if (f.site.isPoint) {
                if (k != 1.0) {
                    // std::cout << " VoronoiDiagramChecker::face_ok() f=" << f << " ERROR:\n";
                    // std::cout << " f = " << f << " site is " << g[f].site->str() << " but k=" <<
                    // k << "\n";
                    // std::cout << " null? " << g[f].null << "\n";
                    return false
                }
            }
        }
        var n = 0
        do {
            if (current_edge.k != k) { // all edges should have the same k-value
                return false
            }
            if (!current_face_equals_next_face(current_edge)) { // all edges should have the same face
                return false
            }
            if (!check_edge(current_edge)) {
                return false
            }
            current_edge = current_edge.next
            n++
            assert(n < 10000) { " n < 10000 " }
        } while (current_edge != start_edge)
        return true
    }

    // check that current edge and next-edge are on the same face
    fun current_face_equals_next_face(e: Edge): Boolean {
        return if (e.face !== e.next.face) {
            false
        } else true
    }

    // sanity-check for edge
    fun check_edge(e: Edge): Boolean {
        val src = e.source
        val trg = e.target
        val twine = e.twin
        if (twine == null) {
            return true
        } else if (e !== twine.twin) {
            return false
        }
        val tw_src = twine.source
        val tw_trg = twine.target
        return src == tw_trg && trg == tw_src
    }
}
