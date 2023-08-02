package kopenvoronoi

import assert
import dataStructures.trees.thirdGenKD.KdTree
import dataStructures.trees.thirdGenKD.SquareEuclideanDistanceFunction
import kopenvoronoi.filter.Filter
import kopenvoronoi.geometry.*
import kopenvoronoi.geometry.*
import kopenvoronoi.site.LineSite
import kopenvoronoi.site.PointSite
import kopenvoronoi.site.Site
import kopenvoronoi.util.KdPoint
import kopenvoronoi.util.Numeric
import kopenvoronoi.util.SplitPointError
import kopenvoronoi.util.VoronoiDiagramChecker
import kopenvoronoi.vertex.Vertex
import kopenvoronoi.vertex.VertexPositioner
import kopenvoronoi.vertex.VertexStatus
import kopenvoronoi.vertex.VertexType
import kotlin.jvm.JvmOverloads
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min
import kotlin.math.sqrt
import kroot.toms748

/**
 * Voronoi diagram. The diagram builds incrementally as sites are added.
 *
 *
 * See http://en.wikipedia.org/wiki/Voronoi_diagram
 *
 *
 * The dual of a Voronoi diagram is the delaunay diagram (triangulation).
 * Voronoi-faces are dual to delaunay-vertices; Voronoi-vertices are dual to
 * delaunay-faces; Voronoi-edges are dual to delaunay-edges.
 */
class VoronoiDiagram @JvmOverloads constructor(farRadius: Double = 5000.0) {
    // HELPER-CLASSES
    /**
     * sanity-checks on the diagram are done by this helper class
     */
    protected var vd_checker: VoronoiDiagramChecker

    /**
     * kd-tree for nearest neighbour search during point Site insertion
     */
    protected var kd_tree: KdTree<KdPoint>

    /**
     * an algorithm for positioning vertices
     */
    protected var vpos: VertexPositioner
    // DATA
    /**
     * priority_queue for vertex for processing sorted by decreasing fabs() of
     * in_circle-predicate, so that the vertices whose IN/OUT status we are 'most
     * certain' about are processed first queue of vertices to be processed
     */
    protected var vertexQueue: PriorityQueue<Pair<Vertex, Double>> =
        PriorityQueue<Pair<Vertex, Double>>(abs_comparison())

    /**
     * the half-edge diagram of the vd
     */
    protected var g: HalfEdgeDiagram = HalfEdgeDiagram()

    /**
     * sites must fall within a circle with radius far_radius
     */
    protected var far_radius // <
            : Double

    /**
     * number of point sites
     */
    protected var num_psites: Int

    /**
     * number of line-segment sites
     */
    protected var num_lsites: Int

    /**
     * number of arc-sites
     */
    protected var num_asites: Int

    /**
     * temporary variable for INCIDENT faces, will be reset to NONINCIDENT after a
     * site has been inserted
     */
    protected var incident_faces: MutableList<Face> = mutableListOf()

    /**
     * temporary variable for in-vertices, out-vertices that need to be reset after
     * a site has been inserted
     */
    protected var modified_vertices: MutableSet<Vertex> = mutableSetOf()

    /**
     * IN-vertices, i.e. to-be-deleted
     */
    protected var v0: MutableList<Vertex> = mutableListOf()
    protected var debug // < turn debug output on/off
            : Boolean

    /**
     * Creates an empty VoronoiDiagram
     *
     * @param farRadius the radius of a circle within which all sites must be
     * located.
     */
    init {
        kd_tree = KdTree<KdPoint>(2)
        vd_checker = VoronoiDiagramChecker(g) // helper-class that checks topology/geometry
        vpos = VertexPositioner(g) // helper-class that positions vertices
        far_radius = farRadius
        initialize()
        num_psites = 3
        num_lsites = 0
        num_asites = 0
        reset_vertex_count()
        debug = false
    }

    /**
     * [.insert_point_site]
     */
    fun insert_point_site(x: Double, y: Double): Vertex? {
        return insert_point_site(Point(x, y))
    }

    /**
     * Insert a PointSite into the voronoi diagram.
     *
     *
     *
     * All PointSites must be inserted before any LineSites or ArcSites are
     * inserted.
     *
     *
     * It is an error to insert duplicate PointSites (i.e. points with the same x,y
     * coordinates). All PointSites must be inserted before any LineSites or
     * ArcSites are inserted.
     *
     *
     *
     * This is roughly "algorithm A" from the Sugihara-Iri 1994 paper, page 15/50 /
     * -# find the face that is closest to the new site, see FaceGrid -# among the
     * vertices on the closest face, find the seed vertex, see find_seed_vertex() -#
     * grow the tree of IN-vertices, see augment_vertex_set() -# add new
     * voronoi-vertices on all IN-OUT edges so they becone IN-NEW-OUT, see
     * add_vertices() -# add new face by splitting each INCIDENT face into two parts
     * by inserting a NEW-NEW edge. see add_edges() -# repair the next-pointers of
     * faces that have been modified. see repair_face() -# remove IN-IN edges and
     * IN-NEW edges, see remove_vertex_set() -# reset vertex/face status to be ready
     * for next incremental operation, see reset_status()
     *
     * @param p position of site
     * @return integer handle to the inserted point. use this integer when inserting
     * lines/arcs with insert_line_site
     */
    fun insert_point_site(p: Point): Vertex {
        if (p.norm() >= far_radius) {
            error(
                "openvoronoi error. All points must lie within unit-circle. You are trying to add p= $p with p.norm()= ${p.norm()}\n",
            )
        }
        assert(p.norm() < far_radius) { " p.norm() < far_radius " }
        val nearest = kd_tree.findNearestNeighbors(doubleArrayOf(p.x, p.y), 1, SquareEuclideanDistanceFunction())
            .max
        if (nearest.p.equals(p)) {
            error(
                "Cannot insert point site $p. A point with the same coordinates already exists."
            )
        }
        if (num_lsites > 0) {
            error("Skipping inserting $p because a line site has already been inserted. Insert all point sites first before inserting any line sites.")
        }
        num_psites++
        val new_vert = g.add_vertex(Vertex(p, VertexStatus.OUT, VertexType.POINTSITE))
        val new_site = PointSite(p)
        new_site.v = new_vert
        val v_seed: Vertex = find_seed_vertex(nearest.face, new_site)
        mark_vertex(v_seed, new_site)
        augment_vertex_set(new_site) // grow the tree to maximum size
        add_vertices(new_site) // insert new vertices on IN-OUT edges
        val newface: Face = add_face(new_site)
        new_vert.face = newface // Vertices that correspond to point-sites have their .face property set!
        for (f in incident_faces) {
            add_edges(newface, f) // add NEW-NEW edges on all INCIDENT faces
        }
        repair_face(newface)
        remove_vertex_set() // remove all IN vertices and adjacent edges
        reset_status() // reset all vertices to UNDECIDED
        assert(vd_checker.face_ok(newface)) { " vd_checker.face_ok( newface ) " }
        assert(vd_checker.is_valid()) { " vd_checker.is_valid() " }
        return new_vert
    }

    /**
     * Inserts a LineSite into the diagram
     *
     *
     * All PointSites must be inserted before any LineSites are inserted. All
     * LineSites should be inserted before any ArcSitess are inserted. It is an
     * error to insert a LineSite that intersects an existing LineSite in the
     * diagram!
     *
     *
     * The basic idea of the incremental diagram update is similar to that in
     * insert_point_site(). The major differences are: - the handling of null-faces
     * at the endpoints of the LineSite. - addition of SEPARATOR edges - addition of
     * SPLIT vertices during augment_vertex_set() - removal of SPLIT vertices at the
     * end
     *
     *
     * The steps of the algorithm are: -# based on \a idx1 and \a idx2, look up the
     * corresponding vertex descriptors. It is an error if these are not found. -#
     * find a seed-vertex -# grow the delete-tree of IN vertices. -# create or
     * modify the null-faces at the startpoint and endpoint of the LineSite -#
     * create and add LINESITE pseudo-edges -# add NEW vertices on all IN-OUT edges.
     * -# add up to four SEPARATOR edges, where applicable -# add non-separator
     * edges by calling add_edges() on all INCIDENT faces -# repair the
     * next-pointers of faces that have been modified. see repair_face() -# remove
     * IN-IN edges and IN-NEW edges, see remove_vertex_set() -# remove SPLIT
     * vertices -# reset vertex/face status to be ready for next incremental
     * operation, see reset_status()
     *
     * @param start startpoint of line-segment
     * @param end   endpoint of line-segment
     * @return
     */
    fun insert_line_site(start: Vertex, end: Vertex): Pair<Edge, Edge> {
        num_lsites++

        // find the vertices corresponding to idx1 and idx2
        start.status = VertexStatus.OUT
        end.status = VertexStatus.OUT
        start.zero_dist()
        end.zero_dist()

        // create a point which is left of src->trg
        // determine k (offset-dir) for this point
        // then we know which site/face is the k==+1 and which is k==-1
        val src_se = start.position
        val trg_se = end.position
        val left = src_se.add(trg_se).mult(0.5).add(trg_se.sub(src_se).xy_perp()) // this is used below and in
        // find_null_face()
        val pos_site = LineSite(end.position, start.position, +1.0)
        val neg_site = LineSite(start.position, end.position, -1.0)
        val seed_face = start.face!! // assumes this point-site has a face!

        // on the face of start-point, find the seed vertex
        val v_seed: Vertex = find_seed_vertex(seed_face, pos_site)
        mark_vertex(v_seed, pos_site)
        augment_vertex_set(pos_site) // it should not matter if we use pos_site or neg_site here
        // todo(?) sanity checks:
        // check that end_face is INCIDENT?
        // check that tree (i.e. v0) includes end_face_seed ?

        // process the null-faces here
        // returns new seg_start/end vertices, new or existing null-faces, and separator
        // endpoints (if separators should be added)
        val dir2 = start.position.sub(end.position)
        val dir1 = end.position.sub(start.position)
        val null_face_start = find_null_face(start, end, left, dir1, pos_site)
        val seg_start = null_face_start.get1!!
        val start_null_face = null_face_start.get2!!
        val pos_sep_start = null_face_start.get3
        val neg_sep_start = null_face_start.get4
        val start_to_null = null_face_start.get5
        val null_face_end = find_null_face(end, start, left, dir2, pos_site)
        val seg_end = null_face_end.get1!!
        val end_null_face = null_face_end.get2!!
        val pos_sep_end = null_face_end.get3
        val neg_sep_end = null_face_end.get4
        val end_to_null = null_face_end.get5

        // now safe to set the zero-face edge
        // in the collinear case, set the edge for the face that "disappears" to a null
        // edge
        if (start_to_null != null) {
            start_to_null.edge = start_null_face.edge
        }
        if (end_to_null != null) {
            end_to_null.edge = end_null_face.edge
        }

        // create LINESITE pseudo edges and faces
        val (pos_edge, neg_edge) = g.add_twin_edges(seg_end, seg_start)
        pos_edge.inserted_direction = false
        neg_edge.inserted_direction = true
        pos_edge.type = EdgeType.LINESITE
        neg_edge.type = EdgeType.LINESITE
        pos_edge.k = +1.0
        neg_edge.k = -1.0
        val pos_face: Face = add_face(pos_site) // this face to the left of start->end edge
        val neg_face: Face = add_face(neg_site) // this face is to the left of end->start edge
        pos_face.edge = pos_edge
        neg_face.edge = neg_edge
        pos_edge.face = pos_face
        neg_edge.face = neg_face

        // associate sites with LINESITE edges
        pos_site.e = pos_edge
        neg_site.e = neg_edge
        add_vertices(pos_site) // add NEW vertices on all IN-OUT edges.

        // add SEPARATORS
        // find SEPARATOR targets first
        val pos_start_target = find_separator_target(start.face!!, pos_sep_start)
        val neg_start_target = find_separator_target(start.face!!, neg_sep_start)

        // add positive separator edge at start
        add_separator(start.face!!, start_null_face, pos_start_target, pos_sep_start, pos_face.site, neg_face.site)

        // add negative separator edge at start
        add_separator(start.face!!, start_null_face, neg_start_target, neg_sep_start, pos_face.site, neg_face.site)
        start.face!!.status = FaceStatus.NONINCIDENT // face is now done.
        assert(vd_checker.face_ok(start.face!!)) { " vd_checker.face_ok( start.face ) " }
        val pos_end_target = find_separator_target(end.face!!, pos_sep_end)
        val neg_end_target = find_separator_target(end.face!!, neg_sep_end)
        // add positive separator edge at end
        add_separator(end.face!!, end_null_face, pos_end_target, pos_sep_end, pos_face.site, neg_face.site)

        // add negative separator edge at end
        add_separator(end.face!!, end_null_face, neg_end_target, neg_sep_end, pos_face.site, neg_face.site)
        end.face!!.status = FaceStatus.NONINCIDENT
        assert(vd_checker.face_ok(end.face!!)) { " vd_checker.face_ok( end.face ) " }

        // add non-separator edges by calling add_edges on all INCIDENT faces
        for (f in incident_faces) {
            if (f.status === FaceStatus.INCIDENT) { // end-point faces already dealt with in add_separator()
                add_edges(pos_face, f, neg_face, Pair<Vertex, Vertex>(seg_start, seg_end)) // each INCIDENT face is
                // split into two parts:
                // newface and f
            }
        }

        // new vertices and edges inserted. remove the delete-set, repair faces.
        remove_vertex_set()
        repair_face(
            pos_face, Pair<Vertex, Vertex>(seg_start, seg_end),
            Pair<Face?, Face?>(start_to_null, end_to_null), Pair<Face, Face>(start_null_face, end_null_face)
        )
        assert(vd_checker.face_ok(pos_face)) { " vd_checker.face_ok( pos_face ) " }
        repair_face(
            neg_face, Pair<Vertex, Vertex>(seg_start, seg_end),
            Pair<Face?, Face?>(start_to_null, end_to_null), Pair<Face, Face>(start_null_face, end_null_face)
        )
        assert(vd_checker.face_ok(neg_face)) { " vd_checker.face_ok( neg_face ) " }

        // we are done and can remove split-vertices
        for (f in incident_faces) {
            remove_split_vertex(f)
        }
        reset_status()
        assert(vd_checker.face_ok(start_null_face)) { " vd_checker.face_ok( start_null_face ) " }
        assert(vd_checker.face_ok(end_null_face)) { " vd_checker.face_ok( end_null_face ) " }
        assert(vd_checker.face_ok(pos_face)) { " vd_checker.face_ok( pos_face ) " }
        assert(vd_checker.face_ok(neg_face)) { " vd_checker.face_ok( neg_face ) " }
        assert(vd_checker.is_valid()) { " vd_checker.is_valid() " }
        return pos_edge to neg_edge
    }

    // return the far radius
    fun get_far_radius(): Double {
        return far_radius
    }

    val faces: ArrayList<Face>
        get() = g.faces
    val edges: ArrayList<Edge>
        get() = g.edges
    val vertices: ArrayList<Vertex>
        get() = g.vertices

    /**
     * Returns the face nearest to the given point.
     */
    fun nearestFace(x: Double, y: Double): Face {
        return nearestFaces(x, y, 1)[0]
    }

    /**
     * Returns the N nearest faces to the given point.
     */
    fun nearestFaces(x: Double, y: Double, n: Int): List<Face> {
        val heap = kd_tree.findNearestNeighbors(doubleArrayOf(x, y), n, SquareEuclideanDistanceFunction())
        val faces = ArrayList<Face>(heap.size())
        while (heap.size() > 0) {
            faces.add(heap.max.face)
            heap.removeMax()
        }
        return faces
    }

    // return number of point sites in diagram
    fun num_point_sites(): Int {
        // the three initial vertices don't count
        return num_psites - 3
    }

    // return number of line-segments sites in diagram
    fun num_line_sites(): Int {
        return num_lsites
    }

    // return number of arc-sites in diagram
    fun num_arc_sites(): Int {
        return num_asites
    }

    // return number of voronoi-vertices
    fun num_vertices(): Int {
        return g.num_vertices() - num_point_sites()
    }

    // return number of faces in graph
    fun num_faces(): Int {
        return g.num_faces()
    }

    // return number of SPLIT vertices
    fun num_split_vertices(): Int {
        var count = 0
        for (v in g.vertices) {
            if (v.type === VertexType.SPLIT) {
                count++
            }
        }
        return count
    }

    val diagram: HalfEdgeDiagram
        // return reference to graph \todo not elegant. only used by vd2svg ?
        get() = g

    // turn on debug output
    fun debug_on() {
        debug = true
    }

    // run topology/geometry check on diagram
    fun check(): Boolean {
        return if (vd_checker.is_valid()) {
            true
        } else {
            false
        }
    }

    // filter the graph using given Filter \a flt
    fun filter(flt: Filter) {
        flt.set_graph(g)
        for (e in g.edges) {
            if (!flt.apply(e)) {
                e.valid = false
            }
        }
    }

    // reset filtering by setting all edges valid
    fun filter_reset() {
        for (e in g.edges) {
            e.valid = true
        }
    }

    // comparison-predicate for VertexQueue
    ///
    // in augment_vertex_set() we grow the delete-tree by processing vertices
    // one-by-one from a priority_queue. This is the priority_queue sort predicate.
    // We handle vertices with a large fabs( in_circle() ) first, since we
    // believe their predicate to be more reliable.
    internal inner class abs_comparison : Comparator<Pair<Vertex, Double>> {
        override fun compare(lhs: Pair<Vertex, Double>, rhs: Pair<Vertex, Double>): Int {
            return -abs(lhs.second).compareTo(abs(rhs.second))
        }
    }

    // data required for adding a new edge
    ///
    // used in add_edge() for storing information related to
    // the new edge.
    internal inner class EdgeData {
        lateinit var v1_prv: Edge  // < edge prior to v1
        lateinit var v1: Vertex  // < NEW edge source
        lateinit var v1_nxt: Edge // < edge following v1
        lateinit var v2_prv: Edge // < edge prior to v2
        lateinit var v2: Vertex // < NEW edge target
        lateinit var v2_nxt: Edge // < edge following v2
        lateinit var f: Face // < face of v1 and v2
        override fun toString(): String {
            return "EdgeData(\n  v1_prv: $v1_nxt\n  v1: $v1\n  v1_nxt: $v1_nxt\n  v2_prv: $v2_prv\n  v2: $v2\n  v2_nxt: $v2_nxt\n)"
        }
    }

    /**
     * Initializes the diagram with three generators. Add one vertex at origin and
     * three vertices at 'infinity' and their associated edges
     */
    private fun initialize() {
        val far_multiplier = 6.0
        // initial generators/sites:
        val gen1 = Point(0.0, 3.0 * far_radius)
        val gen2 = Point(-3.0 * sqrt(3.0) * far_radius / 2.0, -3.0 * far_radius / 2.0)
        val gen3 = Point(+3.0 * sqrt(3.0) * far_radius / 2.0, -3.0 * far_radius / 2.0)
        // initial vd-vertices
        val vd1 = Point(0.0, -3.0 * far_radius * far_multiplier)
        val vd2 = Point(
            +3.0 * sqrt(3.0) * far_radius * far_multiplier / 2.0,
            +3.0 * far_radius * far_multiplier / 2.0
        )
        val vd3 = Point(
            -3.0 * sqrt(3.0) * far_radius * far_multiplier / 2.0,
            +3.0 * far_radius * far_multiplier / 2.0
        )
        // add init vertices
        val v00 = g.add_vertex(Vertex(Point(0.0, 0.0), VertexStatus.UNDECIDED, VertexType.NORMAL, gen1))
        val v01 = g.add_vertex(Vertex(vd1, VertexStatus.OUT, VertexType.OUTER, gen3))
        val v02 = g.add_vertex(Vertex(vd2, VertexStatus.OUT, VertexType.OUTER, gen1))
        val v03 = g.add_vertex(Vertex(vd3, VertexStatus.OUT, VertexType.OUTER, gen2))
        // add initial sites to graph
        val vert1 = g.add_vertex(Vertex(gen1, VertexStatus.OUT, VertexType.POINTSITE))
        val vert2 = g.add_vertex(Vertex(gen2, VertexStatus.OUT, VertexType.POINTSITE))
        val vert3 = g.add_vertex(Vertex(gen3, VertexStatus.OUT, VertexType.POINTSITE))

        // apex-points on the three edges:
        val a1 = g.add_vertex(Vertex(gen2.add(gen3).mult(0.5), VertexStatus.UNDECIDED, VertexType.APEX, gen2))
        val a2 = g.add_vertex(Vertex(gen1.add(gen3).mult(0.5), VertexStatus.UNDECIDED, VertexType.APEX, gen3))
        val a3 = g.add_vertex(Vertex(gen1.add(gen2).mult(0.5), VertexStatus.UNDECIDED, VertexType.APEX, gen1))

        // add face 1: v0-v1-v2 which encloses gen3
        val e1_1 = g.add_edge(v00, a1)
        val e1_2 = g.add_edge(a1, v01)
        val e2 = g.add_edge(v01, v02)
        val e3_1 = g.add_edge(v02, a2)
        val e3_2 = g.add_edge(a2, v00)
        val f1 = g.add_face()
        f1.site = PointSite(gen3, f1, vert3)
        f1.status = FaceStatus.NONINCIDENT
        kd_tree.addPoint(doubleArrayOf(gen3.x, gen3.y), KdPoint(gen3, f1))
        g.set_next_cycle(listOf(e1_1, e1_2, e2, e3_1, e3_2), f1, 1.0)

        // add face 2: v0-v02-v03 which encloses gen1
        val e4_1 = g.add_edge(v00, a2)
        val e4_2 = g.add_edge(a2, v02)
        val e5 = g.add_edge(v02, v03)
        val e6_1 = g.add_edge(v03, a3)
        val e6_2 = g.add_edge(a3, v00)
        val f2 = g.add_face()
        f2.site = PointSite(gen1, f2, vert1)
        f2.status = FaceStatus.NONINCIDENT
        kd_tree.addPoint(doubleArrayOf(gen1.x, gen1.y), KdPoint(gen1, f2))
        g.set_next_cycle(listOf(e4_1, e4_2, e5, e6_1, e6_2), f2, 1.0)

        // add face 3: v0-v3-v1 which encloses gen2
        val e7_1 = g.add_edge(v00, a3)
        val e7_2 = g.add_edge(a3, v03)
        val e8 = g.add_edge(v03, v01)
        val e9_1 = g.add_edge(v01, a1)
        val e9_2 = g.add_edge(a1, v00)
        val f3 = g.add_face()
        f3.site = PointSite(gen2, f3, vert2) // this constructor needs f3...
        f3.status = FaceStatus.NONINCIDENT
        kd_tree.addPoint(doubleArrayOf(gen2.x, gen2.y), KdPoint(gen2, f3))
        g.set_next_cycle(listOf(e7_1, e7_2, e8, e9_1, e9_2), f3, 1.0)

        // set type.
        e1_1.type = EdgeType.LINE
        e1_1.set_parameters(f1.site, f3.site, false)
        e1_2.type = EdgeType.LINE
        e1_2.set_parameters(f1.site, f3.site, true)
        e2.type = EdgeType.OUTEDGE
        e3_1.type = EdgeType.LINE
        e3_1.set_parameters(f2.site, f1.site, true)
        e3_2.type = EdgeType.LINE
        e3_2.set_parameters(f2.site, f1.site, false)
        e4_1.type = EdgeType.LINE
        e4_1.set_parameters(f2.site, f1.site, false)
        e4_2.type = EdgeType.LINE
        e4_2.set_parameters(f2.site, f1.site, true)
        e5.type = EdgeType.OUTEDGE
        e6_1.type = EdgeType.LINE
        e6_1.set_parameters(f2.site, f3.site, false)
        e6_2.type = EdgeType.LINE
        e6_2.set_parameters(f2.site, f3.site, true)
        e7_1.type = EdgeType.LINE
        e7_1.set_parameters(f2.site, f3.site, true)
        e7_2.type = EdgeType.LINE
        e7_2.set_parameters(f2.site, f3.site, false)
        e8.type = EdgeType.OUTEDGE
        e9_1.type = EdgeType.LINE
        e9_1.set_parameters(f1.site, f3.site, true)
        e9_2.type = EdgeType.LINE
        e9_2.set_parameters(f1.site, f3.site, false)

        // twin edges
        g.twin_edges(e1_1, e9_2)
        g.twin_edges(e1_2, e9_1)
        e2.twin = null // the outermost edges have invalid twins
        e5.twin = null
        e8.twin = null
        g.twin_edges(e3_1, e4_2)
        g.twin_edges(e3_2, e4_1)
        g.twin_edges(e6_1, e7_2)
        g.twin_edges(e6_2, e7_1)
        assert(vd_checker.is_valid()) { " vd_checker.is_valid() " }
    }

    /**
     * find amount of clearance-disk violation on all vertices of face f
     *
     * @return vertex with the largest clearance-disk violation
     */
    private fun find_seed_vertex(f: Face, site: Site): Vertex {
        var minPred = 0.0
        var minimalVertex: Vertex? = null
        var first = true
        var current = f.edge
        val start = current
        do {
            val q = current.target
            if (q.status !== VertexStatus.OUT && q.type === VertexType.NORMAL) {
                val h = q.in_circle(site.apex_point(q.position))
                if (first || h < minPred && site.in_region(q.position)) {
                    minPred = h
                    minimalVertex = q
                    first = false
                }
            }
            current = current.next
        } while (!current.equals(start))
        assert(minPred < 0) { " minPred < 0 " }
        return minimalVertex!!
    }

    // find and return IN - OUT edges
    ///
    // given the set v0 of IN vertices, find and return the adjacent IN - OUT
    // edges.
    // Later NEW vertices are inserted into each of the found IN - OUT edges
    private fun find_in_out_edges(): List<Edge> {
        assert(!v0.isEmpty()) { " !v0.isEmpty() " }
        val output: MutableList<Edge> = mutableListOf() // new vertices generated on these edges
        for (v in v0) {
            assert(v.status === VertexStatus.IN) { " v.status == VertexStatus.IN " }
            for (e in v.out_edges) {
                if (e.target.status === VertexStatus.OUT) {
                    output.add(e) // this is an IN-OUT edge
                }
            }
        }
        assert(!output.isEmpty()) { " !output.isEmpty() " }
        return output
    }

    // find EdgeData for a new edge
    ///
    // on a face which has IN and OUT-vertices, find the sequence
    // OUT-OUT-OUT-..-OUT-NEW(v1)-IN-...-IN-NEW(v2)-OUT-OUT
    // and return v1/v2 together with their previous and next edges
    // \param f face on which we search for vertices
    // \param startverts contains NEW-vertices already found, which are not valid
    // for this call to find_edge_data
    // \param segment contains ENDPOINT vertices, when we are inserting a
    // line-segment
    // (these vertices are needed to ensure finding correct points around
    // sites/null-edges)
    private fun find_edge_data(f: Face, startverts: List<Vertex?>, segment: Pair<Vertex?, Vertex?>): EdgeData {
        val ed: EdgeData = EdgeData()
        ed.f = f
        var current_edge = f.edge // start on some edge of the face
        var start_edge = current_edge
        var found = false
        do { // find OUT-NEW-IN vertices in this loop
            val next_edge = current_edge.next
            val previous_vertex = current_edge.source
            val current_vertex = current_edge.target
            val next_vertex = next_edge.target
            val previous_not_endpoint = (!previous_vertex.equals(segment.first)
                    && !previous_vertex.equals(segment.second))
            val next_is_endpoint = next_vertex.equals(segment.first) || next_vertex.equals(segment.second)
            if (current_vertex.status === VertexStatus.NEW && current_vertex.type !== VertexType.SEPPOINT
                && ((previous_vertex.status === VertexStatus.OUT
                        || previous_vertex.status === VertexStatus.UNDECIDED) && previous_not_endpoint
                        || next_is_endpoint)
            ) {
                // slow? linear search through vector. but startverts.size() should not be too
                // large..
                val v_in_startverts = startverts.contains(current_vertex)
                if (!v_in_startverts) {
                    ed.v1 = current_vertex
                    ed.v1_prv = current_edge
                    ed.v1_nxt = current_edge.next
                    found = true
                }
            }
            current_edge = current_edge.next
        } while (!current_edge.equals(start_edge) && !found)
        assert(found) { "found" }

        // now search for v2
        start_edge = current_edge // note that this search starts where we ended in the loop above!
        found = false
        do { // find IN-NEW-OUT vertices in this loop
            val current_vertex = current_edge.target
            if (current_vertex.status === VertexStatus.NEW && current_vertex.type !== VertexType.SEPPOINT) {
                if (!current_vertex.equals(ed.v1)) { // -IN-NEW(v2)
                    ed.v2 = current_vertex
                    ed.v2_prv = current_edge
                    ed.v2_nxt = current_edge.next
                    found = true
                }
            }
            current_edge = current_edge.next
        } while (!current_edge.equals(start_edge) && !found)
        assert(found) { "found" }
        return ed
    }

    // find and return edges on which we potentially need SPLIT vertices
    ///
    // walk around the face f
    // return edges whose endpoints are on separate sides of pt1-pt2 line
    // \todo ?not all edges found like this *need* SPLIT vertices? (but it does not
    // hurt to insert SPLIT-vertices in this case)
    private fun find_split_edges(f: Face, pt1: Point, pt2: Point): List<Edge> {
        assert(vd_checker.face_ok(f)) { " vd_checker.face_ok(f) " }
        val out: MutableList<Edge> = mutableListOf()
        var current_edge = f.edge
        val start_edge = current_edge
        // int count=0;
        do { // FIND ALL! not just one.
            val trg = current_edge.target
            val src = current_edge.source
            val src_is_right = src.position.is_right(pt1, pt2)
            val trg_is_right = trg.position.is_right(pt1, pt2)
            if (src.type === VertexType.NORMAL || src.type === VertexType.APEX || src.type === VertexType.SPLIT) { // ?
                // check
                // edge-type
                // instead?
                if (src_is_right != trg_is_right) {
                    out.add(current_edge)
                    assert(vd_checker.check_edge(current_edge)) { "vd_checker.check_edge(current_edge)" }
                }
            }
            current_edge = current_edge.next
        } while (!current_edge.equals(start_edge))
        return out
    }

    // find a SPLIT vertex on the Face f
    // return true, and set v, if found.
    private fun find_split_vertex(f: Face): Vertex? {
        for (q in g.face_vertices(f)) {
            if (q.type === VertexType.SPLIT) {
                return q
            }
        }
        return null
    }

    // find an adjacent vertex to v, along an edge that is not a NULLEDGE
    // return true if found, otherwise false.
    private fun null_vertex_target(v: Vertex): Vertex? {
        for (e in v.out_edges) {
            if (e.type !== EdgeType.NULLEDGE) {
                return e.target
            }
        }
        return null
    }

    /**
     * Grows the delete-tree of IN vertices by "weighted breadth-first search"
     *
     *
     * we start at the seed and add vertices with detH<0 provided that: - (C4) v
     * should not be adjacent to two or more IN vertices (this would result in a
     * loop/cycle!) - (C5) for an incident face containing v: v is adjacent to an IN
     * vertex on this face
     *
     *
     * C4 and C5 refer to the Sugihara&Iri 1992 "one million" paper. We process
     * UNDECIDED vertices adjacent to known IN-vertices in a "weighted
     * breadth-first-search" manner where vertices with a large fabs(detH) are
     * processed first, since we assume the in-circle predicate to be more reliable
     * the larger fabs(in_circle()) is.
     *
     * @param site
     */
    private fun augment_vertex_set(site: Site) {
        while (!vertexQueue.isEmpty()) {
            val (v, h) = vertexQueue.poll()
            assert(v.status === VertexStatus.UNDECIDED) { " v.status == VertexStatus.UNDECIDED " }
            if (h < 0.0) { // try to mark IN if h<0 and passes (C4) and (C5) tests and in_region().
                // otherwise mark OUT
                if (predicate_c4(v) || !predicate_c5(v) || !site.in_region(v.position)) {
                    v.status = VertexStatus.OUT // C4 or C5 violated, so mark OUT
                } else {
                    mark_vertex(v, site) // h<0 and no violations, so mark IN. push adjacent UNDECIDED vertices onto Q.
                }
            } else {
                v.status = VertexStatus.OUT // detH was positive (or zero), so mark OUT
            }
            modified_vertices.add(v)
        }
        assert(vertexQueue.isEmpty()) { " vertexQueue.isEmpty() " }
        assert(vd_checker.all_in(v0)) { " vd_checker.all_in(v0) " }
        // sanity-check?: for all incident faces the IN/OUT-vertices should be connected
    }

    // adjacent in-count predicate for buildingdelete-tree
    ///
    // number of IN vertices adjacent to given vertex v
    // predicate C4 i.e. "adjacent in-count" from Sugihara&Iri 1992 "one million"
    // paper
    private fun predicate_c4(v: Vertex): Boolean {
        var in_count = 0
        for (e in v.out_edges) {
            val w = e.target
            if (w.status === VertexStatus.IN) {
                in_count++
                if (in_count >= 2) {
                    return true
                }
            }
        }
        return false
    }

    // connectedness-predicate for delete-tree building
    ///
    // do any of the three faces that are adjacent to the given IN-vertex v have an
    // IN-vertex ?
    // predicate C5 i.e. "connectedness" from Sugihara&Iri 1992 "one million" paper
    private fun predicate_c5(v: Vertex): Boolean {
        if (v.type === VertexType.APEX || v.type === VertexType.SPLIT) {
            return true
        } // ?
        val adjacent_incident_faces: MutableList<Face> = mutableListOf()
        for (e in v.out_edges) {
            if (e.face.status === FaceStatus.INCIDENT) {
                adjacent_incident_faces.add(e.face)
            }
        }
        assert(!adjacent_incident_faces.isEmpty()) { " !adjacent_incident_faces.isEmpty() " }
        for (f in adjacent_incident_faces) {
            // check each adjacent face f for an IN-vertex
            var face_ok = false
            var current = f.edge
            val start = current
            do {
                val w = current.target
                if (!w.equals(v) && w.status === VertexStatus.IN && g.has_edge(w, v)) { // v should be adjacent to an IN
                    // vertex on the face
                    face_ok = true
                } else if (!w.equals(v)
                    && (w.type === VertexType.ENDPOINT || w.type === VertexType.APEX || w.type === VertexType.SPLIT)
                ) { // if
                    // we
                    // are
                    // next
                    // to
                    // an
                    // ENDPOINT,
                    // then
                    // ok(?)
                    face_ok = true
                } else if (!w.equals(v) && w.type === VertexType.SEPPOINT && g.has_edge(w, v)) {
                    face_ok = true
                }
                current = current.next
            } while (!current.equals(start))
            if (!face_ok) {
                return false
            }
        }
        return true // if we get here we found all ok
    }

    // mark adjacent faces INCIDENT
    // call this when inserting line-sites
    // since we call add_split_vertex we can't use iterators, because they get
    // invalidated
    // so use the slower adjacent_faces() instead.
    private fun mark_adjacent_faces(v: Vertex, site: Site) {
        assert(v.status === VertexStatus.IN) { "v.status == VertexStatus.IN " }
        val new_adjacent_faces = g.adjacent_faces(v)
        assert(v.type === VertexType.APEX && new_adjacent_faces.size == 2 || v.type === VertexType.SPLIT && new_adjacent_faces.size == 2 || new_adjacent_faces.size == 3)
        for (adj_face in new_adjacent_faces) {
            if (adj_face.status !== FaceStatus.INCIDENT) {
                if (site.isLine) {
                    add_split_vertex(adj_face, site)
                }
                adj_face.status = FaceStatus.INCIDENT
                incident_faces.add(adj_face)
            }
        }
    }

    // mark adjacent faces INCIDENT
    // IN-Vertex v has three adjacent faces, mark nonincident faces incident
    // and push them to the incident_faces queue
    // NOTE: call this only when inserting point-sites
    private fun mark_adjacent_faces_p(v: Vertex) {
        assert(v.status === VertexStatus.IN) { "v.status == VertexStatus.IN " }
        for (e in v.out_edges) {
            val adj_face = e.face
            if (adj_face.status !== FaceStatus.INCIDENT) {
                adj_face.status = FaceStatus.INCIDENT
                incident_faces.add(adj_face)
            }
        }
    }

    /**
     * mark vertex IN and mark adjacent faces INCIDENT; push adjacent UNDECIDED
     * vertices onto queue
     */
    private fun mark_vertex(v: Vertex, site: Site) {
        v.status = VertexStatus.IN
        v0.add(v)
        modified_vertices.add(v)
        if (site.isPoint) {
            mark_adjacent_faces_p(v)
        } else {
            mark_adjacent_faces(v, site)
        }

        // push the v-adjacent vertices onto the queue
        for (e in v.out_edges) {
            val w = e.target
            if (w.status === VertexStatus.UNDECIDED && !w.in_queue) {
                // when pushing onto queue we also evaluate in_circle predicate so that we
                // process vertices in the correct order
                vertexQueue.add(Pair(w, w.in_circle(site.apex_point(w.position))))
                w.in_queue = true
            }
        }
    }

    /**
     * Add NEW vertices on IN-OUT edges (generates new voronoi-vertices on all
     * IN-OUT edges)
     *
     * @implNote used only by [.insert_point_site]
     * @param new_site
     */
    private fun add_vertices(new_site: Site) {
        assert(!v0.isEmpty()) { " !v0.isEmpty() " }
        val q_edges: List<Edge> = find_in_out_edges() // new vertices generated on these IN-OUT edges
        for (e in q_edges) {
            val sl = vpos.position(e, new_site) // vertex_positioner.cpp
            val q = g.add_vertex(
                Vertex(sl.p, VertexStatus.NEW, VertexType.NORMAL, new_site.apex_point(sl.p), sl.k3)
            )
            modified_vertices.add(q)
            g.add_vertex_in_edge(q, e)
            q.max_error = vpos.dist_error(e, sl, new_site)
        }
    }

    // add a new face corresponding to the new Site
    ///
    // call add_new_edge() on all the incident_faces that should be split
    private fun add_face(s: Site): Face {
        val newface = g.add_face()
        newface.site = s
        s.face = newface
        newface.status = FaceStatus.NONINCIDENT
        if (s.isPoint) {
            kd_tree.addPoint(doubleArrayOf(s.position().x, s.position().y), KdPoint(s.position(), newface))
        }
        return newface
    }

    private fun add_edges(newface: Face, f: Face) {
        add_edges(newface, f, null, Pair(null, null))
    }

    // add all NEW-NEW edges
    ///
    // by adding a NEW-NEW edge, split the face f into one part which is newface,
    // and the other part is the old f
    // for linesegment or arc sites we pass in both the k=+1 face newface and the
    // k=-1 face newface2
    // the segment endpoints are passed to find_edge_data()
    private fun add_edges(newface: Face, f: Face, newface2: Face?, segment: Pair<Vertex?, Vertex?>) {
        val new_count = num_new_vertices(f)
        assert(new_count > 0) { " new_count > 0 " }
        assert(new_count % 2 == 0) { " (new_count % 2) == 0 " }
        val new_pairs = new_count / 2 // we add one NEW-NEW edge for each pair found
        val startverts: MutableList<Vertex?> = mutableListOf() // this holds ed.v1 vertices for edges already added
        for (m in 0 until new_pairs) {
            val ed = find_edge_data(f, startverts, segment)
            add_edge(ed, newface, newface2)
            startverts.add(ed.v1)
        }
    }

    // add a new edge to the diagram
    // newface = the k=+1 positive offset face
    // newface2 = the k=-1 negative offset face
    private fun add_edge(ed: EdgeData, newface: Face, newface2: Face?) {
        val new_previous: Edge = ed.v1_prv
        val new_source: Vertex = ed.v1 // -OUT-NEW(v1)-IN-...
        val twin_next: Edge = ed.v1_nxt
        val twin_previous: Edge = ed.v2_prv
        val new_target: Vertex = ed.v2 // -IN-NEW(v2)-OUT-
        val new_next: Edge = ed.v2_nxt
        val f: Face = ed.f
        val f_site = f.site
        val new_site: Site
        val new_face: Face
        if (new_source.k3 == 1.0) { // find out if newface or newface2 should be used
            new_site = newface.site
            new_face = newface
        } else {
            new_site = newface2!!.site
            new_face = newface2
        }
        assert(new_target.k3 == new_source.k3) { " new_target.k3 == new_source.k3 " }

        // f
        // now connect: new_previous -> new_source -> new_target -> new_next
        // and: twin_next <- new_source <- new_target <- twin_previous
        // new_face

        // check for potential apex-split
        // we need an apex-vertex if the source and target are on different branches of
        // the new quadratic edge
        // we can set the src_sign and trg_sign with is_right where we compare to a line
        // through the apex
        var src_sign = true
        var trg_sign = true
        if (f_site.isPoint && new_site.isLine) { // PL or PA
            val pt2 = f_site.position()
            val pt1 = new_site.apex_point(pt2) // projection of pt2 onto LineSite or ArcSite
            src_sign = new_source.position.is_right(pt1, pt2)
            trg_sign = new_target.position.is_right(pt1, pt2)
        } else if (f_site.isPoint && new_site.isArc) {
            val pt2 = f_site.position()
            // project p2 onto circle
            val cen = Point(new_site.x(), new_site.y())
            val cen_pt2 = pt2.sub(cen)
            val pt1 = cen.add(cen_pt2.mult(new_site.r() / cen_pt2.norm()))
            src_sign = new_source.position.is_right(pt1, pt2)
            trg_sign = new_target.position.is_right(pt1, pt2)
        } else if (f_site.isPoint && new_site.isPoint) { // PP
            src_sign = new_source.position.is_right(f_site.position(), new_site.position())
            trg_sign = new_target.position.is_right(f_site.position(), new_site.position())
        } else if (f_site.isLine && new_site.isLine) { // LL
            // a line-line bisector, sign should not matter because there is no sqrt()
            // this is essentially an in-region test
            if (new_source.position !== new_target.position && new_source.position !== f_site.start() && new_source.position !== f_site.end() && new_target.position !== f_site.start() && new_target.position !== f_site.end() && new_source.position.sub(
                    f_site.apex_point(new_source.position)
                ).norm() > 1e-3 && new_target.position.sub(f_site.apex_point(new_target.position))
                    .norm() > 1e-3 // so that the
            // is_right
            // predicate is
            // accurate
            ) {
                assert(
                    !new_source.position.is_right(
                        f_site.start(),
                        f_site.end()
                    )
                ) { " !new_source.position.is_right( f_site.start(), f_site.end() ) " }
                assert(
                    !new_target.position.is_right(
                        f_site.start(),
                        f_site.end()
                    )
                ) { " !new_target.position.is_right( f_site.start(), f_site.end() ) " }
                assert(
                    !new_source.position.is_right(
                        new_site.start(),
                        new_site.end()
                    )
                ) { " !new_source.position.is_right( new_site.start(), new_site.end() ) " }
                assert(
                    !new_target.position.is_right(
                        new_site.start(),
                        new_site.end()
                    )
                ) { " !new_target.position.is_right( new_site.start(), new_site.end() ) " }
            }
        } else if (f_site.isLine && new_site.isArc) { // LA
            val pt2 = Point(new_site.x(), new_site.y())
            val pt1 = f_site.apex_point(pt2)
            src_sign = new_source.position.is_right(pt1, pt2)
            trg_sign = new_target.position.is_right(pt1, pt2)
            // if one vertex is on a null-face, we cannot trust the sign
            if (new_source.dist() == 0.0 || new_target.dist() == 0.0) {
                if (new_source.dist() > new_target.dist()) {
                    src_sign = trg_sign
                } else {
                    trg_sign = src_sign
                }
            }
        } else { // unhandled case!
            throw RuntimeException(" add_edge() WARNING: no code to deremine src_sign and trg_sign!")
        }

        // both src and trg are on the same side of the new site.
        // so no apex-split is required, just add a single edge.
        if (src_sign == trg_sign) { // add a single src-trg edge
            val (e_new, e_twin) = g.add_twin_edges(new_source, new_target)
            e_new.next = new_next
            assert(new_next.k == new_previous.k) { " new_next.k == new_previous.k " }
            e_new.k = new_next.k // the next edge is on the same face, so has the correct k-value
            e_new.face = f // src-trg edge has f on its left
            new_previous.next = e_new
            f.edge = e_new
            e_new.set_parameters(f_site, new_site, !src_sign)
            twin_previous.next = e_twin
            e_twin.next = twin_next
            e_twin.k = new_source.k3
            e_twin.set_parameters(f_site, new_site, !src_sign) // new_site, f_site, src_sign
            e_twin.face = new_face
            new_face.edge = e_twin
            assert(vd_checker.check_edge(e_new) && vd_checker.check_edge(e_twin)) { " vd_checker.check_edge(e_new) && vd_checker.check_edge(e_twin) " }
        } else {
            // need to do apex-split, and add two new edges
            // f f
            // new_prv -> NEW -- e1 ----> APEX --e2 ---> NEW -> new_nxt
            // twn_nxt <- NEW <- e1_tw -- APEX <-e2_tw-- NEW <- twn_prv
            // new1/new2 new1/new2
            //
            val apex = g.add_vertex(Vertex(Point(0.0, 0.0), VertexStatus.NEW, VertexType.APEX))
            val (e1, e1_tw) = g.add_twin_edges(new_source, apex)
            val (e2, e2_tw) = g.add_twin_edges(apex, new_target)
            e1.set_parameters(f_site, new_site, !src_sign)
            e2.set_parameters(f_site, new_site, !trg_sign)
            assert(new_previous.face === f) { " new_previous.face == f " }
            assert(new_next.face === new_previous.face) { " new_next.face == new_previous.face " }
            assert(new_next.k == new_previous.k) { " new_next.k == new_previous.k " }

            // new_previous -> e1 -> e2 -> new_next
            new_previous.next = e1
            e1.next = e2
            e2.next = new_next
            e1.face = f
            e2.face = f
            e1.k = new_next.k
            e2.k = new_next.k
            f.edge = e1
            // twin edges
            e1_tw.set_parameters(new_site, f_site, src_sign)
            e2_tw.set_parameters(new_site, f_site, trg_sign)
            assert(twin_previous.k == twin_next.k) { " twin_previous.k == twin_next.k " }
            assert(twin_previous.face === twin_next.face) { " twin_previous.face == twin_next.face " }
            // twin_prev -> e2_tw -> e1_tw -> twin_next on new_face
            twin_previous.next = e2_tw
            e2_tw.next = e1_tw
            e1_tw.next = twin_next
            e1_tw.k = new_source.k3
            e2_tw.k = new_source.k3
            new_face.edge = e1_tw
            e1_tw.face = new_face
            e2_tw.face = new_face
            assert(vd_checker.check_edge(e1) && vd_checker.check_edge(e1_tw)) { " vd_checker.check_edge(e1) && vd_checker.check_edge(e1_tw) " }
            assert(vd_checker.check_edge(e2) && vd_checker.check_edge(e2_tw)) { " vd_checker.check_edge(e2) && vd_checker.check_edge(e2_tw) " }

            // position the apex
            val min_t = e1.minimum_t(f_site, new_site)
            apex.position = e1.point(min_t)
            apex.init_dist(f_site.apex_point(apex.position))
            modified_vertices.add(apex)
        }
    }

    internal inner class SeparatorTarget(v_previous: Edge?, v_target: Vertex?, v_next: Edge?, out_new_in: Boolean) {
        var v_previous: Edge?
        var v_target: Vertex?
        var v_next: Edge?
        var out_new_in: Boolean

        init {
            this.v_previous = v_previous
            this.v_target = v_target
            this.v_next = v_next
            this.out_new_in = out_new_in
        }

        override fun toString(): String {
            return "SeparatorTarget(${v_previous}, ${v_target}, ${v_next}, ${out_new_in})"
        }
    }

    /**
     * Add SEPARATOR edge on the face f, which contains the endpoint
     *
     * @param f         face of endpoint
     * @param null_face null face of the endpoint
     * @param target    target-data found by ??
     * @param sep_endp
     * @param s1        positive LineSite
     * @param s2        negative LineSite
     */
    private fun add_separator(
        f: Face,
        null_face: Face,
        target: SeparatorTarget,
        sep_endp: Vertex?,
        s1: Site,
        s2: Site
    ) {
        if (sep_endp == null) {
            return  // do nothing!
        }
        assert(sep_endp.k3 == 1.0 || sep_endp.k3 == -1.0) { " (sep_endp.k3==1) || (sep_endp.k3==-1) " }
        sep_endp.zero_dist()
        val (endp_next_tw, endp_prev_tw) = g.find_next_prev(null_face, sep_endp)
        val endp_prev = endp_next_tw.twin!! // NOTE twin!
        val endp_next = endp_prev_tw.twin!! // NOTE twin!
//        assert(endp_next != null) { " endp_next != null " }
//        assert(endp_prev != null) { " endp_prev != null " }

        // find NEW vertex on the old face f
        // this vertex has the correct alfa angle for this endp/separator
        val v_previous = target.v_previous!!
        val v_target = target.v_target!!
        val v_next = target.v_next!!
        val out_new_in = target.out_new_in
        assert(v_target.k3 == 1.0 || v_target.k3 == -1.0) { " (v_target.k3==1) || (v_target.k3==-1) " }
        assert(sep_endp.k3 == v_target.k3) { " sep_endp.k3 == v_target.k3 " }
        // can't assert about in_region - numerical error is always present
        // assert( s1.in_region(v_target.position ) ) : " s1.in_region(v_target.position
        // ) ";
        // assert( s2.in_region(v_target.position ) ) : " s2.in_region(v_target.position
        // ) ";

        // add new separator edge, and its twin
        val (e2, e2_tw) = g.add_twin_edges(sep_endp, v_target)
        e2.type = EdgeType.SEPARATOR
        e2_tw.type = EdgeType.SEPARATOR

        // there are two cases. depending on how v_target (NEW) is found:
        // OUT-NEW-IN, when out_new_in = true
        // IN-NEW-OUT, when out_new_in = false
        // here we set the faces, sites, and next-pointers depending on the case
        if (out_new_in) {
            e2.k = v_target.k3 // e2 is on the segment side
            e2_tw.k = +1.0 // e2_tw is on the point-site side
            e2_tw.face = f // point-site face
            e2_tw.null_face = f
            e2_tw.has_null_face = true
            f.edge = e2_tw
            endp_prev.k = e2.k // endp_prev is on the line-site side
            if (e2.k == -1.0) { // choose either s1 or s2 as the site
                e2.face = s2.face
                s2.face.edge = e2
                endp_prev.face = s2.face
            } else {
                e2.face = s1.face
                s1.face.edge = e2
                endp_prev.face = s1.face
            }
            g.set_next(v_previous, e2_tw)
            g.set_next(e2_tw, endp_next)
            endp_next.face = f // the null-edge
            endp_next.k = 1.0
            g.set_next(e2, v_next)
        } else {
            e2.k = +1.0 // e2 is on the point-site side
            e2_tw.k = v_target.k3 // e2_tw is on the segment side
            e2.face = f // point-site face
            e2.null_face = f
            e2.has_null_face = true
            f.edge = e2
            endp_next.k = e2_tw.k // endp_next is on the linesite-side
            if (e2_tw.k == -1.0) {
                e2_tw.face = s2.face
                s2.face.edge = e2_tw
                endp_next.face = s2.face
            } else {
                e2_tw.face = s1.face
                s1.face.edge = e2_tw
                endp_next.face = s1.face
            }
            g.set_next(v_previous, e2_tw)
            endp_prev.face = f
            endp_prev.k = 1.0
            g.set_next(endp_prev, e2)
            g.set_next(e2, v_next)
        }
        e2.set_sep_parameters(sep_endp.position, v_target.position)
        e2_tw.set_sep_parameters(sep_endp.position, v_target.position)
        assert(vd_checker.check_edge(e2)) { " vd_checker.check_edge(e2) " }
        assert(vd_checker.check_edge(e2_tw)) { " vd_checker.check_edge(e2_tw) " }
    }

    /**
     * Adds one or many SPLIT vertices to the edges of the give face. These are
     * projections/mirrors of the site of f with the new Site s acting as the
     * mirror. SPLIT vertices are inserted to avoid deleting loops during //
     * augment_vertex_set().
     *
     * @param f
     * @param s
     */
    private fun add_split_vertex(f: Face, s: Site) {
        assert(!s.isPoint) { "!s.isPoint" }
        val fs = f.site

        // don't search for split-vertex on the start or end face
        if (fs.isPoint && s.isLine) {
            if (fs.position() === s.start() || fs.position() === s.end()) {
                return
            }
        }
        if (fs.isPoint && s.isLine && s.in_region(fs.position())) {
            // 1) find the correct edge
            val pt1 = fs.position()
            val pt2 = pt1.sub(Point(s.a(), s.b()))
            assert(pt1.sub(pt2).norm() > 0) { " (pt1.sub(pt2)).norm() > 0 " }
            val split_edges: List<Edge> = find_split_edges(f, pt1, pt2)
            // the sought edge should have src on one side of pt1-pt2
            // and trg on the other side of pt1-pt2
            for (split_edge in split_edges) {
                if (split_edge.type === EdgeType.SEPARATOR || split_edge.type === EdgeType.LINESITE) {
                    return  // don't place split points on linesites or separators(?)
                    // find a point = src + u*(trg-src)
                    // with min_t < u < max_t
                    // and minimum distance to the pt1-pt2 line
                }
                var split_pt_pos: Point
                val split_src = split_edge.source
                val split_trg = split_edge.target
                val errFunctr = SplitPointError(g, split_edge, pt1, pt2) // error functor
                val min_t = min(split_src.dist(), split_trg.dist())
                val max_t = max(split_src.dist(), split_trg.dist())
                // require that min_t and max_t bracket the root
                if (errFunctr.value(min_t) * errFunctr.value(max_t) >= 0) {
                    return
                }

                val result = toms748(min_t, max_t, 1E-15, 100, errFunctr::value).x

                split_pt_pos = split_edge.point(result)
                val v =
                    g.add_vertex(Vertex(split_pt_pos, VertexStatus.UNDECIDED, VertexType.SPLIT, fs.position()))
                assert(vd_checker.check_edge(split_edge)) { " vd_checker.check_edge(split_edge) " }
                // 3) insert new SPLIT vertex into the edge
                g.add_vertex_in_edge(v, split_edge)
            }
        }
    }

    /**
     * Utility container
     */
    private inner class FiveTuple1(v1: Vertex?, f1: Face?, v2: Vertex?, v3: Vertex?, f2: Face?) {
        var get1: Vertex?
        var get2: Face?
        var get3: Vertex?
        var get4: Vertex?
        var get5: Face?

        init {
            get1 = v1
            get2 = f1
            get3 = v2
            get4 = v3
            get5 = f2
        }
    }

    // either find an existing null-face, or create a new one.
    // \param start the end of the segment for which we will find/create a
    // null-face
    // \param other the other end of the new segment
    // \param left a point left of the new segment
    // \param dir alfa-direction for positioning endpoint vertex on null-face
    // \param new_site the new Site we are inserting
    ///
    // \return HEVertex ENDPOINT-vertex for the new vertex
    // \return HEFace null-face at endpoint (new or existing)
    // \return HEVertex positive SEPARATOR edge endpoint vertex (if a positive
    // separator should be added)
    // \return HEVertex negative SEPARATOR edge endpoint vertex (if a negative
    // separator should be added)
    // \return HEFace face-to-null. if a PointSite face should disappear, we return
    // it here.
    private fun find_null_face(start: Vertex, other: Vertex, left: Point, dir: Point, new_site: Site): FiveTuple1 {
        var seg_start: Vertex?  // new end-point vertices
        var start_null_face: Face?  // either existing or new null-face at start-vertex
        var pos_sep_start: Vertex?  // optional separator endpoints at start
        var neg_sep_start: Vertex?  // invalid vertices are default
        var face_to_null: Face? = null // invalid face is default

        // this works for LineSite
        var k3_sign: Boolean = if (new_site.isLine) {
            left.is_right(start.position, other.position)
        } else if (new_site.isArc) {
            Point(new_site.x(), new_site.y()).is_right(start.position, other.position)
        } else {
            throw RuntimeException()
        }
        // this is used below and in find_null_face()
        // k3_sign is already calculated in insert_line_segment() ??
        return if (start.null_face != null) { // there is an existing null face
            start_null_face = start.null_face

            // create a new segment ENDPOINT vertex with zero clearance-disk
            seg_start = g.add_vertex(Vertex(start.position, VertexStatus.OUT, VertexType.ENDPOINT, 0.0))
            // find the edge on the null-face where we insert seg_start
            var insert_edge: Edge? = null
            var current2 = start_null_face!!.edge
            val start_edge2 = current2
            seg_start.set_alfa(dir)
            var found = false
            do {
                val face_incident = current2.twin!!.face.status === FaceStatus.INCIDENT
                if (face_incident) { // pick any incident face!
                    insert_edge = current2
                    found = true
                }
                current2 = current2.next
            } while (!current2.equals(start_edge2) && !found) // FIXME end early with !found
            assert(insert_edge != null) { " insert_edge != null" }
            assert(found) { " found " }
            g.add_vertex_in_edge(seg_start, insert_edge!!) // insert endpoint in null-edge

            // "process" the adjacent null-edges
            val (next_edge, prev_edge) = g.find_next_prev(start_null_face, seg_start)
            val neg_null_edge: Pair<Vertex?, Face?> = process_null_edge(dir, next_edge, k3_sign, true)
            neg_sep_start = neg_null_edge.first
            face_to_null = neg_null_edge.second
            val pos_null_edge: Pair<Vertex?, Face?> = process_null_edge(dir, prev_edge, k3_sign, false)
            pos_sep_start = pos_null_edge.first
            face_to_null = pos_null_edge.second
            FiveTuple1(
                seg_start,
                start_null_face,
                pos_sep_start,
                neg_sep_start,
                face_to_null
            )
        } else { // no existing null-face
            // create a new null face at start. the face has three vertices/edges:
            //
            // neg_sep -> seg_endp -> pos_sep
            //
            start_null_face = g.add_face() // this face to the left of start->end edge
            start_null_face.is_null_face = true
            seg_start = g.add_vertex(Vertex(start.position, VertexStatus.OUT, VertexType.ENDPOINT))
            seg_start.zero_dist()
            seg_start.set_alfa(dir)
            seg_start.k3 = 0.0
            pos_sep_start = g.add_vertex(Vertex(start.position, VertexStatus.UNDECIDED, VertexType.SEPPOINT))
            neg_sep_start = g.add_vertex(Vertex(start.position, VertexStatus.UNDECIDED, VertexType.SEPPOINT))
            pos_sep_start.zero_dist()
            neg_sep_start.zero_dist()
            if (k3_sign) {
                pos_sep_start.k3 = +1.0
                neg_sep_start.k3 = -1.0
            } else {
                pos_sep_start.k3 = -1.0
                neg_sep_start.k3 = +1.0
            }
            pos_sep_start.set_alfa(dir.xy_perp().mult(+1.0))
            neg_sep_start.set_alfa(dir.xy_perp().mult(-1.0))
            // null-edges around the face
            val (e1, e1_tw) = g.add_twin_edges(seg_start, pos_sep_start)
            val (e2, e2_tw) = g.add_twin_edges(pos_sep_start, neg_sep_start)
            val (e3, e3_tw) = g.add_twin_edges(neg_sep_start, seg_start)

            // e1 -> e2 -> e3 on start_null_face, k=1
            // e1t <- e2t <- e3t on g[start].face, k=1
            g.set_next_cycle(listOf(e1, e2, e3), start_null_face, 1.0)
            val start_face = start.face!!
            val start_face_edge =
                start_face.edge // crazy workaround, because set_next_cycles sets g[face].edge wrong
            // here!
            g.set_next_cycle(listOf(e3_tw, e2_tw, e1_tw), start.face!!, 1.0)
            start_null_face.edge = e1
            start_face.edge = start_face_edge
            e1.type = EdgeType.NULLEDGE
            e2.type = EdgeType.NULLEDGE
            e3.type = EdgeType.NULLEDGE
            e1_tw.type = EdgeType.NULLEDGE
            e3_tw.type = EdgeType.NULLEDGE
            e2_tw.type = EdgeType.NULLEDGE
            start.null_face = start_null_face
            start_null_face.site = start_face.site
            FiveTuple1(
                seg_start,
                start_null_face,
                pos_sep_start,
                neg_sep_start,
                face_to_null
            )
        }
    }

    /**
     * Finds the target of a new SEPARATOR edge
     *
     * @param f    the HEFace on which we search for the target vertex
     * @param endp the end-point of the null-face with the SEPARATOR source
     * @return
     */
    private fun find_separator_target(f: Face, endp: Vertex?): SeparatorTarget {
        // we want to insert a SEPARATOR edge (endp, target) , on the give face f.
        // find and return the target vertex to which the new SEPARATOR edge should
        // connect
        // also return the adjacent next/prev edges
        ///
        // flag==true when an OUT-NEW-IN vertex was found
        ///
        // flag==false when an IN-NEW-OUT vertex was found
        if (endp == null) {
            return SeparatorTarget(null, null, null, false)
        }
        var current_edge = f.edge // start on some edge of the face
        val start_edge = current_edge
        var found = false
        var v_target: Vertex? = null
        var v_previous: Edge? = null
        var v_next: Edge? = null
        var flag = true
        do {
            val next_edge = current_edge.next
            val previous_vertex = current_edge.source
            val current_vertex = current_edge.target
            val next_vertex = next_edge.target
            val out_new_in =
                (previous_vertex.status === VertexStatus.OUT || previous_vertex.status === VertexStatus.UNDECIDED) && current_vertex.status === VertexStatus.NEW && next_vertex.status === VertexStatus.IN
            val in_new_out =
                previous_vertex.status === VertexStatus.IN && current_vertex.status === VertexStatus.NEW && (next_vertex.status === VertexStatus.OUT || next_vertex.status === VertexStatus.UNDECIDED)
            if (out_new_in || in_new_out) {
                if (endp.k3 === current_vertex.k3 && !endp.equals(current_vertex)) {
                    v_target = current_vertex
                    v_previous = current_edge
                    v_next = next_edge
                    flag = if (out_new_in) true else false
                    found = true
                }
            }
            current_edge = current_edge.next
        } while (!current_edge.equals(start_edge) && !found)
        assert(found) { "found" }
        return SeparatorTarget(v_previous, v_target, v_next, flag)
    }

    // prepare null-face
    // next_edge lies on an existing null face
    // - Here we either insert a NEW NORMAL or SEPPOINT in the edge
    // - OR we push and convert an existing vertex.
    // the next_prev flag indicates if we are dealing with the next edge from the
    // new segment-endpoint next_prev=true
    // or if we are dealing with the previous edge (next_prev=false)
    private fun process_null_edge(dir: Point, next_edge: Edge, k3: Boolean, next_prev: Boolean): Pair<Vertex?, Face?> {
        assert(next_edge.type === EdgeType.NULLEDGE) { " next_edge.type == EdgeType.NULLEDGE " }
        val trg = next_edge.target
        val src = next_edge.source
        val adj = if (next_prev) trg else src // this is the vertex adjacent to the end-point, on the null face
        assert((if (next_prev) src else trg).type === VertexType.ENDPOINT) { " (next_prev ? src : trg).type == VertexType.ENDPOINT " }
        var sep_point: Vertex? = null
        val dir_mult = (if (next_prev) +1 else -1).toDouble()
        val sep_dir = dir.xy_perp().mult(dir_mult)
        val sep_alfa = Numeric.diangle(sep_dir.x, sep_dir.y) // alfa of potential SEPPOINT
        val new_k3: Double // k3-value of a new or pushed vertex
        new_k3 = if (next_prev) {
            (if (k3) +1 else -1).toDouble()
        } else {
            (if (k3) -1 else +1).toDouble()
        }
        return if (adj.type === VertexType.ENDPOINT) { // target is endpoint
            // insert a normal vertex, positioned at mid-alfa between src/trg.
            val new_v = g.add_vertex(Vertex(src.position, VertexStatus.NEW, VertexType.NORMAL, src.position))
            val mid = Numeric.diangle_mid(src.alfa, trg.alfa)
            new_v.alfa = mid
            modified_vertices.add(new_v)
            g.add_vertex_in_edge(new_v, next_edge)
            new_v.k3 = new_k3
            Pair<Vertex?, Face?>(null, null)
        } else {
            // Not an ENDPOINT vertex.
            assert(adj.type !== VertexType.ENDPOINT) { " adj.type != VertexType.ENDPOINT " }
            var mid = 0.0
            var seppoint_pred = false
            var parallel_pred = false

            // set the two predicates
            if (next_prev) {
                val next_next = next_edge.next
                val next_previous = g.previous_edge(next_edge)
                val next_trg = next_next.target // source
                mid = Numeric.diangle_mid(src.alfa, next_trg.alfa) // prev_src, trg
                seppoint_pred = next_trg.type !== VertexType.ENDPOINT
                val next_out_trg: Vertex? = null_vertex_target(next_edge.target)
                val prev_out_trg: Vertex? = null_vertex_target(next_previous.source)
                if (next_out_trg != null && prev_out_trg != null) {
                    parallel_pred =
                        ((next_out_trg.status === VertexStatus.OUT || next_out_trg.status === VertexStatus.NEW || next_out_trg.status === VertexStatus.UNDECIDED)
                                && (prev_out_trg.status === VertexStatus.OUT || prev_out_trg.status === VertexStatus.NEW || prev_out_trg.status === VertexStatus.UNDECIDED))
                }
            } else { // !next_prev
                val prev_prev = g.previous_edge(next_edge)
                val next_next2 = next_edge.next
                val prev_src = prev_prev.source
                mid = Numeric.diangle_mid(prev_src.alfa, trg.alfa)
                seppoint_pred = prev_src.type !== VertexType.ENDPOINT
                val next_out_trg2: Vertex? = null_vertex_target(next_edge.source)
                val prev_out_trg2: Vertex? = null_vertex_target(next_next2.target)
                if (next_out_trg2 != null && prev_out_trg2 != null) {
                    parallel_pred =
                        ((next_out_trg2.status === VertexStatus.OUT || next_out_trg2.status === VertexStatus.NEW || next_out_trg2.status === VertexStatus.UNDECIDED)
                                && (prev_out_trg2.status === VertexStatus.OUT || prev_out_trg2.status === VertexStatus.NEW || prev_out_trg2.status === VertexStatus.UNDECIDED))
                }
            } // predicates now set.
            val adj_out: Vertex? = null_vertex_target(adj)
            assert(adj_out != null) { "adj_out != null" }
            if (adj_out!!.status === VertexStatus.OUT || adj_out!!.status === VertexStatus.UNDECIDED) {
                sep_point = add_separator_vertex(src, next_edge, sep_dir)
                sep_point.k3 = new_k3
                Pair<Vertex, Face?>(sep_point, null)
            } else {
                // target is not endpoint so we push and convert the vertex
                if (seppoint_pred) {
                    // the pushed vertex becomes a SEPPOINT
                    adj.alfa = sep_alfa
                    adj.type = VertexType.SEPPOINT
                    adj.status = VertexStatus.NEW
                    sep_point = adj
                } else {
                    // otherwise it becomes a normal NEW vertex
                    adj.alfa = mid
                    adj.type = VertexType.NORMAL
                    adj.status = VertexStatus.NEW
                }
                adj.k3 = new_k3
                modified_vertices.add(adj)
                Pair<Vertex?, Face?>(sep_point, null)
            }
        }
    }

    // add a SEPPOINT vertex into the given null-edge
    ///
    // \param endp the endpoint corresponding to the null-edge/null-face
    // \param edge the null-edge into which we insert the new vertex
    // \param sep_dir direction for setting alfa of the new vertex
    private fun add_separator_vertex(endp: Vertex, edge: Edge, sep_dir: Point): Vertex {
        val sep = g.add_vertex(Vertex(endp.position, VertexStatus.OUT, VertexType.SEPPOINT))
        sep.set_alfa(sep_dir)
        g.add_vertex_in_edge(sep, edge)
        modified_vertices.add(sep)
        return sep
    }

    /**
     * one-argument version of repair_face() used by insert_point_site()
     *
     * @see .repair_face
     * @param f
     */
    private fun repair_face(f: Face) {
        repair_face(
            f, Pair<Vertex?, Vertex?>(null, null), Pair<Face?, Face?>(null, null),
            Pair<Face?, Face?>(null, null)
        )
    }

    /**
     * Repairs next-pointers of HEFace
     *
     * @param f
     * @param segment
     * @param nulled_faces
     * @param null_face
     */
    private fun repair_face(
        f: Face, segment: Pair<Vertex?, Vertex?>, nulled_faces: Pair<Face?, Face?>,
        null_face: Pair<Face?, Face?>
    ) {
        // start on g[newface].edge, walk around the face and repair the next-pointers
        // this is called on the newly created face after all NEW-NEW edges have been
        // added
        var current_edge = f.edge
        val start_edge = current_edge
        var c = 0
        do {
            assert(vd_checker.check_edge(current_edge)) { " vd_checker.check_edge(current_edge) " }
            val current_target = current_edge.target // an edge on the new face
            val current_source = current_edge.source
            for (e in current_target.out_edges) {
                val out_target = e.target
                if (!out_target.equals(current_source) && (out_target.status === VertexStatus.NEW || out_target.type === VertexType.ENDPOINT || out_target.type === VertexType.SEPPOINT)) { // these
                    // are
                    // the
                    // possible
                    // vertices
                    // we
                    // want
                    // to
                    // go
                    // to

                    // special cases where we do a brute-force face-assignment for a null-edge, or a
                    // separator
                    if ((e.type === EdgeType.NULLEDGE && current_edge.type !== EdgeType.NULLEDGE &&  // only one
                                // null-edge in
                                // succession!
                                ( // from sep to end
                                        current_target.type === VertexType.SEPPOINT && out_target.type === VertexType.ENDPOINT || current_source.type === VertexType.ENDPOINT && current_target.type === VertexType.ENDPOINT
                                                || out_target.equals(segment.first)
                                                || out_target.equals(segment.second))
                                && !e.face.equals(null_face.first) &&  // not along a null-face edge!
                                !e.face.equals(null_face.second)) || e.face.equals(nulled_faces.first) // edge
                        // previously
                        // belonged
                        // to
                        // point-site
                        // that
                        // has
                        // disappeared
                        || e.face.equals(nulled_faces.second)
                    ) {
                        e.face = f // override face-assignment!
                        e.k = current_edge.k // override k-assignment!
                    }

                    // the next vertex should not where we came from
                    // and it should be on the same face.
                    if (e.face === f) {
                        current_edge.next = e // this is the edge we want to take
                        assert(current_edge.k == e.k) { " current_edge.k == e.k " }
                        assert(vd_checker.current_face_equals_next_face(current_edge)) { " vd_checker.current_face_equals_next_face( current_edge ) " }
                    }
                }
            }
            current_edge = current_edge.next // jump to the next edge
            c++
            if (c > 30000) {
                throw AssertionError("c < 30000")
            }
        } while (!current_edge.equals(start_edge))
    }

    /**
     * Removes the IN vertices of the delete-tree; removes the IN vertices stored in
     * v0 (and associated IN-NEW edges)
     */
    private fun remove_vertex_set() {
        for (v in v0) {
            // it should now be safe to delete all IN vertices
            assert(v.status === VertexStatus.IN) { " v.status == VertexStatus.IN " }
            g.delete_vertex(v) // this also removes edges connecting to v
            modified_vertices.remove(v)
        }
    }

    // remove all SPLIT type vertices on the HEFace \a f
    private fun remove_split_vertex(f: Face) {
        assert(vd_checker.face_ok(f)) { " vd_checker.face_ok( f ) " }
        var v: Vertex?
        while (find_split_vertex(f).also { v = it } != null) {
            assert(v!!.type === VertexType.SPLIT) { "v.type == VertexType.SPLIT" }
            g.remove_deg2_vertex(v!!)
            modified_vertices.remove(v)
            assert(vd_checker.face_ok(f)) { " vd_checker.face_ok( f ) " }
        }
        assert(vd_checker.face_ok(f)) { " vd_checker.face_ok( f ) " }
    }

    // reset status of modified_vertices and incident_faces
    ///
    // at the end after an incremental insertion of a new site,
    // reset status of modified_vertices to UNDECIDED and incident_faces to
    // NONINCIDENT,
    // so that we are ready for the next insertion.
    private fun reset_status() {
        for (v in modified_vertices) {
            v.reset_status()
        }
        modified_vertices.clear()
        for (f in incident_faces) {
            f.status = FaceStatus.NONINCIDENT
        }
        incident_faces.clear()
        v0.clear()
    }

    // count number of NEW vertices on the given face \a f
    private fun num_new_vertices(f: Face): Int {
        var current = f.edge
        val start = current
        var count = 0
        do {
            val v = current.target
            if (v.status === VertexStatus.NEW && v.type !== VertexType.SEPPOINT) {
                count++
            }
            current = current.next
        } while (!current.equals(start))
        return count
    }

    companion object {
        // reset vertex index count \todo not very elegant...
        fun reset_vertex_count() {
            Vertex.reset_count()
        }
    }
}
