package kopenvoronoi.vertex

import kopenvoronoi.geometry.Edge
import kopenvoronoi.geometry.Face
import kopenvoronoi.geometry.Point
import kopenvoronoi.util.Numeric

/**
 * A vertex in the voronoi diagram an object of this type is held in the
 * BGL-graph for each vertex.
 */
class Vertex {
    var out_edges: MutableList<Edge> = mutableListOf()
    var in_edges: MutableList<Edge> = mutableListOf()

    /**
     * vertex status. updated/changed during an incremental graph update
     */
    var status: VertexStatus? = null

    /**
     * The type of the vertex. Never(?) changes
     */
    var type: VertexType? = null
    var max_error = 0.0 // < \todo what is this? remove?
    var in_queue = false // < flag for indicating wether vertex is in the vertexQueue
    lateinit var position: Point  // < the position of the vertex.
    var k3 = 0.0 // < the offset-direction {-1,+1} of this vertex to the newly inserted site.
    var alfa = 0.0 // < diangle for a null-vertex. only for debug-drawing
    var null_face: Face? = null // < if this is a null-face, a handle to the null-face
    var face: Face? = null // < the face of this vertex, if the vertex is a point-site
    var r = 0.0 // < clearance-disk radius, i.e. the closest Site is at this distance

    constructor()

    fun degree(): Int {
        return out_edges.size + in_edges.size
    }

    // ctor with given status and type
    constructor(p: Point, st: VertexStatus, t: VertexType) {
        init(p, st, t)
    }

    // ctor with initial apex Point
    constructor(p: Point, st: VertexStatus, t: VertexType, initDist: Point) {
        init(p, st, t, initDist)
    }

    // ctor with initial k3-value
    constructor(p: Point, st: VertexStatus, t: VertexType, initDist: Point, lk3: Double) {
        init(p, st, t, initDist, lk3)
    }

    // ctor with initial clearance-disk radius
    constructor(p: Point, st: VertexStatus, t: VertexType, init_radius: Double) {
        init(p, st, t)
        r = init_radius
    }

    // set index, increase count, initialize in_queue to false.
    fun init() {
        count++
        in_queue = false
        alfa = -1.0 // invalid/non-initialized alfa value
        null_face = null
        type = VertexType.NORMAL
        face = null
        max_error = 0.0
    }

    // set position and status
    fun init(p: Point, st: VertexStatus) {
        init()
        position = p
        status = st
    }

    // set position, status and type
    fun init(p: Point, st: VertexStatus, t: VertexType) {
        init(p, st)
        type = t
    }

    // set position, status, type, and clearance-disk through given apex-point
    fun init(p: Point, st: VertexStatus, t: VertexType, initDist: Point) {
        init(p, st, t)
        init_dist(initDist)
    }

    // set position, status, type, clerance-disk radius, and k3-side
    fun init(p: Point, st: VertexStatus, t: VertexType, initDist: Point, lk3: Double) {
        init(p, st, t, initDist)
        k3 = lk3
    }

    // set in_queue false, and status to ::UNDECIDED
    fun reset_status() {
        in_queue = false
        status = VertexStatus.UNDECIDED
    }

    fun set_alfa(dir: Point) {
        alfa = Numeric.diangle(dir.x, dir.y)
    }

    // initialize clerance-disk
    fun init_dist(p: Point) {
        r = dist(p)
    }

    // return distance to a point from this vertex
    fun dist(p: Point): Double {
        return position.sub(p).norm()
    }

    // set clearance-disk to zero
    fun zero_dist() {
        r = 0.0
    }

    // return clearance disk-radius
    fun dist(): Double {
        return r
    }

    // in-circle predicate
    fun in_circle(p: Point): Double {
        return dist(p) - r
    }

    override fun toString(): String {
        return "V(${position})"
    }

    companion object {
        var count = 0 // < global vertex count TODO hold this in hedigraph instead?

        // A map of this type is used by VoronoiDiagramChecker to check that all
        // vertices
        // have the expected (correct) degree (i.e. number of edges)
        // map for checking topology correctness
        var expected_degree: MutableMap<VertexType, Int> = mutableMapOf()

        init {
            expected_degree.put(VertexType.OUTER, 4) // special outer vertices
            expected_degree.put(VertexType.NORMAL, 6) // normal vertex in the graph
            expected_degree.put(VertexType.POINTSITE, 0) // point site
            expected_degree.put(VertexType.ENDPOINT, 6) // end-point of line or arc
            expected_degree.put(VertexType.SEPPOINT, 6) // end-point of separator
            expected_degree.put(VertexType.SPLIT, 4) // split point, to avoid loops in delete-tree
            expected_degree.put(VertexType.APEX, 4) // apex point on quadratic bisector
        }

        // reset the index count
        fun reset_count() {
            count = 0
        }
    }
}
