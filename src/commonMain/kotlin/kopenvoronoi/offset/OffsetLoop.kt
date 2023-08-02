package kopenvoronoi.offset

//a single offset loop
class OffsetLoop {
    var vertices: MutableList<OffsetVertex> = mutableListOf() // < list of offsetvertices in this loop
    var offset_distance = 0.0
    fun add(v: OffsetVertex) {
        vertices.add(v)
    }
}
