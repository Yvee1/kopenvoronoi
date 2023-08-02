package kopenvoronoi.vertex

import kopenvoronoi.geometry.Point

//\brief a new vertex position solution (position, offset-distance, side)
///
//includes the offset-distamce t, and the offset direction k3
class Solution(pt: Point, tv: Double, k: Double) {
    // position
    var p: Point

    // clearance-disk radius
    var t: Double

    // offset direction to third adjacent Site
    var k3: Double

    // \param pt vertex position
    // \param tv clearance-disk radius
    // \param k offset direction
    init {
        p = pt
        t = tv
        k3 = k
    }

    override fun toString(): String {
        return "Solution(p = ${p}, t = ${t}, k = ${k3.toInt()})"
    }
}
