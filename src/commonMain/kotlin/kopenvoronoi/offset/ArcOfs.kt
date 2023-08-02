package kopenvoronoi.offset

import kopenvoronoi.geometry.Point

//\brief offset-element of PointSite or ArcSite
class ArcOfs(p1: Point, p2: Point, cen: Point, rad: Double) : Ofs() {
    var _start // < start
            : Point
    var _end // < end
            : Point
    var c // < center
            : Point
    var r // < radius
            : Double

    // \param p1 start Point
    // \param p2 end Point
    // \param cen center Point
    // \param rad radius
    init {
        _start = p1
        _end = p2
        c = cen
        r = rad
    }

    override fun toString(): String {
        TODO()
//        return String.format("ArcOfs from %s to %s r=%f\n", _start, _end, r)
    }

    override fun radius(): Double {
        return r
    }

    override fun center(): Point {
        return c
    }

    override fun start(): Point {
        return _start
    }

    override fun end(): Point {
        return _end
    }
}
