package kopenvoronoi.offset

import kopenvoronoi.geometry.Point

//\brief offset-element of LineSite
class LineOfs(p1: Point, p2: Point) : Ofs() {
    protected var _start // < start point
            : Point
    protected var _end // < end point
            : Point

    // \param p1 start point
    // \param p2 end point
    init {
        _start = p1
        _end = p2
    }

    override fun radius(): Double {
        return (-1).toDouble()
    }

    override fun center(): Point {
        return Point(0.0, 0.0)
    }

    override fun start(): Point {
        return _start
    }

    override fun end(): Point {
        return _end
    }
}
