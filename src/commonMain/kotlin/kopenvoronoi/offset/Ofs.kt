package kopenvoronoi.offset

import kopenvoronoi.geometry.Point
//\brief base-class for offset-elements
///
//preliminary offset-prerensentations. experiental...
abstract class Ofs {
    // radius, -1 if line
    abstract fun radius(): Double // {return -1;}

    // center (for arc)
    abstract fun center(): Point // {return Point(0,0);}

    // start point
    abstract fun start(): Point // {return Point(0,0);}

    // end point
    abstract fun end(): Point // {return Point(0,0);}
}
