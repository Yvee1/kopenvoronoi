package kopenvoronoi.geometry

import kopenvoronoi.site.Site

class Face {
    lateinit var edge: Edge
    lateinit var site: Site
    var status: FaceStatus? = null
    var is_null_face = false
    override fun toString(): String {
        val sb = StringBuilder()
        sb.append("F(")
        var current: Edge? = edge
        var c = 0
        do {
            if (current == null) {
                break
            }
            sb.append(current.source.position)
            sb.append(">")
            current = current.next
            c++
        } while (current !== edge && c < 100)
        if (c >= 100) {
            sb.append("...")
        }
        sb.append(")")
        return sb.toString()
    }
}
