package kopenvoronoi.geometry

enum class FaceStatus {
    /**
     * INCIDENT faces contain one or more IN-vertex
     */
    INCIDENT,

    /**
     * NONINCIDENT faces contain only OUT/UNDECIDED-vertices
     */
    NONINCIDENT
}
