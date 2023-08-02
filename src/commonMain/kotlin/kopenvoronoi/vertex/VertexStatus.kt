package kopenvoronoi.vertex

/**
 * As we incrementally construct the diagram voronoi-vertices can have one of
 * these four different states.
 */
enum class VertexStatus {
    /**
     * OUT-vertices will not be deleted
     */
    OUT,

    /**
     * IN-vertices will be deleted
     */
    IN,

    /**
     * UNDECIDED-vertices have not been examied yet
     */
    UNDECIDED,

    /**
     * NEW-vertices are constructed on OUT-IN edges
     */
    NEW
}
