package dataStructures.trees.thirdGenKD

import dataStructures.BinaryHeap
import dataStructures.IntervalHeap
import dataStructures.MinHeap
import kotlin.math.min

/**
 *
 */
class NearestNeighborIterator<T>(
    treeRoot: KdNode<T>, searchPoint: DoubleArray, maxPointsReturned: Int,
    distanceFunction: DistanceFunction
) : MutableIterator<T>, Iterable<T> {
    private val distanceFunction: DistanceFunction
    private val searchPoint: DoubleArray
    private val pendingPaths: MinHeap<KdNode<T>>
    private val evaluatedPoints: IntervalHeap<T>
    private var pointsRemaining: Int
    private var lastDistanceReturned = 0.0

    init {
        this.searchPoint = searchPoint.copyOf(searchPoint.size)
        pointsRemaining = min(maxPointsReturned, treeRoot.size())
        this.distanceFunction = distanceFunction
        pendingPaths = BinaryHeap.Min()
        pendingPaths.offer(0.0, treeRoot)
        evaluatedPoints = IntervalHeap()
    }

    /* -------- INTERFACE IMPLEMENTATION -------- */
    override fun hasNext(): Boolean {
        return pointsRemaining > 0
    }

    override fun next(): T {
        if (!hasNext()) {
            throw IllegalStateException("NearestNeighborIterator has reached end!")
        }
        while (pendingPaths.size() > 0
            && (evaluatedPoints.size() == 0 || pendingPaths.minKey < evaluatedPoints.minKey)
        ) {
            KdTree.nearestNeighborSearchStep(
                pendingPaths, evaluatedPoints, pointsRemaining, distanceFunction,
                searchPoint
            )
        }

        // Return the smallest distance point
        pointsRemaining--
        lastDistanceReturned = evaluatedPoints.minKey
        val value = evaluatedPoints.min
        evaluatedPoints.removeMin()
        return value
    }

    fun distance(): Double {
        return lastDistanceReturned
    }

    override fun remove() {
        throw UnsupportedOperationException()
    }

    override fun iterator(): Iterator<T> {
        return this
    }
}
