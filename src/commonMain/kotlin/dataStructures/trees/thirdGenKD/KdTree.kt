package dataStructures.trees.thirdGenKD

import dataStructures.BinaryHeap
import dataStructures.MaxHeap
import dataStructures.MinHeap
import kotlin.jvm.JvmOverloads
import kotlin.math.min

class KdTree<T> @JvmOverloads constructor(dimensions: Int, bucketCapacity: Int = 24) :
    KdNode<T>(dimensions, bucketCapacity) {
    fun getNearestNeighborIterator(
        searchPoint: DoubleArray, maxPointsReturned: Int,
        distanceFunction: DistanceFunction
    ): NearestNeighborIterator<T> {
        return NearestNeighborIterator<T>(this, searchPoint, maxPointsReturned, distanceFunction)
    }

    fun findNearestNeighbors(
        searchPoint: DoubleArray, maxPointsReturned: Int,
        distanceFunction: DistanceFunction
    ): MaxHeap<T> {
        val pendingPaths = BinaryHeap.Min<KdNode<T>>()
        val evaluatedPoints = BinaryHeap.Max<T>()
        val pointsRemaining: Int = min(maxPointsReturned, size())
        pendingPaths.offer(0.0, this)
        while (pendingPaths.size() > 0 && (evaluatedPoints.size() < pointsRemaining || pendingPaths.minKey < evaluatedPoints.maxKey)) {
            nearestNeighborSearchStep(pendingPaths, evaluatedPoints, pointsRemaining, distanceFunction, searchPoint)
        }
        return evaluatedPoints
    }

    companion object {
        fun <T> nearestNeighborSearchStep(
            pendingPaths: MinHeap<KdNode<T>>, evaluatedPoints: MaxHeap<T>,
            desiredPoints: Int, distanceFunction: DistanceFunction, searchPoint: DoubleArray
        ) {
            // If there are pending paths possibly closer than the nearest evaluated point,
            // check it out
            var cursor = pendingPaths.min
            pendingPaths.removeMin()

            // Descend the tree, recording paths not taken
            while (!cursor.isLeaf) {
                var pathNotTaken: KdNode<T>
                if (searchPoint[cursor.splitDimension] > cursor.splitValue) {
                    pathNotTaken = cursor.left!!
                    cursor = cursor.right!!
                } else {
                    pathNotTaken = cursor.right!!
                    cursor = cursor.left!!
                }
                val otherDistance = distanceFunction.distanceToRect(
                    searchPoint, pathNotTaken.minBound!!,
                    pathNotTaken.maxBound!!
                )
                // Only add a path if we either need more points or it's closer than furthest
                // point on list so far
                if (evaluatedPoints.size() < desiredPoints || otherDistance <= evaluatedPoints.maxKey) {
                    pendingPaths.offer(otherDistance, pathNotTaken)
                }
            }
            if (cursor.singlePoint) {
                val nodeDistance = distanceFunction.distance(cursor.points!![0]!!, searchPoint)
                // Only add a point if either need more points or it's closer than furthest on
                // list so far
                if (evaluatedPoints.size() < desiredPoints || nodeDistance <= evaluatedPoints.maxKey) {
                    for (i in 0 until cursor.size()) {
                        val value = cursor.data!!.get(i) as T

                        // If we don't need any more, replace max
                        if (evaluatedPoints.size() == desiredPoints) {
                            evaluatedPoints.replaceMax(nodeDistance, value)
                        } else {
                            evaluatedPoints.offer(nodeDistance, value)
                        }
                    }
                }
            } else {
                // Add the points at the cursor
                for (i in 0 until cursor.size()) {
                    val point = cursor.points!![i]!!
                    val value = cursor.data!!.get(i) as T
                    val distance = distanceFunction.distance(point, searchPoint)
                    // Only add a point if either need more points or it's closer than furthest on
                    // list so far
                    if (evaluatedPoints.size() < desiredPoints) {
                        evaluatedPoints.offer(distance, value)
                    } else if (distance < evaluatedPoints.maxKey) {
                        evaluatedPoints.replaceMax(distance, value)
                    }
                }
            }
        }
    }
}
