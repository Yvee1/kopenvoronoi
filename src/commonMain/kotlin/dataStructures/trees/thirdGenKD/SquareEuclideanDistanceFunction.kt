package dataStructures.trees.thirdGenKD

/**
 *
 */
class SquareEuclideanDistanceFunction : DistanceFunction {
    override fun distance(p1: DoubleArray, p2: DoubleArray): Double {
        var d = 0.0
        for (i in p1.indices) {
            val diff = p1[i] - p2[i]
            d += diff * diff
        }
        return d
    }

    override fun distanceToRect(point: DoubleArray, min: DoubleArray, max: DoubleArray): Double {
        var d = 0.0
        for (i in point.indices) {
            var diff = 0.0
            if (point[i] > max[i]) {
                diff = point[i] - max[i]
            } else if (point[i] < min[i]) {
                diff = point[i] - min[i]
            }
            d += diff * diff
        }
        return d
    }
}