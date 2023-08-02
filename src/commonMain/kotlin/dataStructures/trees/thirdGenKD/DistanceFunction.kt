package dataStructures.trees.thirdGenKD

interface DistanceFunction {
    fun distance(p1: DoubleArray, p2: DoubleArray): Double
    fun distanceToRect(point: DoubleArray, min: DoubleArray, max: DoubleArray): Double
}
