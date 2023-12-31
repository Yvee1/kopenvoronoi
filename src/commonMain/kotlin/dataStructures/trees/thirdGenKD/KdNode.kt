package dataStructures.trees.thirdGenKD

open class KdNode<T> protected constructor(// All types
    protected var dimensions: Int, protected var bucketCapacity: Int
) {
    private var size = 0

    // Leaf only
    var points: Array<DoubleArray?>?
    var data: Array<Any?>?

    // Stem only
    var left: KdNode<T>? = null
    var right: KdNode<T>? = null
    var splitDimension = 0
    var splitValue = 0.0

    // Bounds
    var minBound: DoubleArray? = null
    var maxBound: DoubleArray? = null
    var singlePoint = true

    init {
        // Init base

        // Init leaf elements
        points = arrayOfNulls(bucketCapacity + 1)
        data = arrayOfNulls(bucketCapacity + 1)
    }

    /* -------- SIMPLE GETTERS -------- */
    fun size(): Int {
        return size
    }

    val isLeaf: Boolean
        get() = points != null

    /* -------- OPERATIONS -------- */
    fun addPoint(point: DoubleArray, value: T) {
        var cursor: KdNode<T>? = this
        while (!cursor!!.isLeaf) {
            cursor.extendBounds(point)
            cursor.size++
            if (point[cursor.splitDimension] > cursor.splitValue) {
                cursor = cursor.right
            } else {
                cursor = cursor.left
            }
        }
        cursor.addLeafPoint(point, value)
    }

    /* -------- INTERNAL OPERATIONS -------- */
    fun addLeafPoint(point: DoubleArray, value: T) {
        // Add the data point
        points!![size] = point
        data!![size] = value
        extendBounds(point)
        size++
        if (size == points!!.size - 1) {
            // If the node is getting too large
            if (calculateSplit()) {
                // If the node successfully had it's split value calculated, split node
                splitLeafNode()
            } else {
                // If the node could not be split, enlarge node
                increaseLeafCapacity()
            }
        }
    }

    private fun checkBounds(point: DoubleArray): Boolean {
        for (i in 0 until dimensions) {
            if (point[i] > maxBound!![i]) {
                return false
            }
            if (point[i] < minBound!![i]) {
                return false
            }
        }
        return true
    }

    private fun extendBounds(point: DoubleArray) {
        if (minBound == null) {
            minBound = point.copyOf(dimensions)
            maxBound = point.copyOf(dimensions)
            return
        }
        for (i in 0 until dimensions) {
            if (point[i].isNaN()) {
                if (!minBound!![i].isNaN() || !maxBound!![i].isNaN()) {
                    singlePoint = false
                }
                minBound!![i] = Double.NaN
                maxBound!![i] = Double.NaN
            } else if (minBound!![i] > point[i]) {
                minBound!![i] = point[i]
                singlePoint = false
            } else if (maxBound!![i] < point[i]) {
                maxBound!![i] = point[i]
                singlePoint = false
            }
        }
    }

    private fun increaseLeafCapacity() {
        points = points!!.copyOf(points!!.size * 2)
        data = data!!.copyOf(data!!.size * 2)
    }

    private fun calculateSplit(): Boolean {
        if (singlePoint) {
            return false
        }
        var width = 0.0
        for (i in 0 until dimensions) {
            var dwidth = maxBound!![i] - minBound!![i]
            if (dwidth.isNaN()) {
                dwidth = 0.0
            }
            if (dwidth > width) {
                splitDimension = i
                width = dwidth
            }
        }
        if (width == 0.0) {
            return false
        }

        // Start the split in the middle of the variance
        splitValue = (minBound!![splitDimension] + maxBound!![splitDimension]) * 0.5

        // Never split on infinity or NaN
        if (splitValue == Double.POSITIVE_INFINITY) {
            splitValue = Double.MAX_VALUE
        } else if (splitValue == Double.NEGATIVE_INFINITY) {
            splitValue = -Double.MAX_VALUE
        }

        // Don't let the split value be the same as the upper value as
        // can happen due to rounding errors!
        if (splitValue == maxBound!![splitDimension]) {
            splitValue = minBound!![splitDimension]
        }

        // Success
        return true
    }

    // I have no idea what the following code does, it's not written by me.
    // Thus, I can only suppress the warnings generated by `(T) oldData` cast.
    private fun splitLeafNode() {
        right = KdNode(dimensions, bucketCapacity)
        left = KdNode(dimensions, bucketCapacity)

        // Move locations into children
        for (i in 0 until size) {
            val oldLocation = points!![i]
            val oldData = data!![i]
            if (oldLocation!![splitDimension] > splitValue) {
                right!!.addLeafPoint(oldLocation, oldData as T)
            } else {
                left!!.addLeafPoint(oldLocation, oldData as T)
            }
        }
        points = null
        data = null
    }
}
