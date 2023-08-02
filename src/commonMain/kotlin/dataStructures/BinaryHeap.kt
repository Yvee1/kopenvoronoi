package dataStructures

/**
 * An implementation of an implicit binary heap. Min-heap and max-heap both
 * supported
 */
abstract class BinaryHeap<T> protected constructor(capacity: Int, private val direction: Int) {
    private var data: Array<Any?>
    private var keys: DoubleArray
    private var capacity: Int
    private var size: Int

    init {
        data = arrayOfNulls(capacity)
        keys = DoubleArray(capacity)
        this.capacity = capacity
        size = 0
    }

    fun offer(key: Double, value: T) {
        // If move room is needed, double array size
        if (size >= capacity) {
            capacity *= 2
            data = data.copyOf(capacity)
            keys = keys.copyOf(capacity)
        }

        // Insert new value at the end
        data[size] = value
        keys[size] = key
        siftUp(size)
        size++
    }

    protected fun removeTip() {
        if (size == 0) {
            throw IllegalStateException()
        }
        size--
        data[0] = data[size]
        keys[0] = keys[size]
        data[size] = null
        siftDown(0)
    }

    protected fun replaceTip(key: Double, value: T) {
        if (size == 0) {
            throw IllegalStateException()
        }
        data[0] = value
        keys[0] = key
        siftDown(0)
    }

    protected val tip: T
        get() {
            if (size == 0) {
                throw IllegalStateException()
            }
            return data[0] as T
        }
    protected val tipKey: Double
        get() {
            if (size == 0) {
                throw IllegalStateException()
            }
            return keys[0]
        }

    private fun siftUp(c: Int) {
        var c = c
        var p = (c - 1) / 2
        while (c != 0 && direction * keys[c] > direction * keys[p]) {
            val pData = data[p]
            val pDist = keys[p]
            data[p] = data[c]
            keys[p] = keys[c]
            data[c] = pData
            keys[c] = pDist
            c = p
            p = (c - 1) / 2
        }
    }

    private fun siftDown(p: Int) {
        var p = p
        var c = p * 2 + 1
        while (c < size) {
            if (c + 1 < size && direction * keys[c] < direction * keys[c + 1]) {
                c++
            }
            if (direction * keys[p] < direction * keys[c]) {
                // Swap the points
                val pData = data[p]
                val pDist = keys[p]
                data[p] = data[c]
                keys[p] = keys[c]
                data[c] = pData
                keys[c] = pDist
            } else {
                break
            }
            p = c
            c = p * 2 + 1
        }
    }

    fun size(): Int {
        return size
    }

    fun capacity(): Int {
        return capacity
    }

    class Max<T> : BinaryHeap<T>, MaxHeap<T> {
        constructor() : super(defaultCapacity, 1)
        constructor(capacity: Int) : super(capacity, 1)

        override fun removeMax() {
            removeTip()
        }

        override fun replaceMax(key: Double, value: T) {
            replaceTip(key, value)
        }

        override val max: T
            get() = tip
        override val maxKey: Double
            get() = tipKey
    }

    class Min<T> : BinaryHeap<T>, MinHeap<T> {
        constructor() : super(defaultCapacity, -1)
        constructor(capacity: Int) : super(capacity, -1)

        override fun removeMin() {
            removeTip()
        }

        override fun replaceMin(key: Double, value: T) {
            replaceTip(key, value)
        }

        override val min: T
            get() = tip
        override val minKey: Double
            get() = tipKey
    }

    companion object {
        protected const val defaultCapacity = 64
    }
}