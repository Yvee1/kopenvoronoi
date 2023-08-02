package dataStructures

import kotlin.jvm.JvmOverloads

/**
 * An implementation of an implicit binary interval heap.
 */
class IntervalHeap<T> @JvmOverloads constructor(capacity: Int = defaultCapacity) : MinHeap<T>, MaxHeap<T> {
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

    override fun offer(key: Double, value: T) {
        // If move room is needed, double array size
        if (size >= capacity) {
            capacity *= 2
            data = data.copyOf(capacity)
            keys = keys.copyOf(capacity)
        }

        // Insert new value at the end
        size++
        data[size - 1] = value
        keys[size - 1] = key
        siftInsertedValueUp()
    }

    override fun removeMin() {
        if (size == 0) {
            throw IllegalStateException()
        }
        size--
        data[0] = data[size]
        keys[0] = keys[size]
        data[size] = null
        siftDownMin(0)
    }

    override fun replaceMin(key: Double, value: T) {
        if (size == 0) {
            throw IllegalStateException()
        }
        data[0] = value
        keys[0] = key
        if (size > 1) {
            // Swap with pair if necessary
            if (keys[1] < key) {
                swap(0, 1)
            }
            siftDownMin(0)
        }
    }

    override fun removeMax() {
        if (size == 0) {
            throw IllegalStateException()
        } else if (size == 1) {
            removeMin()
            return
        }
        size--
        data[1] = data[size]
        keys[1] = keys[size]
        data[size] = null
        siftDownMax(1)
    }

    override fun replaceMax(key: Double, value: T) {
        if (size == 0) {
            throw IllegalStateException()
        } else if (size == 1) {
            replaceMin(key, value)
            return
        }
        data[1] = value
        keys[1] = key
        // Swap with pair if necessary
        if (key < keys[0]) {
            swap(0, 1)
        }
        siftDownMax(1)
    }

    override val min: T
        get() {
            if (size == 0) {
                throw IllegalStateException()
            }
            return data[0] as T
        }
    override val max: T
        get() {
            if (size == 0) {
                throw IllegalStateException()
            } else if (size == 1) {
                return data[0] as T
            }
            return data[1] as T
        }
    override val minKey: Double
        get() {
            if (size == 0) {
                throw IllegalStateException()
            }
            return keys[0]
        }
    override val maxKey: Double
        get() {
            if (size == 0) {
                throw IllegalStateException()
            } else if (size == 1) {
                return keys[0]
            }
            return keys[1]
        }

    private fun swap(x: Int, y: Int): Int {
        val yData = data[y]
        val yDist = keys[y]
        data[y] = data[x]
        keys[y] = keys[x]
        data[x] = yData
        keys[x] = yDist
        return y
    }

    /**
     * Min-side (u % 2 == 0): - leftchild: 2u + 2 - rightchild: 2u + 4 - parent:
     * (x/2-1)&~1
     *
     * Max-side (u % 2 == 1): - leftchild: 2u + 1 - rightchild: 2u + 3 - parent:
     * (x/2-1)|1
     */
    private fun siftInsertedValueUp() {
        var u = size - 1
        if (u == 0) {
            // Do nothing if it's the only element!
        } else if (u == 1) {
            // If it is the second element, just sort it with it's pair
            if (keys[u] < keys[u - 1]) { // If less than it's pair
                swap(u, u - 1) // Swap with it's pair
            }
        } else if (u % 2 == 1) {
            // Already paired. Ensure pair is ordered right
            val p = u / 2 - 1 or 1 // The larger value of the parent pair
            if (keys[u] < keys[u - 1]) { // If less than it's pair
                u = swap(u, u - 1) // Swap with it's pair
                if (keys[u] < keys[p - 1]) { // If smaller than smaller parent pair
                    // Swap into min-heap side
                    u = swap(u, p - 1)
                    siftUpMin(u)
                }
            } else if (keys[u] > keys[p]) { // If larger that larger parent pair
                // Swap into max-heap side
                u = swap(u, p)
                siftUpMax(u)
            }
        } else {
            // Inserted in the lower-value slot without a partner
            val p = u / 2 - 1 or 1 // The larger value of the parent pair
            if (keys[u] > keys[p]) { // If larger that larger parent pair
                // Swap into max-heap side
                u = swap(u, p)
                siftUpMax(u)
            } else if (keys[u] < keys[p - 1]) { // If smaller than smaller parent pair
                // Swap into min-heap side
                u = swap(u, p - 1)
                siftUpMin(u)
            }
        }
    }

    private fun siftUpMin(c: Int) {
        // Min-side parent: (x/2-1)&~1
        var c = c
        var p = c / 2 - 1 and 1.inv()
        while (p >= 0 && keys[c] < keys[p]) {
            swap(c, p)
            c = p
            p = c / 2 - 1 and 1.inv()
        }
    }

    private fun siftUpMax(c: Int) {
        // Max-side parent: (x/2-1)|1
        var c = c
        var p = c / 2 - 1 or 1
        while (p >= 0 && keys[c] > keys[p]) {
            swap(c, p)
            c = p
            p = c / 2 - 1 or 1
        }
    }

    private fun siftDownMin(p: Int) {
        var p = p
        var c = p * 2 + 2
        while (c < size) {
            if (c + 2 < size && keys[c + 2] < keys[c]) {
                c += 2
            }
            if (keys[c] < keys[p]) {
                swap(p, c)
                // Swap with pair if necessary
                if (c + 1 < size && keys[c + 1] < keys[c]) {
                    swap(c, c + 1)
                }
            } else {
                break
            }
            p = c
            c = p * 2 + 2
        }
    }

    private fun siftDownMax(p: Int) {
        var p = p
        var c = p * 2 + 1
        while (c <= size) {
            if (c == size) {
                // If the left child only has half a pair
                if (keys[c - 1] > keys[p]) {
                    swap(p, c - 1)
                }
                break
            } else if (c + 2 == size) {
                // If there is only room for a right child lower pair
                if (keys[c + 1] > keys[c]) {
                    if (keys[c + 1] > keys[p]) {
                        swap(p, c + 1)
                    }
                    break
                }
            } else if (c + 2 < size) {
                // If there is room for a right child upper pair
                if (keys[c + 2] > keys[c]) {
                    c += 2
                }
            }
            if (keys[c] > keys[p]) {
                swap(p, c)
                // Swap with pair if necessary
                if (keys[c - 1] > keys[c]) {
                    swap(c, c - 1)
                }
            } else {
                break
            }
            p = c
            c = p * 2 + 1
        }
    }

    override fun size(): Int {
        return size
    }

    fun capacity(): Int {
        return capacity
    }

//    override fun toString(): String {
//        val twoPlaces = DecimalFormat("0.00")
//        val str: StringBuffer = StringBuffer(IntervalHeap::class.java.getCanonicalName())
//        str.append(", size: ").append(size()).append(" capacity: ").append(capacity())
//        var i = 0
//        var p = 2
//        while (i < size()) {
//            var x = 0
//            str.append("\t")
//            while (i + x < size() && x < p) {
//                str.append(twoPlaces.format(keys[i + x])).append(", ")
//                x++
//            }
//            str.append("\n")
//            i += x
//            p *= 2
//        }
//        return str.toString()
//    }

    private fun validateHeap(): Boolean {
        // Validate left-right
        run {
            var i = 0
            while (i < size - 1) {
                if (keys[i] > keys[i + 1]) {
                    return false
                }
                i += 2
            }
        }
        // Validate within parent interval
        for (i in 2 until size) {
            val maxParent = keys[i / 2 - 1 or 1]
            val minParent = keys[i / 2 - 1 and 1.inv()]
            if (keys[i] > maxParent || keys[i] < minParent) {
                return false
            }
        }
        return true
    }

    companion object {
        private const val defaultCapacity = 64
    }
}