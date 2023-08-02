package dataStructures

/**
 *
 */
interface MinHeap<T> {
    fun size(): Int
    fun offer(key: Double, value: T)
    fun replaceMin(key: Double, value: T)
    fun removeMin()
    val min: T
    val minKey: Double
}
