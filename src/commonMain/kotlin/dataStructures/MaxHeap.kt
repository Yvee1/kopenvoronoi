package dataStructures

/**
 *
 */
interface MaxHeap<T> {
    fun size(): Int
    fun offer(key: Double, value: T)
    fun replaceMax(key: Double, value: T)
    fun removeMax()
    val max: T
    val maxKey: Double
}