package kopenvoronoi.util

class Pair<A, B>(first: A, second: B) {
    var first: A?
        private set
    var second: B?
        private set

    init {
        this.first = first
        this.second = second
    }

    override fun hashCode(): Int {
        val hashFirst = if (first != null) first.hashCode() else 0
        val hashSecond = if (second != null) second.hashCode() else 0
        return (hashFirst + hashSecond) * hashSecond + hashFirst
    }

    override fun equals(other: Any?): Boolean {
        if (other is Pair<*, *>) {
            val otherPair = other
            return ((first === otherPair.first || first != null && otherPair.first != null && first == otherPair.first)
                    && (second === otherPair.second || second != null && otherPair.second != null && second == otherPair.second))
        }
        return false
    }

    override fun toString(): String {
        return "($first, $second)"
    }

    fun setFirst(first: A) {
        this.first = first
    }

    fun setSecond(second: B) {
        this.second = second
    }
}
