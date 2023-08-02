fun assert(b: Boolean, function: () -> String) {
    if (!b) throw AssertionError(function())
}

fun assert(b: Boolean) {
    if (!b) throw AssertionError()
}