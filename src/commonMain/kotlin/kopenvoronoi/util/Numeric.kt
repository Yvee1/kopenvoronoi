package kopenvoronoi.util

import assert
import kotlin.math.abs
import kotlin.math.sqrt

/**
 * Holds general numerical functions that are not specific to voronoi-diagrams
 * and may be useful elsewhere too
 */
object Numeric {
    // solve quadratic eqn: a*x*x + b*x + c = 0
    // returns real roots (0, 1, or 2) as vector
    fun quadratic_roots(a: Double, b: Double, c: Double): List<Double> {
        val roots: MutableList<Double> = mutableListOf()
        if (a == 0.0 && b == 0.0) {
            return roots
        }
        if (a == 0.0) {
            roots.add(-c / b)
            return roots
        }
        if (b == 0.0) {
            val sqr = -c / a
            return if (sqr > 0) {
                roots.add(sqrt(sqr))
                roots.add(-roots[0])
                roots
            } else if (sqr == 0.0) {
                roots.add(0.0)
                roots
            } else {
                // std::cout << " quadratic_roots() b == 0. no roots.\n";
                roots
            }
        }
        val disc = chop(b * b - 4 * a * c) // discriminant, chop!
        if (disc > 0) {
            val q: Double
            q = if (b > 0) {
                (b + sqrt(disc)) / -2
            } else {
                (b - sqrt(disc)) / -2
            }
            roots.add(q / a)
            roots.add(c / q)
            return roots
        } else if (disc == 0.0) {
            roots.add(-b / (2 * a))
            return roots
        }
        // std::cout << " quadratic_roots() disc < 0. no roots. disc= " << disc << "\n";
        return roots
    }

    fun determinant(
        a: Double, b: Double, c: Double, d: Double, e: Double, f: Double, g: Double, h: Double,
        i: Double
    ): Double {
        return a * (e * i - h * f) - b * (d * i - g * f) + c * (d * h - g * e)
    }

    fun sq(a: Double): Double {
        return a * a
    }

    fun chop(`val`: Double, tol: Double): Double {
        return if (abs(`val`) < tol) {
            0.0
        } else {
            `val`
        }
    }

    fun chop(`val`: Double): Double {
        val _epsilon = 1e-10
        return if (abs(`val`) < _epsilon) {
            0.0
        } else {
            `val`
        }
    }

    fun diangle(x: Double, y: Double): Double {
        return if (y >= 0) {
            if (x >= 0) y / (x + y) else 1 - x / (-x + y)
        } else {
            if (x < 0) 2 - y / (-x - y) else 3 + x / (x - y)
        }
    }

    fun diangle_x(a: Double): Double {
        return if (a < 2) 1 - a else a - 3
    }

    fun diangle_y(a: Double): Double {
        return if (a < 3) (if (a > 1) 2 - a else a) else a - 4
    }

    fun diangle_xy(a: Double): Pair<Double, Double> {
        val x = diangle_x(a)
        val y = diangle_y(a)
        val norm: Double = sqrt(x * x + y * y)
        return Pair(x / norm, y / norm)
    }

    // return true if a lies in [less,more]
    fun diangle_bracket(less: Double, a: Double, more: Double): Boolean {
        if (less == more) {
            return false
        } else if (less <= more) { // normal case..
            if (less <= a && a < more) {
                return true
            }
        } else if (less <= a && a <= 4 || 0 <= a && a < more) {
            return true
        }
        return false
    }

    // return average of input angles
    fun diangle_mid(alfa1: Double, alfa2: Double): Double {
        return if (alfa1 <= alfa2) {
            (alfa1 + alfa2) / 2
        } else {
            val opposite_mid = alfa2 + (alfa1 - alfa2) / 2
            var mid = opposite_mid + 2
            if (mid > 4) {
                mid = mid - 4
            }
            assert(0 <= mid && mid <= 4) { " (0<=mid) && (mid<=4) " }
            mid
        }
    }
}
