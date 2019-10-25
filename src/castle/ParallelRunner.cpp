/*
 * ParallelRunner.cpp
 *
 *  Created on: Feb 20, 2015
 *      Author: Euncheon Lim @ Max-Planck-Institute for Developmental Biology
 *            : under the guidance of Prof. Detlef Weigel and Prof. Daniel Huson
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#include "ParallelRunner.hpp"

namespace castle {
function<void(size_t, size_t)> ParallelRunner::empty_ranged_func = [&](size_t, size_t){};
ParallelRunner::ParallelRunner() {

}

ParallelRunner::~ParallelRunner() {
}

} /* namespace castle */
