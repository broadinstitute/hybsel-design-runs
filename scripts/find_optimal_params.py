#!/bin/python3
"""Optimize the parameter values across datasets to construct a probe set.

This chooses the optimal parameter values, with different values for each
dataset. It does so by treating the problem as a constrained nonlinear
optimization problem. In particular, it defines a loss over the parameter
values and seeks to minimize the loss, subject to the constraint that the
total number of probes (using those parameter values) is less than the
maximum number of allowed probes (an argument to the script). It enforces
this constraint using a barrier function added to the loss function, and
uses scipy's optimize module to minimize the sum of the loss and barrier
function.
"""

import argparse
from collections import defaultdict
import math

import numpy as np
from scipy import optimize

import utils


def round_up(x, b):
    """Round float x up to the nearest multiple of int b
    """
    return int(math.ceil(float(x) / b)) * b
def round_down(x, b):
    """Round float x down to the nearest multiple of int b
    """
    return int(math.floor(float(x) / b)) * b


def make_interp_probe_count_for_dataset_fn(probe_counts):

    memoized_bounding_boxes = {dataset: {} for dataset in probe_counts.keys()}
    def immediate_bounding_box(mismatches, cover_extension):
        """Return the box around the values rounded down/up to the nearest int
        to use when memoizing bounding boxes.

        Since the cover_extensions are multiples of 10, we could really round
        down/up cover_extension to the nearest 10. But just do the nearest
        integer here to make it easier to generalize.
        """
        return (round_down(mismatches, 1),
                round_up(mismatches, 1),
                round_down(cover_extension, 1),
                round_up(cover_extension, 1))

    def find_bounding_box_around_point(dataset, mismatches, cover_extension):
        """Return a rectangular bounding box around given parameters.

        Based on probe_counts[dataset], choose parameter values that form the
        smallest box around (mismatches, cover_extension) so that we can
        use actual computed probe counts (based on the box's parameter values)
        to interpolate the probe count for (mismatches, cover_extension).
        """
        # 'mismatches' may be a float between the min/max mismatches for
        # dataset. 'cover_extension' may be a float between the min/max
        # cover_extension for dataset. Find the mismatches and cover_extension
        # parameters for dataset that form the smallest rectangle encompassing
        # 'mismatches' and 'cover_extension', so that we can use the points of
        # this rectangle to perform bilinear interpolation. We brute force over
        # all possible rectangles encompassing 'mismatches' and
        # 'cover_extension' (with some heuristics to not be completely
        # exhaustive), which takes O(n^3) time; there are faster methods, but
        # this should suffice here.
        # Consider mismatches on the x-axis (left to right) and cover_extension
        # on the y-axis (bottom to top)
        points = set(probe_counts[dataset].keys())
        points_topleft = set()
        points_topright = set()
        points_bottomleft = set()
        points_bottomright = set()
        for p in points:
            m, ce = p
            if m == mismatches:
                if ce == cover_extension:
                    points_topleft.add(p)
                    points_topright.add(p)
                    points_bottomleft.add(p)
                    points_bottomright.add(p)
                elif ce > cover_extension:
                    points_topleft.add(p)
                    points_topright.add(p)
                else:
                    points_bottomleft.add(p)
                    points_bottomright.add(p)
            elif m > mismatches:
                if ce == cover_extension:
                    points_topright.add(p)
                    points_bottomright.add(p)
                elif ce > cover_extension:
                    points_topright.add(p)
                else:
                    points_bottomright.add(p)
            else:
                if ce == cover_extension:
                    points_topleft.add(p)
                    points_bottomleft.add(p)
                elif ce > cover_extension:
                    points_topleft.add(p)
                else:
                    points_bottomleft.add(p)

        points_topright_by_y = defaultdict(set)
        for p in points_topright:
            m, ce = p
            points_topright_by_y[ce].add(p)
        points_bottomleft_by_x = defaultdict(set)
        for p in points_bottomleft:
            m, ce = p
            points_bottomleft_by_x[m].add(p)

        min_rectangle, min_area = None, float('inf')
        for p_topleft in points_topleft:
            p_topleft_m, p_topleft_ce = p_topleft
            # find all points in the top-right with the same cover_extension
            # (y-value) as p_topleft
            for p_topright in points_topright_by_y[p_topleft_ce]:
                p_topright_m, p_topright_ce = p_topright
                # find all points in the bottom-left with the same mismatches
                # (x-value) as p_topleft
                for p_bottomleft in points_bottomleft_by_x[p_topleft_m]:
                    p_bottomleft_m, p_bottomleft_ce = p_bottomleft
                    # to form a rectangle, we need the point (p_topright_m,
                    # p_bottomleft_ce); this point should be in the
                    # bottom-right, so check if it exists
                    p_bottomright = (p_topright_m, p_bottomleft_ce)
                    if p_bottomright in points_bottomright:
                        # we found a valid rectangle; now compute its 'area'
                        width = p_topright_m - p_topleft_m
                        height = (p_topright_ce - p_bottomleft_ce) / 10.0
                        # add pseudocounts to width and height because if
                        # width or height is 0, we still want the other
                        # dimension to be accounted for
                        area = (width + 0.001) * (height + 0.001)
                        if area < min_area:
                            min_rectangle = (p_topleft, p_bottomright)
                            min_area = area
        return min_rectangle

    def interp_probe_count_for_dataset(dataset, mismatches, cover_extension):
        """
        Using the given probe counts at particular parameter values, interpolate
        the number of probes for 'dataset' and 'mismatches' mismatches and
        at a cover extension of 'cover_extension', where each of these may be
        floats
        """
        immediate_bb = immediate_bounding_box(mismatches, cover_extension)
        if immediate_bb in memoized_bounding_boxes[dataset]:
            # The bounding box for (mismatches, cover_extension) has been
            # memoized
            min_rectangle = memoized_bounding_boxes[dataset][immediate_bb]
        else:
            # Compute and memoized the bounding box for (mismatches,
            # cover_extension)
            min_rectangle = find_bounding_box_around_point(dataset,
                                                           mismatches,
                                                           cover_extension)
            if min_rectangle is None:
                raise ValueError(("Unable to find rectangular bounding box around "
                                  "(mismatches, cover_extension)=(%f, %f) for "
                                  "dataset %s") % (mismatches, cover_extension,
                                  dataset))
            memoized_bounding_boxes[dataset][immediate_bb] = min_rectangle

        rect_topleft, rect_bottomright = min_rectangle
        mismatches_floor, cover_extension_ceil = rect_topleft
        mismatches_ceil, cover_extension_floor = rect_bottomright

        # At cover_extension_floor and at cover_extension_ceil, interpolate the
        # number of probes at 'mismatches' mismatches (interpolate linearly)
        for ce in [cover_extension_floor, cover_extension_ceil]:
            count_left = probe_counts[dataset][(mismatches_floor, ce)]
            count_right = probe_counts[dataset][(mismatches_ceil, ce)]
            mismatches_diff = mismatches_ceil - mismatches_floor
            if mismatches_diff == 0:
                # count_left should equal count_right
                assert count_left == count_right
                count = count_left
            elif count_left <= count_right:
                count_diff = count_right - count_left
                f = float(mismatches - mismatches_floor) / mismatches_diff
                count = f * count_diff + count_left
            else:
                count_diff = count_left - count_right
                f = float(mismatches - mismatches_floor) / mismatches_diff
                count = count_left - f * count_diff
            if ce == cover_extension_floor:
                count_floor = count
            if ce == cover_extension_ceil:
                count_ceil = count

        # Interpolate between cover_extension_floor and cover_extension_ceil
        # using count_floor and count_ceil (interpolate linearly)
        cover_extension_diff = cover_extension_ceil - cover_extension_floor
        if cover_extension_diff == 0:
            # count_floor should equal count_ceil
            assert count_floor == count_ceil
            final_interp = count_floor
        elif count_floor <= count_ceil:
            count_diff = count_ceil - count_floor
            f = float(cover_extension - cover_extension_floor) / cover_extension_diff
            final_interp = f * count_diff + count_floor
        else:
            count_diff = count_floor - count_ceil
            f = float(cover_extension - cover_extension_floor) / cover_extension_diff
            final_interp = count_floor - f * count_diff

        return final_interp

    return interp_probe_count_for_dataset


def make_total_probe_count_across_datasets_fn(probe_counts):
    interp_probe_count_for_dataset = make_interp_probe_count_for_dataset_fn(
        probe_counts)

    def total_probe_count_across_datasets(x):
        """
        Sum the (interpolated) probe counts across datasets.

        x is a list giving all the parameter values across datasets,
        such that for even i, x_i gives the number of mismatches for the
        (i/2)'th dataset and x_{i+1} gives the cover extension for the
        (i/2)'th dataset
        """
        s = 0
        for i, dataset in enumerate(sorted(probe_counts.keys())):
            mismatches, cover_extension = x[2 * i], x[2 * i + 1]
            s += interp_probe_count_for_dataset(dataset, mismatches,
                                                cover_extension)
        return s

    return total_probe_count_across_datasets


def make_loss_fn(probe_counts, max_probe_count):
    total_probe_count_across_datasets = make_total_probe_count_across_datasets_fn(
        probe_counts)

    def loss(x, *func_args):
        """
        Compute a loss.

        x is a list giving all the parameter values across datasets,
        such that for even i, x_i gives the number of mismatches for the
        (i/2)'th dataset and x_{i+1} gives the cover extension for the
        (i/2)'th dataset
        """
        # First compute a loss over the parameters by taking their L2-norm (and
        # down-weighting cover_extension by a factor of 10.0)
        # This is the function we really want to minimize
        opt_val = 0
        for i, dataset in enumerate(sorted(probe_counts.keys())):
            mismatches, cover_extension = x[2 * i], x[2 * i + 1]
            opt_val += np.power(mismatches, 2.0) + np.power(cover_extension / 10.0, 2.0)

        # We also have the constraint that the total probe count be less than
        # max_probe_count
        # We add a barrier function to enforce this constraint and weight the
        # barrier by eps
        eps = func_args[0]
        total_probe_count = total_probe_count_across_datasets(x)
        if total_probe_count >= max_probe_count:
            # Since the count is beyond the barrier, we should in theory
            # return infinity. But if the optimizer does indeed try parameters
            # that put the probe count here, it would be unable to compute
            # an approximate gradient and may get stuck. So help it out
            # by giving a value such that the negative gradient points toward
            # a direction outside the barrier.
            barrier_val = 9999 + 10.0 * np.log(total_probe_count)
        else:
            # The barrier function is -log(max_probe_count - total_probe_count), to
            # enforce the constraint that total_probe_count be less than
            # max_probe_count.
            barrier_val = -1.0 * eps * np.log(max_probe_count - total_probe_count)

        return opt_val + barrier_val

    return loss


def make_param_bounds(probe_counts, step_size=0.001):
    bounds = []
    for dataset in sorted(probe_counts.keys()):
        params = probe_counts[dataset].keys()

        # bound cover_extensions by the lowest and highest value for
        # which we have a probe count result
        cover_extensions = [k[1] for k in params]
        cover_extensions_lo = min(cover_extensions)
        cover_extensions_hi = max(cover_extensions)

        # to ensure we can find a rectangular bounding box around an
        # arbitrary point ((mismatches, cover_extension)), our lower
        # bound on mismatches should have a corresponding cover extension
        # of min(cover_extensions) and of max(cover_extensions);
        # so should our upper bound on mismatches
        mismatches = [k[0] for k in params]
        mismatches_with_valid_cover_extension = \
            [m for m in mismatches if ((m, cover_extensions_lo) in params and
             (m, cover_extensions_hi) in params)]
        mismatches_lo = min(mismatches_with_valid_cover_extension)
        mismatches_hi = max(mismatches_with_valid_cover_extension)

        bounds += [(mismatches_lo, mismatches_hi - step_size)]
        bounds += [(min(cover_extensions), max(cover_extensions) - step_size)]
    return bounds


def make_initial_guess(probe_counts, max_probe_count):
    # Guess mismatches=5, cover_extension=30 for all datasets
    x0 = np.array([5, 30] * len(probe_counts))

    # Verify that this yields fewer probes than the maximum allowed
    # (i.e., is not beyond the barrier)
    guess_probe_count = make_total_probe_count_across_datasets_fn(probe_counts)(x0)
    if guess_probe_count >= max_probe_count:
        raise ValueError(("Initial guess yields too many probes (%d, but the max "
                         "is %d)") % (guess_probe_count, max_probe_count))

    return x0


def optimize_loss(probe_counts, loss_fn, bounds, x0,
                  initial_eps=10.0, step_size=0.001):
    # Keep minimizing loss_fn while decreasing eps (so that the weight
    # of the barrier function is decreased until it is very small).
    # On each iteration, start the initial guess/position at the solution
    # of the previous iteration.

    eps = initial_eps
    while eps >= 0.01:
        x0_probe_count = make_total_probe_count_across_datasets_fn(probe_counts)(x0)
        print("Starting an iteration with eps=%f, with x0 yielding %d probes" % \
              (eps, x0_probe_count))

        sol, nfeval, rc = optimize.fmin_tnc(loss_fn, x0, bounds=bounds,
                                            args=(eps,),
                                            approx_grad=True,
                                            epsilon=step_size, disp=1, maxfun=2500)

        if rc in [0, 1, 2]:
            # rc == 0 indicates reaching the local minimum, and rc == 1 or
            # rc == 2 indicates the function value converged
            print("  Iteration was successful")
        else:
            print("  Iteration failed to converge!")

        x0 = sol
        eps = 0.1 * eps

    return sol


def total_probe_count_without_interp(params, probe_counts):
    """
    The result of make_total_probe_count_across_datasets_fn should give
    the same count as this function, assuming that params are keys in
    the datasets of probe_counts. But this uses probe_counts directly
    as a sanity check (i.e., does not do any interpolation).
    """
    s = 0
    for i, dataset in enumerate(sorted(probe_counts.keys())):
        mismatches, cover_extension = params[2 * i], params[2 * i + 1]
        s += probe_counts[dataset][(mismatches, cover_extension)]
    return s


def round_params(params, probe_counts, max_probe_count,
        mismatches_eps=0.01, cover_extension_eps=0.1):
    # Params are floats. We want the mismatches parameters to be integers
    # and the cover_extension parameters to be multiples of 10.
    #
    # The floats, as given in params, should satisfy the constraint (i.e.,
    # the interpolated total number of probes is less than max_probe_count).
    # Thus, we can round them up, because (generally) increasing the parameter
    # values will decrease the number of probes; therefore, after rounding up
    # they should still satisfy the constraint.
    #
    # But we also check if the parameter values are within eps of their
    # rounded-down value. The loss optimizer has a tendency to make this happen
    # for some parameters (e.g., finding an optimal mismatches parameter value
    # of 1.00001). The reason likely has to do with the fact that, because
    # we are linearly interpolating total probe counts, the gradient of the
    # barrier function changes greatly around certain values (i.e., around
    # the actual data values). That is, the barrier function is not at all
    # smooth around actual data values. This may cause the optimizer to
    # yield parameter values that are very close to parameter values for
    # which probe counts have actually been computed.
    #
    # After rounding up, some parameters are decreased; we repeatedly
    # choose to decrease the parameter whose reduction yields the smallest
    # loss while still yielding a number of probes that is less than
    # max_probe_count.

    params_rounded = []
    for i, dataset in enumerate(sorted(probe_counts.keys())):
        mismatches, cover_extension = params[2 * i], params[2 * i + 1]

        if mismatches - round_down(mismatches, 1) < mismatches_eps:
            # Round mismatches down
            mismatches = round_down(mismatches, 1)
        else:
            # Round mismatches up
            mismatches = round_up(mismatches, 1)

        if cover_extension - round_down(cover_extension, 10) < cover_extension_eps:
            # Round cover_extension down
            cover_extension = round_down(cover_extension, 10)
        else:
            # Round cover_extension up
            cover_extension = round_up(cover_extension, 10)

        params_rounded += [mismatches, cover_extension]

    total_probe_count = make_total_probe_count_across_datasets_fn(probe_counts)
    # Verify that the probe count satisfies the constraint
    # Note that this assertion may fail if we are dealing with datasets
    # for which few actual probe counts have been computed; in these
    # cases, the interpolation may severely underestimate the number
    # of probes at a particular parameter choice
    assert total_probe_count(params_rounded) < max_probe_count

    # Keep decreasing parameters while satisfying the constraint.
    # In particular, choose to decrease the parameter whose reduction
    # yields the smallest loss while still satisfying the constraint.
    loss_fn = make_loss_fn(probe_counts, max_probe_count)
    while True:
        curr_loss = loss_fn(params_rounded, 0)
        # Find a parameter to decrease
        min_loss, min_loss_new_params = curr_loss, None
        for i in range(len(params_rounded)):
            params_tmp = list(params_rounded)
            if params_tmp[i] == 0:
                # This cannot be decreased
                continue
            if i % 2 == 0:
                # This is a mismatch, so decrease by 1
                params_tmp[i] -= 1
            else:
                # This is a cover_extension, so decrease by 10
                params_tmp[i] -= 10
            if total_probe_count(params_tmp) >= max_probe_count:
                # This change yields too many probes, so skip it
                continue
            new_loss = loss_fn(params_tmp, 0)
            if new_loss < min_loss:
                min_loss = new_loss
                min_loss_new_params = params_tmp

        if min_loss_new_params != None:
            # There was a change that led to a better loss, so
            # update params_rounded
            params_rounded = min_loss_new_params
        else:
            # No parameter change satisfies the constraint and
            # yields an improvement in the loss
            break

    return params_rounded


def print_params_by_dataset(params, probe_counts, type="float"):
    for i, dataset in enumerate(sorted(probe_counts.keys())):
        mismatches, cover_extension = params[2 * i], params[2 * i + 1]
        if type == "float":
            print("%s: (%f, %f)" % (dataset, mismatches, cover_extension))
        elif type == "int":
            print("%s: (%d, %d)" % (dataset, mismatches, cover_extension))
        else:
            raise ValueError("Unknown type %s", type)


def write_params_to_file(params, probe_counts, path, type="int"):
    lines = []
    for i, dataset in enumerate(sorted(probe_counts.keys())):
        mismatches, cover_extension = params[2 * i], params[2 * i + 1]
        if type == "float":
            lines += ["%s\t(%f, %f)" % (dataset, mismatches, cover_extension)]
        elif type == "int":
            lines += ["%s\t(%d, %d)" % (dataset, mismatches, cover_extension)]
        else:
            raise ValueError("Unknown type %s", type)

    with open(path, 'w') as f:
        for line in lines:
            f.write(line + '\n')


def main(args):
    probe_counts = utils.read_probe_counts(args)

    loss_fn = make_loss_fn(probe_counts, args.max_probe_count)
    bounds = make_param_bounds(probe_counts)
    x0 = make_initial_guess(probe_counts, args.max_probe_count)

    x_sol = optimize_loss(probe_counts, loss_fn, bounds, x0)

    print("##############################")
    print("Continuous parameter values:")
    print_params_by_dataset(x_sol, probe_counts, "float")
    x_sol_count = make_total_probe_count_across_datasets_fn(probe_counts)(x_sol)
    print("TOTAL INTERPOLATED PROBE COUNT: %d" %  x_sol_count)
    print("##############################")
    print()

    opt_params = round_params(x_sol, probe_counts, args.max_probe_count)

    print("##############################")
    print("Rounded parameter values:")
    print_params_by_dataset(opt_params, probe_counts, "int")
    opt_params_count = make_total_probe_count_across_datasets_fn(probe_counts)(opt_params)
    print("TOTAL PROBE COUNT: %d" % opt_params_count)
    print("##############################")

    if args.verify_without_interp:
        # As a sanity check, verify that we get the same total probe count
        # without using interpolation
        assert opt_params_count == total_probe_count_without_interp(opt_params, probe_counts)

    if args.output_params:
        write_params_to_file(opt_params, probe_counts, args.output_params)


if __name__ == "__main__":
    argparse = argparse.ArgumentParser()
    argparse.add_argument('--results_dir', '-i', required=True)
    argparse.add_argument('--limit_datasets', '-d', nargs='+')
    argparse.add_argument('--max_probe_count', '-n', type=int, default=90000)
    argparse.add_argument('--output_params', '-o')
    argparse.add_argument('--verify_without_interp',
                          dest='verify_without_interp',
                          action='store_true')
    args = argparse.parse_args()

    main(args)
