# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 12:44:30 2018

@author: bimta
"""
#The filter that grey_erosion is actually calling:
import numpy 

def _min_or_max_filter(input, size, footprint, structure, output, mode,
                       cval, origin, minimum):
    if structure is None:
        if footprint is None:
            if size is None:
                raise RuntimeError("no footprint provided")
            separable = True
        else:
            footprint = numpy.asarray(footprint, dtype=bool)
            if not footprint.any():
                raise ValueError("All-zero footprint is not supported.")
            if footprint.all():
                size = footprint.shape
                footprint = None
                separable = True
            else:
                separable = False
    else: #my case
        structure = numpy.asarray(structure, dtype=numpy.float64)
        separable = False
        if footprint is None: #my case
            footprint = numpy.ones(structure.shape, bool)
        else:
            footprint = numpy.asarray(footprint, dtype=bool)
    
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    output, return_value = _ni_support._get_output(output, input)
    origins = _ni_support._normalize_sequence(origin, input.ndim)

    if separable:
        sizes = _ni_support._normalize_sequence(size, input.ndim)
        modes = _ni_support._normalize_sequence(mode, input.ndim)
        axes = list(range(input.ndim))
        axes = [(axes[ii], sizes[ii], origins[ii], modes[ii])
                for ii in range(len(axes)) if sizes[ii] > 1]
        if minimum:
            filter_ = minimum_filter1d
        else:
            filter_ = maximum_filter1d
        if len(axes) > 0:
            for axis, size, origin, mode in axes:
                filter_(input, int(size), axis, output, mode, cval, origin)
                input = output
        else:
            output[...] = input[...]

    else:
        fshape = [ii for ii in footprint.shape if ii > 0]
        if len(fshape) != input.ndim:
            raise RuntimeError('footprint array has incorrect shape.')
        for origin, lenf in zip(origins, fshape):
            if (lenf // 2 + origin < 0) or (lenf // 2 + origin >= lenf):
                raise ValueError('invalid origin')
        if not footprint.flags.contiguous:
            footprint = footprint.copy()
        if structure is not None:
            if len(structure.shape) != input.ndim:
                raise RuntimeError('structure array has incorrect shape')
            if not structure.flags.contiguous:
                structure = structure.copy()
        mode = _ni_support._extend_mode_to_code(mode)
        _nd_image.min_or_max_filter(input, footprint, structure, output,
                                    mode, cval, origins, minimum)
    return return_value