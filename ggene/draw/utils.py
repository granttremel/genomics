

def quantize(data, bit_depth, maxval=None, minval=None, mid = False):
    if not maxval:
        maxval = max(data)
    if not minval:
        minval = min(data)
    rng = maxval-minval
    c = (minval+maxval)/2
    off = 0 if mid else 0.5
    return [int(bit_depth * (((d-c)/rng) + off)) for d in data]

