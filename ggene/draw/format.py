

def format_genomic(index):
    index = int(index)
    fstr = "0.1f"
    if index > 0.5e6:
        div = 1e6
        unit = "M"
    elif index > 0.5e3:
        div = 1e3
        unit = "k"
    else:
        div = 1
        unit = "b"
        fstr = ""
    return format(index/div, fstr) + unit


