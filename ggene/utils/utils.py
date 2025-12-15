







def format_genomic(index):
    
    if index > 0.5e6:
        div = 1e6
        unit = "M"
    elif index > 0.5e3:
        div = 1e3
        unit = "k"
    else:
        div = 1
        unit = "b"
    return format(index/div, ".1f") + unit
