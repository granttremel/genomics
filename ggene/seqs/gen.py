

import random





def get_random_sequence(seq_len, ab = "ATGC"):
    return "".join(random.choice(ab) for n in range(seq_len))

