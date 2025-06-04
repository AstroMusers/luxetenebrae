def is_ms_bh_pair(st1, st2):
    return ((st1 <= 9 and st2 == 14) or (st2 <= 9 and st1 == 14))