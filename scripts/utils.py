def is_ms_bh_pair(st1, st2):
    return ((st1 <= 9 and st2 == 14) or (st2 <= 9 and st1 == 14))

def is_bh_bh_pair(st1, st2):
    return (st1 == 14 and st2 == 14)

def check_merger_MT_hist(mt_history):
    for x in mt_history:
        if x == 6:
            return True
        else:
            return False

def check_merger_manual(r1, r2, sa):
    if sa >= 0:
        return ((r1 + r2) >= sa)
    else:
        return False
    
def which_bh(stellar_type_1, stellar_type_2):
    if stellar_type_1 == 14:
        return 1
    elif stellar_type_2 == 14:
        return 0
