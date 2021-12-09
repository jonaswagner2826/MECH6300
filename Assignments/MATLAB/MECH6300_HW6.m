% MECH 6300 - HW 6

A = blkdiag([2,1;0,2],2,2,[1,1;0,1],1)
B = [2  1   0
    2   1   1
    1   1   1
    3   2   1
    -1  0   1
    1   0   1
    1   0   0]
C = [2  2   1   3   -1  1   1
    1   1   1   2   0   0   0
    0   1   1   1   1   1   0]

U = ctrb(A,B)

ctrlRank = rank(U)

V = obsv(A,C)

obsvRank = rank(V)