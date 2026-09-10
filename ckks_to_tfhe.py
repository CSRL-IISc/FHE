#!/usr/bin/env python3
"""
=====================================================================
 CKKS  ->  FHEW/TFHE  SCHEME SWITCHING, FROM SCRATCH  (no OpenFHE)
=====================================================================

This file re-implements, with toy parameters, exactly the pipeline that
OpenFHE runs inside:

    cc->EvalCKKStoFHEWSetup(...)      -> Section 0/7  (parameters, LWE ctx)
    cc->EvalCKKStoFHEWKeyGen(...)     -> Section 8    (RLWE->LWE switching key)
    cc->EvalCKKStoFHEWPrecompute(s)   -> Section 6    (scaled SlotsToCoeffs matrix)
    cc->EvalCKKStoFHEW(ct, n)         -> Section 6+8+9+10

Parameters are deliberately tiny (N = 8) so every intermediate value can
be printed and checked by hand.  There is NO security here.

Conventions used throughout
---------------------------
Ring          R_Q = Z_Q[x] / (x^N + 1)
CKKS ct       (c0, c1),  phase(ct) = c0 + c1 * s   mod Q
LWE  ct       (a, b),    phase(ct) = b + <a, s>    mod q
(TFHE papers usually write b - <a,s>; a one-line converter is provided
 at the very bottom -- it just negates a.)
"""

import math, cmath, random

random.seed(2024)

# =====================================================================
# 0.  PARAMETERS
# =====================================================================
N        = 8                 # RLWE ring dimension (x^N + 1)
M        = 2 * N             # order of the root of unity  (16)
NS       = N // 2            # number of CKKS slots        (4)

LOGDELTA = 35
DELTA    = 1 << LOGDELTA     # CKKS scaling factor  2^35
LOGQ0    = 45
Q0       = 1 << LOGQ0        # CKKS *last* modulus (level 0)   2^45
Q1       = Q0 * DELTA        # CKKS *fresh* modulus (level 1)  2^80

BKS      = 1 << 5            # key-switching digit base (32).
                             # Q0 and Q1 are exact powers of BKS -> the
                             # balanced digit decomposition is exact mod Q.
SIGMA    = 3.2               # error stddev

n_lwe    = 4                 # LWE dimension  (must be <= N)
q_lwe    = 1 << 10           # LWE modulus    (1024)
p_lwe    = 4                 # LWE plaintext modulus -> messages {0,1,2,3}

# The magic constant folded into the SlotsToCoeffs matrix.
# After S2C + rescale the coefficient holds  DELTA * CONST * m.
# We want it to hold  (Q0/p) * m  so that mod-switching Q0 -> q_lwe
# lands exactly on (q_lwe/p) * m, the standard LWE encoding.
CONST = Q0 // (p_lwe * DELTA)          # = 2^45 / (4 * 2^35) = 256


def banner(t):
    print("\n" + "=" * 68); print(t); print("=" * 68)


# =====================================================================
# 1.  RING ARITHMETIC IN  Z_Q[x]/(x^N+1)
# =====================================================================
def center(x, Q):
    """representative of x in (-Q/2, Q/2]"""
    x %= Q
    return x - Q if x > Q // 2 else x

def round_div(x, d):
    """nearest integer to x/d, exact integer arithmetic (x may be negative)"""
    return (2 * x + d) // (2 * d)

def poly_add(a, b, Q): return [(x + y) % Q for x, y in zip(a, b)]
def poly_sub(a, b, Q): return [(x - y) % Q for x, y in zip(a, b)]

def poly_mul(a, b, Q):
    """schoolbook negacyclic convolution: x^N = -1"""
    r = [0] * N
    for i, ai in enumerate(a):
        if ai == 0: continue
        for j, bj in enumerate(b):
            k = i + j
            if k < N: r[k] += ai * bj
            else:     r[k - N] -= ai * bj
    return [x % Q for x in r]

def poly_scal(a, c, Q): return [(x * c) % Q for x in a]
def poly_mod (a, Q):    return [x % Q for x in a]

def sample_ternary():  return [random.choice([-1, 0, 1]) for _ in range(N)]
def sample_gauss():    return [int(round(random.gauss(0, SIGMA))) for _ in range(N)]
def sample_unif(Q):    return [random.randrange(Q) for _ in range(N)]


# =====================================================================
# 2.  CKKS ENCODING / DECODING  (canonical embedding)
# =====================================================================
# zeta = exp(2*pi*i/M).  The Galois group of Q(zeta)/Q is {+-5^j}.
# Slot j of a plaintext polynomial p(x) is  p(zeta^(5^j)).
ROOT = [cmath.exp(2j * math.pi * t / M) for t in range(M)]
GAL  = [pow(5, j, M) for j in range(NS)]          # [1, 5, 9, 13]

def ckks_encode(zvec, scale, Q):
    """complex slot vector (length NS) -> integer polynomial mod Q"""
    out = []
    for k in range(N):
        acc = 0j
        for j in range(NS):
            acc += zvec[j] * ROOT[(-GAL[j] * k) % M]   # conj(zeta^{G_j k})
        out.append(int(round(scale * 2 * acc.real / N)) % Q)
    return out

def ckks_decode(poly, scale, Q):
    """integer polynomial mod Q -> complex slot vector"""
    c = [center(x, Q) for x in poly]
    res = []
    for j in range(NS):
        acc = 0j
        for k in range(N):
            acc += c[k] * ROOT[(GAL[j] * k) % M]
        res.append(acc / scale)
    return res


# =====================================================================
# 3.  CKKS KEYGEN / ENCRYPT / DECRYPT
# =====================================================================
def ckks_keygen():
    return sample_ternary()

def ckks_encrypt(pt, s, Q):
    a  = sample_unif(Q)
    e  = poly_mod(sample_gauss(), Q)
    c0 = poly_add(poly_sub(pt, poly_mul(a, poly_mod(s, Q), Q), Q), e, Q)
    return (c0, a)                       # phase = c0 + c1*s = pt + e

def ckks_phase(ct, s, Q):
    """raw noisy plaintext polynomial  c0 + c1*s  mod Q"""
    return poly_add(ct[0], poly_mul(ct[1], poly_mod(s, Q), Q), Q)

def ckks_decrypt(ct, s, Q, scale):
    return ckks_decode(ckks_phase(ct, s, Q), scale, Q)

def ct_pt_mul(ct, pt, Q):
    return (poly_mul(ct[0], pt, Q), poly_mul(ct[1], pt, Q))

def ct_add(ct1, ct2, Q):
    return (poly_add(ct1[0], ct2[0], Q), poly_add(ct1[1], ct2[1], Q))

def rescale(ct, Q, factor):
    """divide ciphertext (and modulus) by `factor`, rounding"""
    Qn = Q // factor
    f = lambda c: [round_div(center(x, Q), factor) % Qn for x in c]
    return (f(ct[0]), f(ct[1])), Qn


# =====================================================================
# 4.  GADGET (DIGIT-DECOMPOSITION) KEY SWITCHING
# =====================================================================
def ks_len(Q, B):
    l, t = 0, 1
    while t < Q:
        t *= B; l += 1
    return l

def gen_ksk(s_from, s_to, Q, B):
    """
    KSK_i = (b_i, a_i) with  b_i = -a_i*s_to + B^i * s_from + e_i
    so that  b_i + a_i*s_to = B^i*s_from + e_i.
    """
    l, ksk = ks_len(Q, B), []
    sf, st = poly_mod(s_from, Q), poly_mod(s_to, Q)
    for i in range(l):
        a = sample_unif(Q)
        e = poly_mod(sample_gauss(), Q)
        b = poly_add(poly_sub(e, poly_mul(a, st, Q), Q),
                     poly_scal(sf, pow(B, i), Q), Q)
        ksk.append((b, a))
    return ksk

def decompose(c, Q, B, l):
    """balanced base-B digits: c = sum_i d_i B^i (mod Q), d_i in [-B/2, B/2)"""
    digits = [[0] * N for _ in range(l)]
    for k in range(N):
        x = c[k] % Q
        for i in range(l):
            d = x % B
            x //= B
            if d >= B // 2:
                d -= B; x += 1
            digits[i][k] = d
    return digits

def key_switch(ct, ksk, Q, B):
    """(c0,c1) under s_from  ->  (c0',c1') under s_to, same phase"""
    l  = len(ksk)
    dg = decompose(ct[1], Q, B, l)
    c0, c1 = ct[0][:], [0] * N
    for i in range(l):
        b, a = ksk[i]
        c0 = poly_add(c0, poly_mul(dg[i], b, Q), Q)
        c1 = poly_add(c1, poly_mul(dg[i], a, Q), Q)
    return (c0, c1)


# =====================================================================
# 5.  GALOIS AUTOMORPHISMS = SLOT ROTATIONS
# =====================================================================
def auto_poly(a, g, Q):
    """p(x) -> p(x^g),  g odd, using x^N = -1"""
    r = [0] * N
    for k in range(N):
        idx = (k * g) % M
        if idx < N: r[idx] = (r[idx] + a[k]) % Q
        else:       r[idx - N] = (r[idx - N] - a[k]) % Q
    return r

def gen_rot_keys(s, Q, B):
    """rotation key for step r: switches from  s(x^{5^r})  back to  s(x)"""
    keys = {}
    for r in range(1, NS):
        g = pow(5, r, M)
        keys[r] = (g, gen_ksk(auto_poly(poly_mod(s, Q), g, Q), s, Q, B))
    return keys

def rotate(ct, r, rotkeys, Q, B):
    """slots  ->  [z_r, z_{r+1}, ..., ] (cyclic left rotation by r)"""
    if r == 0: return ct
    g, ksk = rotkeys[r]
    ct_g = (auto_poly(ct[0], g, Q), auto_poly(ct[1], g, Q))
    return key_switch(ct_g, ksk, Q, B)


# =====================================================================
# 6.  SlotsToCoeffs   (== EvalCKKStoFHEWPrecompute + the linear transform)
# =====================================================================
# Goal: the values sitting in the *slots* must be moved into the
# *coefficients* of the plaintext polynomial, because LWE extraction can
# only reach coefficients.
#
# We want a plaintext p' with coeff vector d = (z0,z1,z2,z3, 0,0,0,0).
# Its slot vector is  slots(p') = U d  where U[j][k] = zeta^{G_j k}.
# Since d is zero above index NS,  slots(p') = A z  with the NS x NS matrix
#        A[j][k] = zeta^{G_j * k},   j,k = 0..NS-1
# So SlotsToCoeffs is just the homomorphic linear map  z -> CONST * A z.
def s2c_diagonals(const):
    A = [[ROOT[(GAL[j] * k) % M] for k in range(NS)] for j in range(NS)]
    # diag_r(A)[j] = A[j][(j+r) mod NS]   ->   A z = sum_r diag_r . rot(z,r)
    return [[const * A[j][(j + r) % NS] for j in range(NS)] for r in range(NS)]

def slots_to_coeffs(ct, diags, rotkeys, Q, B, scale):
    acc = ([0] * N, [0] * N)
    for r in range(NS):
        ct_r = rotate(ct, r, rotkeys, Q, B)
        pt_r = ckks_encode(diags[r], scale, Q)
        acc  = ct_add(acc, ct_pt_mul(ct_r, pt_r, Q), Q)
    return acc


# =====================================================================
# 7.  LWE / TFHE SIDE
# =====================================================================
def lwe_keygen(n):
    return [random.randrange(2) for _ in range(n)]     # binary secret

def lwe_phase(ct, s, q):
    a, b = ct
    return (b + sum(ai * si for ai, si in zip(a, s))) % q

def lwe_decrypt(ct, s, q, p):
    return round_div(lwe_phase(ct, s, q) * p, q) % p

def lwe_add(c1, c2, q):
    return ([(x + y) % q for x, y in zip(c1[0], c2[0])], (c1[1] + c2[1]) % q)


# =====================================================================
# 8.  RLWE -> LWE KEY SWITCHING KEY   (== EvalCKKStoFHEWKeyGen)
# =====================================================================
# The FHEW/TFHE secret s_lwe (dimension n) is *embedded* into a ring
# element  z(x) = s_lwe[0] + s_lwe[1] x + ... + s_lwe[n-1] x^{n-1}
# (zero-padded up to degree N-1).  We then build a normal RLWE key-
# switching key from s_ckks to z.  After switching, every coefficient of
# the ciphertext is already an LWE sample under s_lwe.
def embed_lwe_key(s_lwe):
    return [s_lwe[i] if i < len(s_lwe) else 0 for i in range(N)]


# =====================================================================
# 9.  COEFFICIENT EXTRACTION  (RLWE sample extraction)
# =====================================================================
# For c1 * z in R_Q with x^N = -1:
#   (c1*z)_k = sum_{j<=k} c1[k-j] z[j]  -  sum_{j>k} c1[N+k-j] z[j]
# so the LWE mask for coefficient k is
#   a[j] =  c1[k-j]        for j <= k
#   a[j] = -c1[N+k-j]      for j >  k
# truncated to the first n entries (the rest multiply zeros of z).
def extract_lwe(ct, k, n, Q):
    c0, c1 = ct
    a = [c1[k - j] % Q if j <= k else (-c1[N + k - j]) % Q for j in range(n)]
    return (a, c0[k] % Q)


# =====================================================================
# 10. MODULUS SWITCHING  Q0 -> q_lwe
# =====================================================================
def lwe_modswitch(ct, Qf, Qt):
    f = lambda x: round_div(center(x, Qf) * Qt, Qf) % Qt
    a, b = ct
    return ([f(x) for x in a], f(b))


# =====================================================================
# ---------------------------  SELF TESTS  ----------------------------
# =====================================================================
def self_tests():
    banner("SELF TESTS (encode/decode, encryption, rotation, key switch)")
    z = [1.5, -2.25, 0.75, 3.0]
    pt = ckks_encode(z, DELTA, Q1)
    back = ckks_decode(pt, DELTA, Q1)
    err = max(abs(a - b) for a, b in zip(z, back))
    print(f"encode/decode roundtrip error : {err:.3e}")

    s = ckks_keygen()
    ct = ckks_encrypt(pt, s, Q1)
    dec = ckks_decrypt(ct, s, Q1, DELTA)
    print(f"encrypt/decrypt error         : "
          f"{max(abs(a-b) for a,b in zip(z,dec)):.3e}")

    rk = gen_rot_keys(s, Q1, BKS)
    for r in [1, 2, 3]:
        cr = rotate(ct, r, rk, Q1, BKS)
        got = ckks_decrypt(cr, s, Q1, DELTA)
        exp = [z[(j + r) % NS] for j in range(NS)]
        print(f"rotate by {r}: {[round(v.real,4) for v in got]}  "
              f"(expected {exp})")


# =====================================================================
# ------------------------  THE FULL PIPELINE  ------------------------
# =====================================================================
def main():
    self_tests()

    MSG = [3, 1, 0, 2]          # integers in [0, p_lwe)
    banner("PARAMETERS")
    print(f"  N (ring dim)      = {N}     slots = {NS}")
    print(f"  DELTA             = 2^{LOGDELTA}")
    print(f"  Q1 (fresh CKKS)   = 2^{Q1.bit_length()-1}")
    print(f"  Q0 (level 0 CKKS) = 2^{LOGQ0}")
    print(f"  n_lwe             = {n_lwe}")
    print(f"  q_lwe             = {q_lwe}   p_lwe = {p_lwe}")
    print(f"  CONST = Q0/(p*DELTA) = {CONST}")
    print(f"  message vector    = {MSG}")

    # ---------------- STEP 1: CKKS keys ------------------------------
    banner("STEP 1  -  CKKS key generation")
    s_ckks  = ckks_keygen()
    rotkeys = gen_rot_keys(s_ckks, Q1, BKS)
    print("  s_ckks  =", s_ckks)
    print("  rotation keys for r = 1..3 generated "
          f"({ks_len(Q1,BKS)} gadget digits each)")

    # ---------------- STEP 2: FHEW/TFHE key --------------------------
    banner("STEP 2  -  FHEW/TFHE (LWE) key generation")
    s_lwe = lwe_keygen(n_lwe)
    z     = embed_lwe_key(s_lwe)
    print("  s_lwe        =", s_lwe)
    print("  z(x) embed   =", z, "   <- s_lwe zero-padded to degree N-1")

    # ---------------- STEP 3: encrypt in CKKS ------------------------
    banner("STEP 3  -  CKKS encryption of the message vector")
    pt = ckks_encode(MSG, DELTA, Q1)
    ct = ckks_encrypt(pt, s_ckks, Q1)
    print("  plaintext poly (centered) :", [center(x, Q1) for x in pt])
    dec = ckks_decrypt(ct, s_ckks, Q1, DELTA)
    print("  decrypted slots           :", [round(v.real, 6) for v in dec])

    # ---------------- STEP 4: SlotsToCoeffs --------------------------
    banner("STEP 4  -  SlotsToCoeffs  (homomorphic linear transform)")
    diags = s2c_diagonals(CONST)
    ct_s2c = slots_to_coeffs(ct, diags, rotkeys, Q1, BKS, DELTA)
    print(f"  applied {NS} plaintext diagonals + {NS-1} rotations")
    print(f"  scale is now DELTA^2 = 2^{2*LOGDELTA}, modulus still 2^{Q1.bit_length()-1}")

    # ---------------- STEP 5: rescale --------------------------------
    banner("STEP 5  -  Rescale  Q1 -> Q0   (drop one level)")
    ct_lvl0, Qcur = rescale(ct_s2c, Q1, DELTA)
    assert Qcur == Q0
    ph = ckks_phase(ct_lvl0, s_ckks, Q0)
    expected = [(Q0 // p_lwe) * (MSG[k] if k < NS else 0) % Q0 for k in range(N)]
    print("  coefficients of the noisy plaintext, and target (Q0/p)*m_k :")
    print(f"    {'k':>2} {'phase coeff':>16} {'expected':>16} {'noise':>10} {'log2':>7}")
    for k in range(N):
        d = center(ph[k] - expected[k], Q0)
        lg = f"{math.log2(abs(d)):.1f}" if d else "-inf"
        print(f"    {k:>2} {center(ph[k],Q0):>16} {center(expected[k],Q0):>16} "
              f"{d:>10} {lg:>7}")
    print(f"  correctness threshold Q0/(2p) = {Q0//(2*p_lwe)} = 2^{(Q0//(2*p_lwe)).bit_length()-1}")

    # ---------------- STEP 6: key switch CKKS key -> embedded LWE key -
    banner("STEP 6  -  Key switch  s_ckks -> z  (the embedded LWE key)")
    swk    = gen_ksk(s_ckks, z, Q0, BKS)
    ct_z   = key_switch(ct_lvl0, swk, Q0, BKS)
    ph_z   = ckks_phase(ct_z, z, Q0)
    print(f"  switching key: {len(swk)} gadget digits at modulus 2^{LOGQ0}")
    print("  phase under z, and drift vs. the phase under s_ckks :")
    for k in range(N):
        print(f"    k={k}  {center(ph_z[k],Q0):>16}   drift = "
              f"{center(ph_z[k]-ph[k],Q0):>8}")

    # ---------------- STEP 7: LWE extraction -------------------------
    banner("STEP 7  -  Coefficient extraction  ->  LWE ciphertexts mod Q0")
    lwes_bigQ = [extract_lwe(ct_z, k, n_lwe, Q0) for k in range(NS)]
    for k, c in enumerate(lwes_bigQ):
        ph_l = lwe_phase(c, s_lwe, Q0)
        tgt  = (Q0 // p_lwe) * MSG[k]
        print(f"    slot {k}: phase = {center(ph_l,Q0):>16}  target = {center(tgt,Q0):>16}"
              f"  noise = {center(ph_l-tgt,Q0):>8}")
        assert ph_l == ph_z[k], "extraction must reproduce the k-th coefficient!"
    print("  OK: every extracted LWE phase equals the corresponding RLWE coefficient.")

    # ---------------- STEP 8: modulus switch -------------------------
    banner(f"STEP 8  -  Modulus switch  2^{LOGQ0} -> {q_lwe}")
    lwes = [lwe_modswitch(c, Q0, q_lwe) for c in lwes_bigQ]
    print(f"    {'slot':>4} {'a (mask)':>22} {'b':>6} {'phase':>7} "
          f"{'target':>7} {'noise':>6} {'msg':>4}")
    ok = True
    for k, c in enumerate(lwes):
        ph_l = lwe_phase(c, s_lwe, q_lwe)
        tgt  = (q_lwe // p_lwe) * MSG[k]
        m    = lwe_decrypt(c, s_lwe, q_lwe, p_lwe)
        ok  &= (m == MSG[k])
        print(f"    {k:>4} {str(c[0]):>22} {c[1]:>6} {ph_l:>7} {tgt:>7} "
              f"{center(ph_l-tgt,q_lwe):>6} {m:>4}")
    print(f"  decoding threshold q/(2p) = {q_lwe//(2*p_lwe)}")

    # =================================================================
    banner("VERIFICATION")
    got = [lwe_decrypt(c, s_lwe, q_lwe, p_lwe) for c in lwes]
    print(f"  [V1] recovered messages   : {got}")
    print(f"       original messages    : {MSG}")
    print(f"       MATCH                : {got == MSG}")

    noise = [abs(center(lwe_phase(c, s_lwe, q_lwe) - (q_lwe//p_lwe)*MSG[k], q_lwe))
             for k, c in enumerate(lwes)]
    print(f"  [V2] |noise| per LWE ct   : {noise}   (must be < {q_lwe//(2*p_lwe)})")

    # V3: they are genuine LWE ciphertexts -> additively homomorphic
    csum = lwe_add(lwes[0], lwes[1], q_lwe)
    print(f"  [V3] LWE(m0)+LWE(m1) decrypts to "
          f"{lwe_decrypt(csum,s_lwe,q_lwe,p_lwe)}  "
          f"(expected {(MSG[0]+MSG[1])%p_lwe})")

    # V4: wrong key must give garbage
    bad = [1 - b for b in s_lwe]
    print(f"  [V4] decrypting with a WRONG LWE key : "
          f"{[lwe_decrypt(c,bad,q_lwe,p_lwe) for c in lwes]}  (should be junk)")

    # V5: TFHE sign convention
    tfhe = [([(-x) % q_lwe for x in c[0]], c[1]) for c in lwes]
    ph_t = [(t[1] - sum(a*s for a, s in zip(t[0], s_lwe))) % q_lwe for t in tfhe]
    print(f"  [V5] same cts in TFHE convention (b - <a,s>) : phases {ph_t}")

    print("\n  RESULT:", "SCHEME SWITCHING CORRECT" if got == MSG else "FAILED")


if __name__ == "__main__":
    main()
