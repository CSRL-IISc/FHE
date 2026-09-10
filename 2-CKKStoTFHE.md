# CKKS → TFHE Scheme Switching

Step-by-step guide to converting a CKKS ciphertext into TFHE LWE ciphertexts so you can
evaluate non-polynomial functions (sign, ReLU, max, comparison, sigmoid, modular reduction)
via programmable bootstrapping.

Worked with toy parameters **N = 8, Q = 256** for illustration; production notes throughout.

> **Reference implementations of this idea:**
> - **CHIMERA** — Boura, Gama, Georgieva, Jetchev (2018/2020). Common-ring framework
>   bridging TFHE / BFV / HEAAN.
> - **PEGASUS** — Lu, Huang, Hong, Ma, Qian (IEEE S&P 2021). The CKKS ↔ FHEW/TFHE
>   pipeline described here, including repacking back to CKKS.
> - **CDKS repacking** — Chen, Dai, Kim, Song (2021). LWEs-to-RLWE packing.
> - Available in **OpenFHE** (`SchemeSwitching` / `EvalCKKStoFHEW`) and **Lattigo**.
uu
---
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

Companion document for `ckks_to_tfhe.py` — a from-scratch, dependency-free
replication of OpenFHE's `EvalCKKStoFHEW` pipeline with toy parameters
(`N = 8`) so that every intermediate value can be printed and checked by hand.

---

## 1. The one idea behind the whole thing

CKKS and TFHE build ciphertexts out of the same raw material: an inner product
plus noise.

| | ciphertext | decryption ("phase") |
|---|---|---|
| CKKS (RLWE) | `(c0, c1)` polynomials in `Z_Q[x]/(x^N+1)` | `c0 + c1·s mod Q` |
| TFHE (LWE) | `(a, b)`, `a ∈ Z_q^n`, `b ∈ Z_q` | `b + ⟨a, s⟩ mod q` |

Scheme switching is **not** a decrypt-and-re-encrypt. It is a sequence of
reformattings that never touch the plaintext. After the CKKS side, the message
is wrong in exactly four ways, and each one gets its own fix:

| What's wrong | Fix | Section |
|---|---|---|
| Message is in the wrong **place** — slots, not coefficients | homomorphic linear transform (SlotsToCoeffs) | §6 |
| Message is under the wrong **key** — `s_ckks`, not `s_lwe` | gadget key switching | §4, §8 |
| Message has the wrong **shape** — a polynomial, not integers | RLWE sample extraction | §9 |
| Message has the wrong **modulus/scale** — `Q0`, scale Δ | modulus switching | §10 |

Every function in the file serves one of these four, or is plumbing underneath
them.

---

## 2. Pipeline at a glance

```
CKKS ciphertext at Q1
   message lives in 4 slots
            │
            ▼   SlotsToCoeffs (4 pt-mults + 3 rotations), then rescale Q1 → Q0
   message now lives in coefficients 0..3
            │
            ▼   key switch  s_ckks → z(x)
   same plaintext, now under the embedded LWE key
            │
            ▼   coefficient extraction (pure bookkeeping, exact)
   4 LWE ciphertexts mod Q0
            │
            ▼   modulus switch Q0 → q
   4 standard TFHE ciphertexts mod q = 1024
```

### Mapping to the OpenFHE API

| OpenFHE call | What it really does | Section in the code |
|---|---|---|
| `EvalCKKStoFHEWSetup` | pick CKKS + LWE parameters, build the FHEW context | §0 |
| `EvalCKKStoFHEWKeyGen` | embed `s_lwe` into a ring element `z(x)`, build the RLWE key-switching key `s_ckks → z` | §8 |
| `EvalCKKStoFHEWPrecompute(scale)` | fold the constant into the SlotsToCoeffs matrix | §6 |
| `EvalCKKStoFHEW(ct, n)` | SlotsToCoeffs → rescale → key switch → extract → mod switch | §6, §3, §4, §9, §10 |

---

## 3. Parameters

| Symbol | Value | Role |
|---|---|---|
| `N` | 8 | RLWE ring dimension, `x^N + 1` |
| `M` | 16 | `2N`, order of the root of unity ζ |
| `NS` | 4 | CKKS slots = `N/2` |
| `DELTA` | 2³⁵ | CKKS scaling factor |
| `Q1` | 2⁸⁰ | fresh CKKS modulus (level 1) |
| `Q0` | 2⁴⁵ | CKKS modulus after one rescale (level 0) |
| `BKS` | 2⁵ | gadget / digit-decomposition base |
| `n_lwe` | 4 | LWE dimension (must be ≤ N) |
| `q_lwe` | 1024 | LWE modulus |
| `p_lwe` | 4 | LWE plaintext modulus → messages in {0,1,2,3} |
| `CONST` | 256 | `Q0 / (p·Δ)`, folded into the transform |

Two constraints are load-bearing, not cosmetic:

* **`Q0` and `Q1` are exact powers of `BKS`.** The balanced digit
  decomposition produces a carry out of the top digit; it only vanishes mod `Q`
  when `B^l = Q` exactly. Break this and key switching silently returns garbage.
* **`n_lwe ≤ N`.** The LWE key is zero-padded into a degree-`N` ring element.

There is no security at these sizes. `N = 8` is for readability only.

---

## 4. Section-by-section walkthrough

### §1 — Ring arithmetic in `Z_Q[x]/(x^N+1)`

`poly_mul` is schoolbook convolution with one twist: when `i + j ≥ N` the term
is **subtracted** into position `i + j − N`. That is `x^N = −1`. This sign flip
reappears in `auto_poly` (§5) and `extract_lwe` (§9), and it is the single most
common place to get scheme switching wrong.

`center(x, Q)` maps a residue into `(−Q/2, Q/2]`. Noise is a small *signed*
quantity, so a raw residue of `Q − 3` must be read as `−3`.

`round_div(x, d)` is exact-integer rounded division. At `Q1 = 2⁸⁰` a float
division would silently drop bits — this must never be `round(x/d)`.

### §2 — CKKS encoding / decoding (canonical embedding)

Slot `j` of a plaintext polynomial `p` is `p(ζ^{5^j})` with `ζ = e^{2πi/M}`.
The exponents `GAL = [1, 5, 9, 13]` are the powers of 5 mod 16.

* `ckks_decode(p)` literally evaluates `p` at those four points.
* `ckks_encode(z)` inverts it:

  ```
  m_k = (2/N) · Re( Σ_j  z_j · conj(ζ^{G_j · k}) )
  ```

The `Re(·)` is what forces the coefficients to be real, and it is exactly why
the whole switch only works on **real-valued** slots. Complex slots would need
coefficients that cannot exist.

The self-test at the top of the script confirms round-trip error ≈ 3e-11.

### §3 — CKKS keygen / encrypt / decrypt / rescale

`ckks_encrypt` sets `c1 = a` (uniform) and `c0 = pt − a·s + e`, so

```
phase(ct) = c0 + c1·s = pt + e
```

`rescale(ct, Q, Δ)` divides both components by Δ with rounding and drops the
modulus to `Q/Δ`. This is what stops the scale from squaring on every
multiplication. It is a purely modular operation — it stays correct even when
the plaintext is large enough to wrap around `Q`, which matters here because the
final coefficient is `(Q0/p)·m` and reaches `3Q0/4`.

### §4 — Gadget (digit-decomposition) key switching

The engine, used twice for two different purposes (rotations in §5, and the
CKKS→LWE switch in §8).

Goal: given a ciphertext under `s_from`, produce one under `s_to` with the same
phase. The switching key is

```
KSK_i = (b_i, a_i)   with   b_i + a_i·s_to = B^i · s_from + e_i
```

Decompose `c1` into balanced base-`B` digits `d_i` and take the digit-weighted
combination:

```
c0' = c0 + Σ d_i·b_i
c1' =      Σ d_i·a_i

c0' + c1'·s_to = c0 + Σ d_i(B^i·s_from + e_i) = c0 + c1·s_from + small
```

Digits are *balanced* (in `[−B/2, B/2)`), which roughly halves the noise
compared to plain non-negative digits.

### §5 — Galois automorphisms = slot rotations

Applying `x ↦ x^{5^r}` to the plaintext maps slot `j` to slot `j + r`, because

```
p(x^{5^r}) evaluated at ζ^{5^j}  =  p(ζ^{5^{j+r}})  =  z_{j+r}
```

Applied to a ciphertext it also transforms the key into `s(x^{5^r})`, so
`rotate` immediately key-switches back to `s`. A "rotation key" is nothing more
than a KSK from the twisted key to the original.

`auto_poly` implements the permutation: coefficient `k` moves to index
`k·g mod 2N`, with a sign flip if that index lands in `[N, 2N)`.

### §6 — SlotsToCoeffs (the crux)

Extraction can only reach **coefficients**, but CKKS puts the message in
**slots**. So the values have to be moved first.

We want a plaintext `p'` whose coefficient vector is
`(z0, z1, z2, z3, 0, 0, 0, 0)`. Its slot vector is then

```
slots(p') = A · z      with      A[j][k] = ζ^{5^j · k},   j,k = 0..3
```

So SlotsToCoeffs is just a 4×4 matrix applied homomorphically. It is evaluated
with the diagonal identity

```
A·z = Σ_r  diag_r(A) ⊙ rot(z, r)        diag_r(A)[j] = A[j][(j+r) mod NS]
```

which costs 4 plaintext multiplications and 3 rotations.

**The constant.** `CONST = Q0/(p·Δ) = 256` is folded into every diagonal. After
the transform and the rescale, the coefficient holds

```
Δ · CONST · m  =  (Q0/p) · m
```

This is what `EvalCKKStoFHEWPrecompute` does in OpenFHE, and it is what makes
the final modulus switch land *exactly* on the LWE encoding grid instead of at
some arbitrary scale. Without it the FHEW bootstrapper would see garbage.

### §7 — LWE / TFHE side

Plain LWE: `phase = b + ⟨a,s⟩ mod q`, decode by rounding `phase·p/q`. Also
`lwe_add`, used later as a correctness check.

### §8 — RLWE → LWE switching key

The FHEW/TFHE secret `s_lwe` (dimension `n`) is **embedded** into a ring
element:

```
z(x) = s_lwe[0] + s_lwe[1]·x + … + s_lwe[n-1]·x^{n-1}     (zero-padded to degree N-1)
```

Then a completely ordinary RLWE key switch from `s_ckks` to `z`. Once that is
done, every coefficient of the ciphertext **already is** an LWE sample under
`s_lwe`. No further cryptographic work is needed — only bookkeeping.

### §9 — Coefficient extraction (the bookkeeping)

Expand `(c1·z)_k` with the negacyclic rule and collect terms by `z[j]`:

```
(c1·z)_k = Σ_{j ≤ k} c1[k−j]·z[j]  −  Σ_{j > k} c1[N+k−j]·z[j]
```

so the LWE sample for coefficient `k` is

```
a[j] =  c1[k − j]        for j ≤ k
a[j] = −c1[N + k − j]    for j >  k
b    =  c0[k]
```

truncated to the first `n` entries, since `z[j] = 0` above that. This step is
**exact** — the assertion in Step 7 of the script checks that
`lwe_phase(extracted)` equals the RLWE coefficient bit for bit.

### §10 — Modulus switching

Scale everything by `q/Q0` and round:

```
a'[j] = round(a[j] · q / Q0)   mod q
b'    = round(b    · q / Q0)   mod q
```

Noise shrinks by the same factor (here from ~2¹⁸ down to below 1), and the
rounding introduces a fresh error of order `‖s‖₁/2`. The output is
`(q/p)·m + tiny` — exactly what a FHEW bootstrapper expects.

---

## 5. Reading the program output

Three numbers are worth watching whenever you change parameters.

**Step 5 — coefficients after rescale.** Prints `k = 0..7` against the target
`(Q0/p)·m_k`. Coefficients 4–7 must be pure noise near zero. If they are not,
SlotsToCoeffs is wrong.

```
   k      phase coeff         expected      noise    log2
   0   -8796093261073   -8796093022208    -238865    17.9
   1    8796092928075    8796093022208     -94133    16.5
   2           -55758                0     -55758    15.8
   3   17592185760758   17592186044416    -283658    18.1
   4           252681                0     252681    17.9
   ...
   correctness threshold Q0/(2p) = 2^42
```

**Step 6 — drift.** Noise added purely by the key switch to `z`. It should be
tiny (hundreds) compared to the transform noise (~2¹⁸). If drift dominates,
the gadget base `BKS` is too large.

**Step 8 — final noise.** Measured against the decoding threshold
`q/(2p) = 128`. The reference run ends at `|noise| = 1`, i.e. about 7 bits of
margin. If that number creeps toward 128, raise Δ or lower `BKS`.

```
   slot               a (mask)      b   phase  target  noise  msg
      0   [610, 579, 680, 472]     62     769     768      1    3
      1   [703, 610, 579, 680]    436     257     256      1    1
      2   [810, 703, 610, 579]    157       1       0      1    0
      3   [374, 810, 703, 610]    438     513     512      1    2
```

### Where the noise actually comes from

The dominant term, by a wide margin, is the **rotation key-switch error
multiplied by the plaintext diagonal magnitude** before the rescale divides it
back down:

```
final_noise  ≈  e_ks · (Δ · CONST) · N · NS / Δ  =  e_ks · CONST · N · NS
```

Relative to the threshold `Q0/(2p)`, this is roughly `64·e_ks / Δ`, so Δ must be
comfortably larger than `64·e_ks`. That is why Δ = 2³⁵ rather than something
smaller.

---

## 6. How to verify the switch is correct

Five independent checks, in increasing strength. Passing only the last one is
not enough — a decrypt can be accidentally correct.

| # | Check | Why it matters |
|---|---|---|
| **V0** | Coefficients after rescale ≈ `(Q0/p)·m_k`, and ≈ 0 for `k ≥ 4` | isolates SlotsToCoeffs from everything downstream |
| **V1** | `assert lwe_phase(extract(ct,k)) == rlwe_phase(ct)[k]` | exact, not approximate — catches index and sign errors immediately |
| **V2** | Recovered messages equal the originals | end-to-end correctness |
| **V3** | `|phase − (q/p)·m| ≪ q/(2p)` | a passing decrypt with noise near the threshold means you got lucky, not correct — **always print this number** |
| **V4** | `LWE(m0) + LWE(m1)` decrypts to `m0+m1 mod p`; a wrong key gives junk | proves these are genuine LWE ciphertexts, rules out "accidentally correct because noise happened to be zero" |

The script also prints the same ciphertexts in the TFHE sign convention
(`b − ⟨a,s⟩`) — the converter is a one-liner that negates `a`.

---

## 7. Common failure modes

| Symptom | Likely cause |
|---|---|
| Rotations return the wrong slots | direction: `rot(z,r)[j] = z[j+r]`, not `z[j−r]` |
| Extraction assertion fails | sign flip from `x^N = −1` missed in `extract_lwe` or `auto_poly` |
| Key switching returns noise | `B^l ≠ Q` — the top-digit carry does not vanish mod Q |
| Decryption off by a constant factor | `CONST` computed as something other than `Q0/(p·Δ)` |
| Works for small messages, fails for `m = p−1` | plaintext wrapping — check that rescale uses the *centered* representative |
| Phase looks negated | `b + ⟨a,s⟩` vs `b − ⟨a,s⟩` convention mismatch |
| Precision loss at large moduli | a float `round(x/d)` somewhere instead of integer `round_div` |

---

## 8. QUESTIONS 

1. **S2C & NTT Question**For the S2C evaluation, we have to compute the homomorphic matrix multiplication. Since this requires many ciphertext-plaintext polynomial multiplications, should I build a dedicated NTT (Number Theoretic Transform) pipeline for this? Or given the structure of the S2C matrix, can we compute it more efficiently using a systolic array architecture to minimize memory reads? 
2. **Galois Key Streaming** S2C requires many slot rotations, which means we need to pull in massive Galois evaluation keys. Since these won't fit entirely in on-chip SRAM/BRAM, how do you recommend architecting the memory interface to stream these from off-chip DDR without stalling the polynomial multiplication pipelines?
3. **RNS** CKKS relies heavily on Residue Number System (RNS) limbs to handle the large $Q$. At what exact stage in the hardware pipeline should we break out of RNS and drop down to the single-limb TFHE representation to minimize routing congestion?
4. How should I handle the CKKS scaling factor in hardware? Should I build a dedicated fixed-point rounding unit right before the modulus switch, or can we just truncate the lower bits to save area?
5. How deep we have to understand the code or mathematics part?
6. What will be the broad architecture design, Should we build a direct pipeline where every step—S2C, Mod-Switch, and Blind Rotation—gets its own dedicated hardware? Or should we use a central controller that reuses the same polynomial math units for everything?

