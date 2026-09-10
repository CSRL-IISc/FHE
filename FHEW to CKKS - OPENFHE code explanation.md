# OpenFHE FHEW → CKKS Scheme Switching — Function Call Tree and Code Explanation

## 1. Scope

This document covers **only the FHEW → CKKS direction** in `ckksrns-schemeswitching.cpp`.

The central function is:

```cpp
SWITCHCKKSRNS::EvalFHEWtoCKKS(...)
```

The key idea is that OpenFHE does **not decrypt the FHEW/LWE ciphertexts and then encrypt the plaintexts into CKKS**. Instead, it reconstructs the LWE decryption expression homomorphically inside CKKS.

For an LWE ciphertext `(a, b)`:

```text
b = a · s + Δm + e   (mod q)
```

so the decryption quantity is:

```text
b - a · s = Δm + e   (mod q)
```

OpenFHE represents `a` and `b` as CKKS plaintext data, encrypts the FHEW secret key `s` as a CKKS ciphertext, homomorphically computes `a · s`, and then evaluates the remainder of the LWE-to-CKKS conversion.

---

# 2. High-Level FHEW → CKKS Tree

```text
EvalFHEWtoCKKS()
|
+-- determine number of LWE ciphertexts / CKKS slots
|
+-- determine FHEW modulus-dependent parameters
|   |
|   +-- choose K
|   +-- choose Chebyshev coefficients
|
+-- Step 1: build A and b from LWE ciphertexts
|   |
|   +-- GetA()
|   +-- GetB()
|   +-- convert to complex<double>
|   +-- prescale by 1/(q*K)
|
+-- Step 2: compute A*s homomorphically
|   |
|   +-- EvalPartialHomDecryption()
|       |
|       +-- pad A to power-of-two columns
|       |
|       +-- EvalLTRectPrecomputeSwitch()
|       |       |
|       |       +-- extract matrix diagonals
|       |       +-- arrange diagonals for BSGS
|       |       +-- scale diagonals
|       |
|       +-- EvalLTRectWithPrecomputeSwitch()
|               |
|               +-- EvalFastRotationPrecompute()
|               +-- EvalFastRotationExt()
|               +-- KeySwitchExt()
|               +-- EvalMultExt()
|               +-- EvalAddExtInPlace()
|               +-- giant-step rotations
|               +-- KeySwitchDown()
|
+-- Step 3: compute B - A*s
|   |
|   +-- MakeCKKSPackedPlaintext(b)
|   +-- EvalNegate(AdotS)
|   +-- EvalAdd(BPlain, -AdotS)
|
+-- Step 4: approximate modular reduction
|   |
|   +-- EvalChebyshevSeries()
|   +-- iterative sine refinement
|       +-- EvalMult()
|       +-- EvalAddInPlace()
|       +-- EvalSubInPlace()
|
+-- Step 5: post-scale and post-bias
|   |
|   +-- create postScale plaintext
|   +-- EvalMult()
|   +-- create postBias plaintext
|   +-- EvalAddInPlace()
|
+-- restore sparse encoding if required
|
+-- optional ModReduce()
|
+v
CKKS ciphertext containing the converted FHEW result
```

The actual source explicitly describes the major stages as: forming `A` and `b`, homomorphically computing `A*skLWE`, obtaining `B - A*s`, performing modular reduction with a sine approximation, and applying the final scaling/bias. fileciteturn13file0L32-L62 fileciteturn15file0L10-L57

---

# 3. What Is Being Converted?

Suppose we have several LWE/FHEW ciphertexts:

```text
c0 = (a0, b0)
c1 = (a1, b1)
c2 = (a2, b2)
...
```

Each one satisfies approximately:

```text
bi = ai · s + Δmi + ei   (mod q)
```

The goal is to produce one CKKS ciphertext containing the messages:

```text
[m0, m1, m2, ...]
```

The crucial operation is therefore:

```text
bi - ai · s
```

for every input ciphertext.

Because `s` must remain secret, OpenFHE does not compute this in plaintext. It encodes/encrypts the secret key into CKKS and evaluates the dot products homomorphically.

---

# 4. `EvalSchemeSwitchingKeyGen()` — Preparing the Key for FHEW → CKKS

Before `EvalFHEWtoCKKS()` can run, a special switching key is generated.

The relevant source first obtains the LWE secret key and rounds its entries to a representation suitable for CKKS. fileciteturn16file0L10-L34

Conceptually:

```text
LWE secret key s
      |
      | convert coefficients
      v
CKKS plaintext representation of s
      |
      | encrypt with CKKS public key
      v
m_FHEWtoCKKSswk
```

The code does:

```cpp
uint32_t n_po2 = 1 << static_cast<uint32_t>(std::ceil(std::log2(n)));
```

so the LWE dimension is padded to the nearest power of two.

Then:

```cpp
skLWEDouble[i] = std::complex<double>(tmp == neg ? -1.0 : tmp, 0);
```

This maps the LWE representation into the real-valued CKKS representation. In particular, a coefficient equal to `q-1` is interpreted as `-1`.

The secret key is then encoded as a CKKS packed plaintext and encrypted:

```cpp
m_FHEWtoCKKSswk = ccCKKS->Encrypt(publicKey, skLWEPlainswk);
```

The source comment explicitly describes this object as the **CKKS encryption of the FHEW secret key**. fileciteturn16file0L26-L44

This is the central trick that makes homomorphic partial LWE decryption possible.

---

# 5. The LWE Equation Being Evaluated

The source eventually wants the CKKS ciphertext corresponding to:

```text
B - A*s
```

The mathematical mapping is:

```text
b = a · s + Δm + e   (mod q)
```

therefore:

```text
b - a · s = Δm + e   (mod q)
```

For multiple LWE ciphertexts, stack their `a` vectors into a matrix:

```text
A =
[ a0
  a1
  a2
  ... ]
```

and stack their `b` values into:

```text
B = [b0, b1, b2, ...]
```

Then:

```text
A*s =
[ a0 · s
  a1 · s
  a2 · s
  ... ]
```

and:

```text
B - A*s =
[ b0 - a0 · s
  b1 - a1 · s
  b2 - a2 · s
  ... ]
```

This is exactly what the source implements in two stages: compute `AdotS`, then compute `B - AdotS`. fileciteturn13file0L42-L62

---

# 6. `EvalFHEWtoCKKS()`

Signature:

```cpp
Ciphertext<DCRTPoly> SWITCHCKKSRNS::EvalFHEWtoCKKS(
    std::vector<std::shared_ptr<LWECiphertextImpl>>& LWECiphertexts,
    uint32_t numCtxts,
    uint32_t numSlots,
    uint32_t p,
    double pmin,
    double pmax,
    uint32_t dim1) const
```

The function first rejects an empty input vector and determines how many LWE ciphertexts will be packed into the CKKS result. It caps this number by the available CKKS slots. fileciteturn12file4L279-L300

Conceptually:

```text
input LWE ciphertexts
        |
        v
choose numValues
        |
        v
pack numValues results into CKKS slots
```

---

# 7. Selecting `K` and the Polynomial Approximation

The code selects a parameter `K` based on the LWE dimension:

```cpp
if (n == 32) {
    K = 16.0;
    coefficientsFHEW.assign(g_coefficientsFHEW16);
}
else {
    K = 128.0;
    if (p <= 4)
        coefficientsFHEW.assign(g_coefficientsFHEW128_8);
    else
        coefficientsFHEW.assign(g_coefficientsFHEW128_9);
}
```

The source comments associate `K = 128` with the desired failure probability and distinguish the polynomial used for bit messages from the polynomial used for larger plaintext moduli. fileciteturn13file0L15-L30

The practical purpose is to normalize the value before the homomorphic modular-reduction approximation.

---

# 8. Step 1 — Construct `A` and `b`

The source creates:

```cpp
std::vector<std::vector<std::complex<double>>> A(numValues);
std::vector<std::complex<double>> b(b_size);
```

`A[i]` stores the `a` vector of the `i`-th LWE ciphertext. `b[i]` stores its `b` value after prescaling. fileciteturn13file0L32-L49

The code computes:

```cpp
const double prescale =
    (1.0 / LWECiphertexts[0]->GetModulus().ConvertToDouble()) / K;
```

Therefore:

```text
prescale = 1 / (q*K)
```

where `q` is the LWE modulus.

Then:

```cpp
A[i][j] = std::complex<double>(a[j].ConvertToDouble(), 0);
```

and:

```cpp
b[i] = std::complex<double>(prescale * B_i, 0);
```

So the source is preparing the numerical data so that the subsequent CKKS computation operates on a normalized representation. fileciteturn15file1L97-L114

---

# 9. Why Is `A` a Matrix Here?

This `A` is different from the `A` used in the earlier CKKS → FHEW discussion.

Here:

```text
A[i][j] = a_j of the i-th LWE ciphertext
```

Therefore, if there are `r` LWE ciphertexts and LWE dimension `n`:

```text
A has approximately r rows × n columns
```

For example, with three LWE ciphertexts:

```text
A =
[ a00 a01 a02 a03
  a10 a11 a12 a13
  a20 a21 a22 a23 ]
```

Then:

```text
A · s =
[ a00*s0 + a01*s1 + a02*s2 + a03*s3
  a10*s0 + a11*s1 + a12*s2 + a13*s3
  a20*s0 + a21*s1 + a22*s2 + a23*s3 ]
```

That vector is exactly the collection of inner products required by LWE decryption.

---

# 10. Step 2 — `EvalPartialHomDecryption()`

The source calls:

```cpp
auto AdotS = EvalPartialHomDecryption(
    *ccCKKS, A, m_FHEWtoCKKSswk, dim1, prescale, 0);
```

The source comment is explicit:

```text
Step 2. Perform the homomorphic linear transformation of A*skLWE
```

This is the heart of the FHEW → CKKS scheme switch. fileciteturn13file0L51-L54

Conceptually:

```text
A
 |
 | plaintext matrix
 v
A*s
 |
 | evaluated homomorphically in CKKS
 v
AdotS ciphertext
```

The switching key `m_FHEWtoCKKSswk` supplies the encrypted secret key `s`.

---

# 11. `EvalPartialHomDecryption()` — Why Padding Is Needed

The function first copies `A` and computes the next power of two for its number of columns:

```cpp
size_t cols_po2 =
    1 << static_cast<uint32_t>(std::ceil(std::log2(A[0].size())));
```

If necessary, each row is resized to that power-of-two length. fileciteturn11file0L30-L45

Conceptually:

```text
original LWE dimension n
        |
        v
next power of two
        |
        v
pad each a vector with zeros
```

This is required because the rectangular linear-transform machinery is designed around power-of-two dimensions.

The function then calls:

```cpp
auto Apre = EvalLTRectPrecomputeSwitch(Acopy, dim1, scale);
```

followed by:

```cpp
return EvalLTRectWithPrecomputeSwitch(
    cc, Apre, ct,
    (Acopy.size() < A[0].size()), dim1, L);
```

The source comment states that the result is repeated every `Acopy.size()` slots. fileciteturn11file0L34-L45

Thus:

```text
EvalPartialHomDecryption()
    = prepare A for the rectangular linear transform
      and evaluate A*s homomorphically
```

---

# 12. `EvalLTRectPrecomputeSwitch()`

This function prepares the rectangular matrix for the linear-transform engine.

For power-of-two matrix dimensions it computes:

```cpp
const uint32_t n = std::min(A.size(), A[0].size());
```

and creates a vector of diagonal representations:

```cpp
std::vector<std::vector<std::complex<double>>> diags(n);
```

For the case where the matrix has at least as many rows as columns, it computes BSGS parameters:

```cpp
bStep = ...
gStep = ceil(n / bStep)
```

and then extracts the shifted diagonals with:

```cpp
ExtractShiftedDiagonal(A_slices[k], bStep * j + i)
```

The extracted diagonal is multiplied by `scale` and stored in `diags[...]`. fileciteturn13file5L369-L404

So conceptually:

```text
matrix A
  |
  +-- shifted diagonal D0
  +-- shifted diagonal D1
  +-- shifted diagonal D2
  +-- ...
  |
  v
precomputed diagonal representation
```

This is the same diagonal linear-transform technique used by CKKS matrix-vector multiplication.

---

# 13. Why Does the Linear Transform Compute `A*s`?

Suppose:

```text
A =
[ a00 a01 a02 a03
  a10 a11 a12 a13
  a20 a21 a22 a23 ]
```

and the encrypted secret key is:

```text
s = [s0, s1, s2, s3]
```

Then the desired result is:

```text
A*s =
[ a00*s0 + a01*s1 + a02*s2 + a03*s3
  a10*s0 + a11*s1 + a12*s2 + a13*s3
  a20*s0 + a21*s1 + a22*s2 + a23*s3 ]
```

The diagonal method rewrites the same computation as a sequence of:

```text
rotate s
multiply by a plaintext diagonal of A
add the results
```

Conceptually:

```text
A*s
  = D0 ⊙ Rot0(s)
  + D1 ⊙ Rot1(s)
  + D2 ⊙ Rot2(s)
  + ...
```

Here the `Dk` values come from `ExtractShiftedDiagonal()`.

The important difference from ordinary matrix multiplication is that `s` is not available as plaintext. In OpenFHE, `s` is already encrypted by `m_FHEWtoCKKSswk`, so all of these operations are performed homomorphically.

---

# 14. `EvalLTRectWithPrecomputeSwitch()`

This function actually evaluates the precomputed linear transform.

It first determines the BSGS parameters:

```cpp
uint32_t bStep = (dim1 == 0) ? getRatioBSGSLT(n) : dim1;
uint32_t gStep = std::ceil(static_cast<double>(n) / bStep);
```

Then it computes the hoisted-rotation preparation:

```cpp
auto digits = cc.EvalFastRotationPrecompute(ct);
```

and creates:

```cpp
std::vector<Ciphertext<DCRTPoly>> fastRotation(bStep - 1);
```

The source then generates the baby-step rotations:

```cpp
for (uint32_t j = 1; j < bStep; ++j)
    fastRotation[j - 1] = cc.EvalFastRotationExt(ct, j, digits, true);
```

This is the concrete implementation of the repeated rotations required by the diagonal linear transform. fileciteturn14file0L15-L30 fileciteturn14file0L69-L72

---

# 15. `EvalFastRotationPrecompute()`

The source says:

```cpp
auto digits = cc.EvalFastRotationPrecompute(ct);
```

Conceptually this prepares common NTT-related information for several subsequent automorphisms/rotations.

It is an optimization rather than an additional mathematical term.

Think:

```text
one expensive common preparation
          |
          +--> Rot1(s)
          +--> Rot2(s)
          +--> Rot3(s)
          +--> ...
```

The source comment identifies this as computing the NTTs used for hoisted automorphisms. fileciteturn14file0L27-L30

---

# 16. `EvalFastRotationExt()`

The source creates the baby-step rotations with:

```cpp
cc.EvalFastRotationExt(ct, j, digits, true)
```

Conceptually:

```text
EvalFastRotationExt(ct, j, ...)
              |
              v
           Rot_j(ct)
```

In this FHEW → CKKS path, `ct` contains the encrypted representation of the FHEW secret key.

Therefore these rotations are effectively rotations of encrypted `s`.

---

# 17. `KeySwitchExt()`

The core accumulation begins with:

```cpp
auto inner = FHECKKSRNS::EvalMultExt(
    cc.KeySwitchExt(ctxt, true),
    A[bStep * j]);
```

The `KeySwitchExt()` call prepares the ciphertext representation needed by the extended-RNS multiplication path.

It is not a separate mathematical term in `A*s`.

The mathematical operation is still:

```text
Dk ⊙ Rotk(s)
```

---

# 18. `EvalMultExt()`

The source computes terms such as:

```cpp
FHECKKSRNS::EvalMultExt(
    fastRotation[i - 1],
    A[bStep * j + i])
```

Conceptually this is:

```text
Dk ⊙ Rotk(s)
```

where:

```text
Dk = corresponding plaintext diagonal of A
```

and:

```text
Rotk(s) = rotated encrypted secret key
```

This is the homomorphic matrix-vector multiplication term.

---

# 19. `EvalAddExtInPlace()`

The code accumulates the diagonal products using:

```cpp
FHECKKSRNS::EvalAddExtInPlace(
    inner,
    FHECKKSRNS::EvalMultExt(...))
```

Conceptually:

```text
inner = inner + Dk ⊙ Rotk(s)
```

After all relevant terms are included:

```text
inner ≈ A*s
```

The OpenFHE source performs exactly this baby-step accumulation inside the giant-step loop. fileciteturn14file1L92-L115

---

# 20. BSGS: Why `A[bStep*j+i]` Appears

The linear transform is organized using:

```text
k = j*bStep + i
```

where:

```text
j = giant-step index
i = baby-step index
```

Therefore the source uses:

```cpp
A[bStep * j + i]
```

which identifies the diagonal associated with rotation index `k`.

The important conceptual expression is:

```text
A*s
  = Σk Dk ⊙ Rotk(s)
```

while the implementation groups those terms into BSGS blocks.

---

# 21. Giant-Step Handling

After forming the baby-step sum for a giant step, the source performs:

```cpp
inner = cc.KeySwitchDown(inner);
```

then determines the automorphism corresponding to the giant-step rotation:

```cpp
uint32_t autoIndex =
    FindAutomorphismIndex2nComplex(bStep * j, M);
```

and prepares the corresponding automorphism map:

```cpp
PrecomputeAutoMap(N, autoIndex, &map);
```

The source then transforms the first element and applies a fast rotation for the giant-step amount. fileciteturn14file1L99-L115

Conceptually:

```text
baby-step products
       |
       v
sum for this block
       |
       v
giant-step rotation
       |
       v
add to overall result
```

Thus BSGS reduces the cost of evaluating the complete diagonal sum.

---

# 22. `KeySwitchDown()`

The source finishes the linear-transform accumulation with:

```cpp
result = cc.KeySwitchDown(result);
result->GetElements()[0] += first;
```

This returns the extended-RNS result to the ordinary representation and restores the part accumulated separately in `first`. fileciteturn14file1L113-L120

This is implementation machinery around the linear transform; mathematically the output is still the ciphertext of:

```text
A*s
```

---

# 23. Return to `EvalFHEWtoCKKS()` — Build `BPlain`

Once `AdotS` has been generated, the source encodes the vector `b` as a CKKS plaintext:

```cpp
Plaintext BPlain = ccCKKS->MakeCKKSPackedPlaintext(
    b,
    AdotS->GetNoiseScaleDeg(),
    AdotS->GetLevel(),
    nullptr,
    N / 2);
```

This is important because `b` is not encrypted yet. It is known data extracted from the LWE ciphertexts.

So now:

```text
BPlain = plaintext containing prescaled b values
AdotS = ciphertext containing prescaled A*s values
```

---

# 24. Computing `B - A*s`

The source does:

```cpp
auto BminusAdotS = ccCKKS->EvalAdd(
    ccCKKS->EvalNegate(AdotS),
    BPlain);
```

Mathematically:

```text
BminusAdotS = B - AdotS
```

which corresponds to:

```text
b - a · s
```

for each LWE ciphertext.

Using the LWE equation:

```text
b = a · s + Δm + e   (mod q)
```

we therefore get:

```text
b - a · s = Δm + e   (mod q)
```

This is the central correctness equation for the FHEW → CKKS switch.

The source explicitly labels this as **Step 3: Get the ciphertext of B - A*s**. fileciteturn13file0L56-L62

---

# 25. Why Isn't `B - A*s` Already the Final CKKS Plaintext?

It is not yet in the desired CKKS representation.

The value is still a quantity represented modulo the LWE modulus and normalized by the earlier scaling.

The source therefore performs a homomorphic modular-reduction/approximation stage before applying the final output scaling.

The source calls this:

```text
Step 4. Do the modulus reduction: homomorphically evaluate modular function.
```

and implements it using a sine approximation through a Chebyshev polynomial followed by iterative refinement. fileciteturn13file0L56-L62 fileciteturn14file3L216-L234

---

# 26. `EvalChebyshevSeries()`

The source performs:

```cpp
auto BminusAdotS3 = ccCKKS->EvalChebyshevSeries(
    BminusAdotS,
    coefficientsFHEW,
    -1.0,
    1.0);
```

The earlier `prescale` included the factor `1/K`, so the input is normalized before this approximation.

The polynomial coefficients selected earlier approximate the required sine-based modular operation.

Conceptually:

```text
B - A*s
     |
     | normalized
     v
approximately wrapped/modular value
```

This is not ordinary plaintext arithmetic. `EvalChebyshevSeries()` evaluates the polynomial **homomorphically** on the ciphertext.

---

# 27. Iterative Sine Refinement

After the Chebyshev approximation, the code executes three iterations:

```cpp
const int32_t BT_ITER = 3;
for (int32_t j = 1; j <= BT_ITER; ++j) {
    BminusAdotS3 = ccCKKS->EvalMult(
        BminusAdotS3,
        BminusAdotS3);

    ccCKKS->EvalAddInPlace(
        BminusAdotS3,
        BminusAdotS3);

    double scalar = 1.0 / std::pow(
        (2.0 * M_PI),
        std::pow(2.0, j - BT_ITER));

    ccCKKS->EvalSubInPlace(BminusAdotS3, scalar);
}
```

Mathematically, if the current approximation is `x`, the update is of the form:

```text
x -> 2*x^2 - scalar
```

with the scalar chosen according to the iteration.

This is the repeated nonlinear refinement used by the implementation to obtain the desired periodic/modular behavior.

The source comments explicitly identify this stage as part of the sine approximation. fileciteturn14file3L224-L234

---

# 28. Why Use a Polynomial at All?

The CKKS evaluator can evaluate additions and multiplications, but it cannot directly execute an arbitrary operation such as:

```text
x mod q
```

as a native operation on encrypted data.

Therefore OpenFHE approximates the needed periodic function with a polynomial. Polynomial evaluation can be implemented with homomorphic additions and multiplications.

So the chain is:

```text
desired modular reduction
        |
        v
approximate with polynomial
        |
        v
EvalChebyshevSeries()
        |
        v
homomorphic multiplications/additions
```

---

# 29. Post-Scaling

After the modular approximation, the code determines:

```cpp
double postScale =
    (p >= 1 && p <= 4) ? (2.0 * M_PI) : static_cast<double>(p);
```

and, if a nonzero `pmin` is supplied:

```cpp
postScale *= (pmax - pmin) / 4.0;
postBias = (pmax - pmin) / 4.0;
```

The source comments explain that this accounts for the plaintext modulus and, when needed, the desired output interval. fileciteturn15file0L10-L26

Conceptually:

```text
normalized recovered value
        |
        | × postScale
        v
correct CKKS numerical scale
```

---

# 30. Post-Scale Plaintext Multiplication

The code constructs a plaintext containing `postScale` in the first
`numValues` positions:

```cpp
auto postScalePlain =
    ccCKKS->MakeCKKSPackedPlaintext(...);
```

and then:

```cpp
auto BminusAdotSres =
    ccCKKS->EvalMult(BminusAdotS3, postScalePlain);
```

Mathematically:

```text
BminusAdotSres
    = postScale × BminusAdotS3
```

This converts the normalized approximation into the desired CKKS output range. fileciteturn15file0L24-L37

---

# 31. Post-Bias

If the output interval requires a bias, the code creates:

```cpp
postBiasPlain
```

and performs:

```cpp
ccCKKS->EvalAddInPlace(
    BminusAdotSres,
    postBiasPlain);
```

Mathematically:

```text
output = scaled_value + postBias
```

This is what converts the normalized interval to the requested `[pmin, pmax]`-style output encoding. fileciteturn15file0L38-L42

---

# 32. Restoring Sparse Encoding

The code checks:

```cpp
if (isSparse) {
    for (uint32_t j = 1; j < N / (2 * slots); j <<= 1) {
        auto temp = ccCKKS->EvalAtIndex(
            BminusAdotSres,
            j * slots);
        ccCKKS->EvalAddInPlace(BminusAdotSres, temp);
    }
    BminusAdotSres->SetSlots(slots);
}
```

Conceptually:

```text
partially packed result
        |
        v
combine repeated blocks
        |
        v
restore sparse CKKS layout
```

The source calls this “Go back to the sparse encoding if needed.” fileciteturn15file0L44-L51

---

# 33. Final Modulus Reduction

For `FIXEDMANUAL` scaling, the source finally performs:

```cpp
ccCKKS->ModReduceInPlace(BminusAdotSres);
```

This consumes a modulus level and returns the ciphertext at the expected level/scaling state. fileciteturn15file0L53-L57

---

# 34. Complete Mathematical Flow

The entire FHEW → CKKS conversion can be summarized as:

```text
LWE ciphertexts:

(ai, bi)

where

bi = ai · s + Δmi + ei   (mod q)

        |
        | encode ai as matrix A
        | encode bi as vector B
        v

A = [a0; a1; ...]
B = [b0, b1, ...]

        |
        | CKKS encrypt s
        v

Enc(s)

        |
        | homomorphic linear transform
        v

Enc(A*s)

        |
        | plaintext B - ciphertext A*s
        v

Enc(B - A*s)

        |
        | LWE equation
        v

Enc(Δm + e)

        |
        | polynomial / sine-based modular reduction
        v

normalized message representation

        |
        | postScale + postBias
        v

CKKS ciphertext containing the FHEW results
```

---

# 35. Complete Function Tree With Mathematical Meaning

```text
EvalFHEWtoCKKS()
|
+-- Step 1: build A and B
|      |
|      +-- LWE.GetA()
|      |      -> a_i
|      |
|      +-- LWE.GetB()
|             -> b_i
|
+-- Step 2: EvalPartialHomDecryption()
|      |
|      +-- pad A to power-of-two columns
|      |
|      +-- EvalLTRectPrecomputeSwitch()
|      |      |
|      |      +-- ExtractShiftedDiagonal()
|      |      |      -> D_k
|      |      |
|      |      +-- scale D_k
|      |      +-- arrange diagonals for BSGS
|      |
|      +-- EvalLTRectWithPrecomputeSwitch()
|             |
|             +-- EvalFastRotationPrecompute()
|             |      -> rotation preparation
|             |
|             +-- EvalFastRotationExt()
|             |      -> Rot_k(Enc(s))
|             |
|             +-- KeySwitchExt()
|             |      -> representation management
|             |
|             +-- EvalMultExt()
|             |      -> D_k ⊙ Rot_k(Enc(s))
|             |
|             +-- EvalAddExtInPlace()
|             |      -> accumulate Σ D_k ⊙ Rot_k(s)
|             |
|             +-- giant-step automorphism / rotation
|             |
|             +-- KeySwitchDown()
|                    -> normal CKKS representation
|
|      -> AdotS = Enc(A*s)
|
+-- Step 3: B - A*s
|      |
|      +-- MakeCKKSPackedPlaintext(B)
|      +-- EvalNegate(AdotS)
|      +-- EvalAdd()
|             -> Enc(B - A*s)
|
+-- Step 4: modular reduction / sine approximation
|      |
|      +-- EvalChebyshevSeries()
|      |      -> polynomial approximation
|      |
|      +-- EvalMult()
|      +-- EvalAddInPlace()
|      +-- EvalSubInPlace()
|             -> iterative refinement
|
+-- Step 5: output scaling
|      |
|      +-- EvalMult(postScalePlain)
|      +-- EvalAddInPlace(postBiasPlain)
|
+-- sparse-encoding restoration if needed
|
+-- optional ModReduce()
|
+v
CKKS ciphertext
```

---

# 36. What `m_FHEWtoCKKSswk` Really Is

The most important object to understand is:

```cpp
m_FHEWtoCKKSswk
```

It is a **CKKS ciphertext encrypting the FHEW/LWE secret key**.

So the data flow is:

```text
FHEW secret key s
       |
       | encode as CKKS plaintext
       v
CKKS plaintext s
       |
       | Encrypt(publicKey, ...)
       v
m_FHEWtoCKKSswk = Enc(s)
```

Then:

```text
A plaintext
     ×
Enc(s)
     |
     v
Enc(A*s)
```

This is why the process is called **partial homomorphic decryption**: the decryption formula is evaluated homomorphically, but the final result is not produced by explicitly decrypting the LWE ciphertext in ordinary software.

The source explicitly creates the CKKS encryption of the FHEW secret key during switching-key generation. fileciteturn16file0L26-L44

---

# 37. The Most Important Distinction: `A` vs `Dk`

There are three different objects to keep separate:

```text
A
=
original matrix formed from the LWE a-vectors
```

then:

```text
A
 |
 | ExtractShiftedDiagonal()
 v
D0, D1, D2, ...
```

and finally:

```text
Dk
 |
 | plaintext representation
 v
precomputed diagonal data
```

Thus:

```text
LWE ciphertexts
      |
      v
A = matrix of a-vectors
      |
      v
Dk = diagonals of A
      |
      v
Σ Dk ⊙ Rotk(Enc(s))
      |
      v
Enc(A*s)
```

This is the cleanest way to connect the source code to the diagonal linear-transform formulation.

---

# 38. Why This Is Not Ordinary Decryption

Ordinary LWE decryption would be:

```text
receive (a,b)
      |
      v
compute b - a*s using secret key s
      |
      v
round / decode
```

But in scheme switching, the goal is to keep the result encrypted in CKKS.

OpenFHE therefore performs:

```text
receive encrypted LWE ciphertexts
          |
          v
extract a and b as numerical data
          |
          +--------------------------+
          |                          |
          v                          v
       A plaintext              B plaintext
          |                          |
          | × Enc(s)                 |
          v                          |
       Enc(A*s)                      |
          |                          |
          +---------- B - A*s -------+
                     |
                     v
          polynomial modular reduction
                     |
                     v
              CKKS ciphertext
```

That is the essential architecture of OpenFHE's FHEW → CKKS scheme switch.

---

# 39. One-Line Meaning of Each Important Function

```text
EvalSchemeSwitchingKeyGen()
    = generate the key material needed for scheme switching, including
      the CKKS encryption of the FHEW secret key

EvalFHEWtoCKKS()
    = perform the complete FHEW/LWE -> CKKS conversion

EvalPartialHomDecryption()
    = homomorphically evaluate A*s using the encrypted LWE secret key

EvalLTRectPrecomputeSwitch()
    = convert the LWE coefficient matrix A into a diagonal/BSGS form

EvalLTRectWithPrecomputeSwitch()
    = evaluate that diagonal linear transform on Enc(s)

EvalFastRotationPrecompute()
    = prepare reusable data for several rotations

EvalFastRotationExt()
    = perform encrypted CKKS slot rotations

KeySwitchExt()
    = prepare extended-RNS ciphertext representation

EvalMultExt()
    = diagonal plaintext × rotated encrypted secret key

EvalAddExtInPlace()
    = accumulate the diagonal products

KeySwitchDown()
    = return from the extended representation

MakeCKKSPackedPlaintext()
    = encode B as a CKKS plaintext

EvalNegate()
    = negate Enc(A*s)

EvalAdd()
    = compute Enc(B - A*s)

EvalChebyshevSeries()
    = homomorphically evaluate the polynomial approximation used for
      modular/sine reduction

EvalMult()
    = perform final numerical post-scaling

EvalAddInPlace()
    = apply final bias or combine repeated sparse blocks

ModReduceInPlace()
    = consume a CKKS modulus level when required
```

---

# 40. Final Mental Model

When reading `EvalFHEWtoCKKS()`, think:

```text
I have LWE ciphertexts:

    (a_i, b_i)

with:

    b_i = a_i · s + Δm_i + e_i   (mod q)

I want the messages m_i in CKKS.

1. Put all a_i vectors into matrix A.
2. Put all b_i values into B.
3. Encrypt the FHEW secret key s as a CKKS ciphertext.
4. Homomorphically evaluate A*s using a diagonal linear transform.
5. Compute B - A*s.
6. Apply the polynomial/sine-based modular reduction.
7. Apply post-scaling and post-bias.
8. Return one CKKS ciphertext containing the converted values.
```

The central mathematical identity is therefore:

```text
b - a · s = Δm + e   (mod q)
```

and the central implementation identity is:

```text
A*s
  = Σk Dk ⊙ Rotk(s)
```

with `s` represented by the encrypted CKKS switching key
`m_FHEWtoCKKSswk`.
