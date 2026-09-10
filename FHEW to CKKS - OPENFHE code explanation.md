# OpenFHE FHEW → CKKS Scheme Switching

This document explains **the FHEW → CKKS path** in OpenFHE's
`ckksrns-schemeswitching.cpp`, including the function-call tree,
important code operations, and their mathematical interpretation.

The main runtime function is:

```cpp
SWITCHCKKSRNS::EvalFHEWtoCKKS(...)
```

---

# 1. Big-Picture Function Call Tree

```text
EvalFHEWtoCKKS()
|
+-- EvalPartialHomDecryption()
|   |
|   +-- EvalLTRectPrecomputeSwitch()
|   |
|   +-- EvalLTRectWithPrecomputeSwitch()
|
+-- MakeCKKSPackedPlaintext()
|
+-- EvalNegate()
|
+-- EvalAdd()
|
+-- EvalChebyshevSeries()
|
+-- EvalMult()
|
+-- EvalAddInPlace()
|
+-- EvalSubInPlace()
|
+-- MakeCKKSPackedPlaintext()
|
+-- EvalMult()
|
+-- EvalAddInPlace()
|
+-- EvalAtIndex()
|
+-- EvalAddInPlace()
|
+-- ModReduceInPlace()
|
+-- return CKKS ciphertext
```

A useful conceptual view is:

```text
FHEW/LWE ciphertexts
        |
        v
extract LWE a-vectors and b-values
        |
        v
form A and b
        |
        v
homomorphically evaluate A * s
        |
        v
construct B - A*s
        |
        v
remove the LWE modulus periodicity
using a polynomial / sine approximation
        |
        v
post-scale and post-bias
        |
        v
CKKS ciphertext
```

---

# 2. Why Is FHEW → CKKS Needed?

An LWE/FHEW ciphertext is of the form:

```text
(a, b)
```

where the encrypted message is contained in the LWE relation

```text
b = a · s + Δm + e   (mod q)
```

Therefore:

```text
b - a · s = Δm + e   (mod q)
```

The problem is that `s` is secret.

In an ordinary LWE decryption procedure, one would explicitly compute:

```text
b - a · s
```

using the secret key.

During scheme switching, we do **not** decrypt the FHEW ciphertext first.
Instead, OpenFHE moves the computation into CKKS and evaluates the
decryption expression homomorphically.

The high-level idea is:

```text
LWE ciphertexts
     |
     v
A and b
     |
     v
homomorphically compute A*s
     |
     v
B - A*s
     |
     v
recover the encoded message values
     |
     v
CKKS ciphertext
```

This is the essential purpose of `EvalFHEWtoCKKS()`.

---

# 3. The LWE Equation

For each FHEW/LWE ciphertext:

```text
(a, b)
```

the LWE relation is:

```text
b = a · s + Δm + e   (mod q)
```

where:

```text
a = LWE ciphertext vector
s = LWE secret key
m = plaintext message
Δ = message scaling factor
e = LWE error
q = LWE modulus
```

Therefore:

```text
b - a · s = Δm + e   (mod q)
```

The left-hand side is the LWE decryption expression.

For multiple LWE ciphertexts, OpenFHE collects the `a` vectors into a
matrix `A` and the `b` values into a vector `b`.

Conceptually:

```text
A =
[ a0^T
  a1^T
  a2^T
  ...
  ak-1^T ]

b =
[ b0
  b1
  b2
  ...
  bk-1 ]
```

Then:

```text
A*s =
[ a0 · s
  a1 · s
  a2 · s
  ...
  ak-1 · s ]
```

and therefore:

```text
b - A*s
```

contains the corresponding LWE decryption expressions.

---

# 4. What Does `A` Mean Here?

This is an important distinction from the CKKS → FHEW direction.

In `EvalFHEWtoCKKS()`, the code creates:

```cpp
std::vector<std::vector<std::complex<double>>> A(numValues);
```

Here `A` really is the matrix formed from the **LWE `a` vectors**.

For each LWE ciphertext:

```cpp
auto& a = LWECiphertexts[i]->GetA();
```

and then:

```cpp
A[i][j] = std::complex<double>(a[j].ConvertToDouble(), 0);
```

So conceptually:

```text
LWE ciphertext 0: (a0, b0)
LWE ciphertext 1: (a1, b1)
LWE ciphertext 2: (a2, b2)
...

          |
          v

A =
[ a0
  a1
  a2
  ... ]
```

Here `A` is therefore a genuine matrix of LWE coefficients.

This is different from the `A` parameter used by
`EvalLTWithPrecomputeSwitch()` in the CKKS → FHEW direction, where the
parameter is a collection of precomputed plaintext diagonals.

In FHEW → CKKS:

```text
A before linear-transform precomputation
    = LWE coefficient matrix
```

After:

```text
EvalLTRectPrecomputeSwitch(A, ...)
```

the matrix has effectively been converted into a diagonalized plaintext
representation suitable for homomorphic matrix multiplication.

---

# 5. `EvalFHEWtoCKKS()`

## Purpose

This is the main FHEW → CKKS scheme-switching function.

The source follows three major conceptual steps before the final modular
reduction:

```text
Step 1:
construct A and b from LWE ciphertexts

Step 2:
homomorphically compute A*s

Step 3:
homomorphically compute B - A*s
```

The code explicitly labels these stages.

The overall mathematical flow is:

```text
LWE ciphertexts
        |
        v
A, b
        |
        v
A*s
        |
        v
b - A*s
        |
        v
encoded message
```

---

# 6. `numValues`, `slots`, and `n`

The function first determines how many LWE ciphertexts will be packed into
the resulting CKKS ciphertext.

Conceptually:

```text
numValues
    = number of LWE ciphertexts that will be switched
```

and:

```text
slots
    = number of CKKS slots available for the result
```

The source also obtains:

```cpp
uint32_t n = LWECiphertexts[0]->GetA().GetLength();
```

Here:

```text
n = LWE lattice parameter
  = length of each LWE a-vector
```

Thus if each LWE ciphertext is:

```text
(a, b)

a = [a0, a1, ..., an-1]
```

then:

```text
n = length(a)
```

---

# 7. Prescaling

The code computes:

```cpp
const double prescale =
    (1.0 / LWECiphertexts[0]->GetModulus().ConvertToDouble()) / K;
```

So conceptually:

```text
prescale = 1 / (q * K)
```

where:

```text
q = LWE modulus
K = scaling factor used by the modular-reduction approximation
```

The purpose is to put the LWE decryption quantity into the numerical range
expected by the CKKS polynomial approximation used later.

Conceptually:

```text
b - A*s
      |
      v
divide by q
      |
      v
normalized value
      |
      v
divide by K
      |
      v
range suitable for polynomial approximation
```

The source uses different `K` values depending on the LWE dimension.

---

# 8. Step 1: Form Matrix `A` and Vector `b`

The code creates:

```cpp
std::vector<std::vector<std::complex<double>>> A(numValues);
std::vector<std::complex<double>> b(b_size);
```

Then, for every LWE ciphertext:

```cpp
auto& a = LWECiphertexts[i]->GetA();
A[i].resize(a.GetLength());

for (uint32_t j = 0; j < a.GetLength(); ++j)
    A[i][j] =
        std::complex<double>(a[j].ConvertToDouble(), 0);

b[i] =
    std::complex<double>(
        prescale * LWECiphertexts[i]->GetB().ConvertToDouble(),
        0);
```

Conceptually:

```text
LWE ciphertext i:

(ai, bi)

        |
        +----> A[i] = ai
        |
        +----> b[i] = prescale * bi
```

After all LWE ciphertexts are processed:

```text
A =
[ a0^T
  a1^T
  ...
  a(numValues-1)^T ]

b =
[ prescale*b0
  prescale*b1
  ...
  prescale*b(numValues-1) ]
```

---

# 9. Why Is `A*s` a Linear Transform?

The product:

```text
A*s
```

has one output per LWE ciphertext.

For example:

```text
A =
[ a00 a01 a02 a03
  a10 a11 a12 a13
  a20 a21 a22 a23 ]

s =
[ s0
  s1
  s2
  s3 ]
```

Then:

```text
A*s =
[ a00*s0 + a01*s1 + a02*s2 + a03*s3
  a10*s0 + a11*s1 + a12*s2 + a13*s3
  a20*s0 + a21*s1 + a22*s2 + a23*s3 ]
```

The secret key `s` is encoded/encrypted in the CKKS-side switching key.

Therefore the computation becomes:

```text
encrypted secret-key representation
                |
                v
homomorphic matrix-vector multiplication
                |
                v
                  A*s
```

This is why `EvalPartialHomDecryption()` uses the same diagonalized
linear-transform machinery discussed in the FHE Textbook's matrix
multiplication formulation.

---

# 10. `EvalPartialHomDecryption()`

This function is the bridge between:

```text
LWE coefficient matrix A
```

and:

```text
homomorphic CKKS computation of A*s
```

The source first makes a copy:

```cpp
auto Acopy = A;
```

Then computes:

```cpp
size_t cols_po2 =
    1 << static_cast<uint32_t>(
        std::ceil(std::log2(A[0].size()))
    );
```

This finds the next power of two at least as large as the number of
columns.

If necessary, every row is padded:

```cpp
Acopy[i].resize(cols_po2);
```

Conceptually:

```text
original A:

[a00 a01 a02 a03 a04]

        |
        v

power-of-two padding

[a00 a01 a02 a03 a04 0 0 0]
```

---

# 11. Why Is `A` Padded to a Power of Two?

The CKKS linear-transform implementation is designed around power-of-two
dimensions and rotation structures.

The function comment explicitly states that the LWE lattice parameter is
padded to a power of two.

Conceptually:

```text
LWE dimension n
      |
      v
next power of two
      |
      v
matrix dimension accepted by the linear-transform routine
```

This makes the rotation/BSGS decomposition regular.

---

# 12. `EvalLTRectPrecomputeSwitch()`

The source then calls:

```cpp
auto Apre =
    EvalLTRectPrecomputeSwitch(Acopy, dim1, scale);
```

This converts the matrix into the diagonalized representation used by the
homomorphic linear-transform evaluator.

Conceptually:

```text
LWE matrix A
      |
      v
extract shifted diagonals
      |
      v
scale diagonals
      |
      v
organize for BSGS
      |
      v
Apre
```

The important difference from the square version is that this function is
specifically designed for the rectangular matrix that occurs here.

---

# 13. `EvalLTRectPrecomputeSwitch()`: Diagonal Extraction

For the relevant branch, the source does conceptually:

```cpp
auto tmp =
    ExtractShiftedDiagonal(A_slices[k], bStep*j + i);
```

Then:

```cpp
diag.insert(diag.end(), tmp.begin(), tmp.end());
```

and scales:

```cpp
elem * scale
```

So the process is:

```text
A matrix
 |
 +-- shifted diagonal 0 -> D0
 |
 +-- shifted diagonal 1 -> D1
 |
 +-- shifted diagonal 2 -> D2
 |
 +-- ...
 |
 v
diags[]
```

These diagonal vectors are the matrix-multiplication representation used by
the CKKS linear-transform engine.

---

# 14. Mathematical Form of the Linear Transform

For a matrix-vector multiplication:

```text
y = A*s
```

the diagonal decomposition is conceptually:

```text
y = D0 ⊙ Rot0(s)
  + D1 ⊙ Rot1(s)
  + D2 ⊙ Rot2(s)
  + ...
```

where:

```text
Dk
    = a cyclic/shifted diagonal of A

Rotk(s)
    = rotation of the CKKS-packed secret-key vector

⊙
    = slot-wise multiplication
```

The exact sign/direction of the rotation and the corresponding shifted
diagonal follows the convention implemented by OpenFHE's
`ExtractShiftedDiagonal()` and rotation helpers.

---

# 15. `m_FHEWtoCKKSswk`

A crucial object is:

```text
m_FHEWtoCKKSswk
```

This is a CKKS ciphertext containing an encoding of the FHEW/LWE secret
key.

The setup code constructs the LWE secret-key vector, maps the LWE
representation into floating-point CKKS values, and encrypts it under the
CKKS public key.

Conceptually:

```text
LWE secret key s
      |
      v
encode s as CKKS plaintext
      |
      v
CKKS Encrypt(pk, s)
      |
      v
m_FHEWtoCKKSswk
```

This is the key reason the FHEW → CKKS direction can evaluate LWE
decryption without possessing the plaintext secret key during evaluation.

---

# 16. `EvalLTRectWithPrecomputeSwitch()`

This function actually evaluates:

```text
A*s
```

using the precomputed diagonals.

Its inputs are conceptually:

```text
Apre  = diagonalized representation of A
ct    = CKKS ciphertext containing the encoded/encrypted LWE secret key
```

Its output is:

```text
encrypted A*s
```

The call from `EvalPartialHomDecryption()` is:

```cpp
return EvalLTRectWithPrecomputeSwitch(
    cc,
    Apre,
    ct,
    (Acopy.size() < A[0].size()),
    dim1,
    L
);
```

---

# 17. BSGS in `EvalLTRectWithPrecomputeSwitch()`

The function calculates:

```cpp
uint32_t bStep =
    (dim1 == 0) ? getRatioBSGSLT(n) : dim1;

uint32_t gStep =
    std::ceil(static_cast<double>(n) / bStep);
```

where:

```text
bStep = baby-step size
gStep = number of giant steps
```

The diagonal index is conceptually decomposed as:

```text
k = j*bStep + i
```

where:

```text
j = giant-step index
i = baby-step index
```

This is the same BSGS strategy used by the CKKS linear-transform machinery.

---

# 18. `EvalFastRotationPrecompute()`

The function computes:

```cpp
auto digits = cc.EvalFastRotationPrecompute(ct);
```

The comment explains that this computes the NTT representation needed for
the hoisted automorphisms used later.

Conceptually:

```text
encrypted secret-key vector
        |
        v
precompute common rotation information
        |
        v
digits
```

This is an optimization and does not change the mathematical operation:

```text
A*s
```

---

# 19. `EvalFastRotationExt()`

The function prepares baby-step rotations:

```cpp
fastRotation[j - 1] =
    cc.EvalFastRotationExt(
        ct,
        j,
        digits,
        true
    );
```

Conceptually:

```text
ct = encrypted s

EvalFastRotationExt(ct, j, ...)
            |
            v
        Rotj(s)
```

So the mathematical mapping is:

```text
FHE Textbook / linear transform      OpenFHE
------------------------------------------------
Rotj(s)                              EvalFastRotationExt(...)
```

---

# 20. `KeySwitchExt()`

The first giant-step computation contains:

```cpp
cc.KeySwitchExt(ctxt, true)
```

This is not another matrix term.

It is part of the extended-RNS ciphertext machinery used by
`EvalMultExt()`.

Conceptually:

```text
CKKS ciphertext
      |
      v
extended representation
      |
      v
extended plaintext multiplication
```

The mathematical operation remains:

```text
Dk ⊙ Rotk(s)
```

---

# 21. `EvalMultExt()`

The core operation is:

```cpp
auto inner =
    FHECKKSRNS::EvalMultExt(
        cc.KeySwitchExt(ctxt, true),
        A[bStep * j]
    );
```

For baby-step terms:

```cpp
FHECKKSRNS::EvalAddExtInPlace(
    inner,
    FHECKKSRNS::EvalMultExt(
        fastRotation[i - 1],
        A[bStep * j + i]
    )
);
```

Conceptually:

```text
A[k]
   |
   v
Dk

fastRotation[i-1]
   |
   v
Rotk(s)

EvalMultExt(...)
   |
   v
Dk ⊙ Rotk(s)
```

---

# 22. `EvalAddExtInPlace()`

The individual diagonal contributions are accumulated:

```cpp
EvalAddExtInPlace(
    inner,
    EvalMultExt(...)
);
```

Conceptually:

```text
inner =
    inner + Dk ⊙ Rotk(s)
```

After all terms:

```text
A*s =
D0 ⊙ Rot0(s)
+ D1 ⊙ Rot1(s)
+ D2 ⊙ Rot2(s)
+ ...
```

Therefore:

```text
EvalLTRectWithPrecomputeSwitch()
    =
evaluate A*s homomorphically
```

---

# 23. Giant-Step Rotation

After the inner baby-step sum has been formed, the code handles the giant
step.

It identifies the automorphism associated with:

```cpp
bStep * j
```

using:

```cpp
FindAutomorphismIndex2nComplex(
    bStep * j,
    M
)
```

Then it evaluates:

```cpp
cc.EvalFastRotationExt(
    inner,
    bStep * j,
    innerDigits,
    false
)
```

Conceptually:

```text
baby-step accumulation
        |
        v
rotate by the giant-step amount
        |
        v
giant-step contribution
```

This is how the BSGS decomposition reconstructs the full set of required
rotations.

---

# 24. `KeySwitchDownFirstElement()` and `KeySwitchDown()`

The implementation separates one polynomial component:

```cpp
first =
    cc.KeySwitchDownFirstElement(inner);
```

and later uses:

```cpp
inner =
    cc.KeySwitchDown(inner);
```

and finally:

```cpp
result =
    cc.KeySwitchDown(result);
```

Conceptually:

```text
extended-RNS intermediate result
        |
        +-- KeySwitchDownFirstElement()
        |
        +-- KeySwitchDown()
        |
        v
normal CKKS representation
```

These are implementation-level operations required by the extended
ciphertext machinery.

They are not additional terms in:

```text
A*s
```

---

# 25. Result of `EvalPartialHomDecryption()`

After:

```cpp
EvalLTRectWithPrecomputeSwitch(...)
```

the result is approximately:

```text
AdotS = A*s
```

but evaluated **homomorphically in CKKS**.

This is the critical conceptual transition:

```text
LWE secret key s
      |
      | encrypted under CKKS
      v
CKKS ciphertext of s
      |
      | linear transform A
      v
CKKS ciphertext of A*s
```

No LWE plaintext is decrypted during this process.

---

# 26. Step 3: Construct `BPlain`

After `AdotS` has been computed, the code constructs:

```cpp
Plaintext BPlain =
    ccCKKS->MakeCKKSPackedPlaintext(
        b,
        AdotS->GetNoiseScaleDeg(),
        AdotS->GetLevel(),
        nullptr,
        N / 2
    );
```

Conceptually:

```text
b vector
  |
  v
CKKS plaintext
  |
  v
BPlain
```

The important point is that `BPlain` is **not encrypted**.

The LWE `b` values are ciphertext components, not the plaintext messages.
Once the `b` values have been extracted from the LWE ciphertexts, they can
be inserted into a CKKS plaintext.

---

# 27. Compute `B - A*s`

The source performs:

```cpp
auto BminusAdotS =
    ccCKKS->EvalAdd(
        ccCKKS->EvalNegate(AdotS),
        BPlain
    );
```

Mathematically:

```text
BminusAdotS = B - A*s
```

Using the LWE equation:

```text
b = a · s + Δm + e   (mod q)
```

we obtain:

```text
b - a · s = Δm + e   (mod q)
```

For all ciphertexts together:

```text
B - A*s
    =
[ Δm0 + e0
  Δm1 + e1
  ...
]
```

up to the normalization/scaling introduced by `prescale`.

This is the central mathematical step of the FHEW → CKKS switch.

---

# 28. Why Is Another Step Needed After `B - A*s`?

The quantity:

```text
B - A*s
```

is still fundamentally an LWE-modulus quantity.

The LWE modulus is periodic:

```text
value mod q
```

while CKKS operates over approximate real/complex values.

OpenFHE therefore applies a polynomial approximation to recover a
representative of the desired message from the modular quantity.

Conceptually:

```text
B - A*s
      |
      v
normalize
      |
      v
approximate modular reduction
      |
      v
message-related value
```

---

# 29. `EvalChebyshevSeries()`

The code uses:

```cpp
auto BminusAdotS3 =
    ccCKKS->EvalChebyshevSeries(
        BminusAdotS,
        coefficientsFHEW,
        -1.0,
        1.0
    );
```

The coefficient array is chosen according to the LWE parameter and desired
precision.

Conceptually:

```text
x = normalized (B - A*s)

        |
        v

Chebyshev polynomial P(x)

        |
        v

approximation to the required periodic/modular function
```

So:

```text
EvalChebyshevSeries()
    =
evaluate an encrypted polynomial approximation
```

This is performed entirely homomorphically in CKKS.

---

# 30. Why Use a Sine-Based Function?

The source comment describes the modular-reduction stage as:

```text
homomorphically evaluate modular function
using sine approximation
```

The coefficient tables correspond to approximations of a periodic function
related to:

```text
sin(2*pi*x)
```

The purpose is to map values that differ by multiples of the LWE modulus
back toward the representative corresponding to the desired plaintext.

Conceptually:

```text
q-periodic LWE value
        |
        v
normalized periodic value
        |
        v
approximate periodic function
        |
        v
recover a bounded representative
```

---

# 31. Repeated Squaring / Sine-Refinement Loop

The source then executes:

```cpp
const int32_t BT_ITER = 3;

for (int32_t j = 1; j <= BT_ITER; ++j) {
    BminusAdotS3 =
        ccCKKS->EvalMult(
            BminusAdotS3,
            BminusAdotS3
        );

    ccCKKS->EvalAddInPlace(
        BminusAdotS3,
        BminusAdotS3
    );

    double scalar =
        1.0 /
        std::pow(
            (2.0 * M_PI),
            std::pow(2.0, j - BT_ITER)
        );

    ccCKKS->EvalSubInPlace(
        BminusAdotS3,
        scalar
    );
}
```

The conceptual recurrence is:

```text
x(j) -> 2*x(j)^2 - c(j)
```

with a stage-dependent constant:

```text
c(j) =
1 / (2*pi)^(2^(j-BT_ITER))
```

The purpose is to repeatedly refine the approximation to the desired
periodic/modular transformation.

So the code is not simply "doing another matrix multiplication." This
section is the **modular-reduction / message-recovery stage**.

---

# 32. Why Does `K` Appear?

The code selects:

```text
K = 16
```

for one LWE parameter choice and:

```text
K = 128
```

for another.

The prescale is:

```text
prescale = 1 / (q*K)
```

Conceptually:

```text
B - A*s
      |
      | divide by q
      v
message-scale quantity
      |
      | divide by K
      v
small normalized interval
      |
      v
Chebyshev / sine approximation
```

The reason for doing this before the polynomial approximation is to place the
input into the interval for which the stored polynomial approximation was
generated.

---

# 33. Post-Scaling

After the modular-reduction polynomial has been evaluated, the source
computes:

```cpp
double postScale =
    (p >= 1 && p <= 4) ? (2.0 * M_PI) : static_cast<double>(p);
```

and, when needed:

```cpp
postScale *= (pmax - pmin) / 4.0;
```

Conceptually:

```text
normalized recovered value
        |
        v
multiply by postScale
        |
        v
restore the desired message scale
```

The exact factor depends on the plaintext modulus and output encoding.

---

# 34. Post-Bias

If a nonzero `pmin` is used, the code also computes:

```cpp
postBias = (pmax - pmin) / 4.0;
```

Then creates:

```cpp
postBiasPlain =
    ccCKKS->MakeCKKSPackedPlaintext(
        postBiasVec,
        ...
    );
```

and adds:

```cpp
ccCKKS->EvalAddInPlace(
    BminusAdotSres,
    postBiasPlain
);
```

Conceptually:

```text
scaled value
    +
offset
    =
final encoded CKKS value
```

This is needed when the target output interval is not centered at zero.

---

# 35. Restore Sparse Encoding

If the CKKS context uses sparse packing:

```cpp
if (isSparse) {
    for (...) {
        auto temp =
            ccCKKS->EvalAtIndex(
                BminusAdotSres,
                j * slots
            );

        ccCKKS->EvalAddInPlace(
            BminusAdotSres,
            temp
        );
    }

    BminusAdotSres->SetSlots(slots);
}
```

Conceptually:

```text
full / repeated intermediate encoding
        |
        v
sum the repeated blocks
        |
        v
original sparse-slot organization
```

This is encoding-layout restoration, not message computation.

---

# 36. Complete FHEW → CKKS Tree

```text
SWITCHCKKSRNS::EvalFHEWtoCKKS()
|
+-- determine:
|     +-- numValues
|     +-- slots
|     +-- LWE dimension n
|     +-- K
|     +-- prescale
|
+-- form A
|     |
|     +-- copy every LWE a-vector
|
+-- form b
|     |
|     +-- prescale every LWE b-value
|
+-- EvalPartialHomDecryption()
|   |
|   +-- pad A columns to power of two
|   |
|   +-- EvalLTRectPrecomputeSwitch()
|   |   |
|   |   +-- ExtractShiftedDiagonal()
|   |   +-- scale
|   |   +-- arrange BSGS diagonals
|   |
|   +-- EvalLTRectWithPrecomputeSwitch()
|       |
|       +-- calculate bStep / gStep
|       |
|       +-- EvalFastRotationPrecompute()
|       |
|       +-- EvalFastRotationExt()
|       |       -> rotate encrypted LWE secret key
|       |
|       +-- KeySwitchExt()
|       |
|       +-- EvalMultExt()
|       |       -> Dk ⊙ Rotk(s)
|       |
|       +-- EvalAddExtInPlace()
|       |       -> sum diagonal terms
|       |
|       +-- giant-step automorphisms
|       |
|       +-- KeySwitchDownFirstElement()
|       |
|       +-- KeySwitchDown()
|       |
|       v
|       AdotS = A*s
|
+-- MakeCKKSPackedPlaintext(b)
|       |
|       v
|       BPlain
|
+-- EvalNegate(AdotS)
|
+-- EvalAdd(BPlain, -AdotS)
|       |
|       v
|       B - A*s
|
+-- EvalChebyshevSeries(...)
|       |
|       v
|       polynomial approximation
|
+-- repeated:
|     +-- EvalMult(x, x)
|     +-- EvalAddInPlace(x, x)
|     +-- EvalSubInPlace(x, scalar)
|     |
|     v
|     refined modular reduction
|
+-- postScale
|
+-- postBias
|
+-- restore sparse encoding if needed
|
+-- final ModReduce
|
v
CKKS ciphertext
```

---

# 37. Mathematical View of the Same Tree

The whole FHEW → CKKS path can be summarized mathematically as:

```text
LWE ciphertexts:

bi = ai · s + Δmi + ei   (mod q)

        |
        v

A =
[ a0^T
  a1^T
  ...
]

b =
[ b0
  b1
  ...
]

        |
        v

homomorphically evaluate:

A*s

        |
        v

construct:

b - A*s

        |
        v

approximately obtain:

Δm + e   (mod q)

        |
        v

normalize by q and K

        |
        v

evaluate polynomial approximation
to the required periodic/modular function

        |
        v

post-scale + post-bias

        |
        v

CKKS ciphertext containing
the recovered message values
```

---

# 38. Exact Conceptual Mapping of the Important Operations

```text
Mathematical operation                       OpenFHE operation
====================================================================

collect LWE ai vectors                      GetA()

collect LWE bi values                       GetB()

construct LWE coefficient matrix A          A[i][j] = a[j]

homomorphically encode secret key s          m_FHEWtoCKKSswk

matrix A times secret key s                  EvalPartialHomDecryption()

pad A for power-of-two transform             Acopy[i].resize(cols_po2)

diagonalize A                                EvalLTRectPrecomputeSwitch()

extract shifted diagonal                    ExtractShiftedDiagonal()

encode/use transform diagonal                Apre / plaintext diagonal

rotate encrypted s                           EvalFastRotationExt()

plaintext diagonal × rotated s               EvalMultExt()

sum diagonal contributions                  EvalAddExtInPlace()

giant-step automorphism                     FindAutomorphismIndex2nComplex()

return from extended representation          KeySwitchDown*()

construct CKKS plaintext for B              MakeCKKSPackedPlaintext()

negate A*s                                   EvalNegate()

compute B - A*s                             EvalAdd(...)

approximate modular/periodic function        EvalChebyshevSeries()

polynomial refinement                        EvalMult + EvalAdd + EvalSub

restore message scale                       postScale

restore shifted output range                 postBias

restore sparse packing                       EvalAtIndex + EvalAddInPlace

final CKKS level reduction                   ModReduce...
```

---

# 39. The Most Important Distinction: `A*s` Is Homomorphic

The most important conceptual point in this whole function is:

```text
The code does NOT decrypt the LWE ciphertexts.
```

It knows:

```text
A = public LWE a-vectors
B = public LWE b-values
```

but `s` remains secret.

Instead, OpenFHE has:

```text
encrypted CKKS representation of s
```

in:

```text
m_FHEWtoCKKSswk
```

Therefore:

```text
A*s
```

is evaluated homomorphically.

The flow is:

```text
             encrypted s
                 |
                 v
       +---------------------+
       |   linear transform  |
       |                     |
       |       A*s           |
       +---------------------+
                 |
                 v
          encrypted A*s
```

This is why `EvalLTRectWithPrecomputeSwitch()` is required in the
FHEW → CKKS path.

---

# 40. The Core Equation to Keep in Mind

For the FHEW → CKKS direction, keep this equation in mind:

```text
b = A*s + Δm + e   (mod q)
```

Therefore:

```text
b - A*s = Δm + e   (mod q)
```

OpenFHE implements this as:

```text
LWE a-vectors
      |
      v
      A
      |
      | homomorphic linear transform
      v
   encrypted A*s
      |
      +------------------+
                         |
b --------------------> subtract
                         |
                         v
                    B - A*s
                         |
                         v
              polynomial / sine-based
                 modular reduction
                         |
                         v
                  recovered message
                         |
                         v
                    CKKS ciphertext
```

---

# 41. One-Line Explanation of Each Function

```text
EvalFHEWtoCKKS()
    = perform the complete FHEW -> CKKS scheme switch

EvalPartialHomDecryption()
    = homomorphically evaluate A*s

EvalLTRectPrecomputeSwitch()
    = convert the LWE coefficient matrix A into
      a diagonalized linear-transform representation

EvalLTRectWithPrecomputeSwitch()
    = evaluate that linear transform on the encrypted LWE secret key

EvalFastRotationPrecompute()
    = prepare reusable data for multiple CKKS rotations

EvalFastRotationExt()
    = homomorphically rotate the encrypted LWE secret key

KeySwitchExt()
    = prepare/use the extended-RNS ciphertext representation

EvalMultExt()
    = multiply a rotated secret-key ciphertext by a transform diagonal

EvalAddExtInPlace()
    = accumulate the diagonal contributions

KeySwitchDownFirstElement()
    = extract/restore the first component from the extended representation

KeySwitchDown()
    = return the result from the extended representation

MakeCKKSPackedPlaintext()
    = encode the LWE b-values as a CKKS plaintext

EvalNegate()
    = negate the encrypted A*s result

EvalAdd()
    = compute B - A*s

EvalChebyshevSeries()
    = evaluate the polynomial approximation used for modular reduction

EvalMult()
    = perform encrypted polynomial multiplication

EvalAddInPlace()
    = accumulate/add encrypted polynomial terms

EvalSubInPlace()
    = subtract the stage-dependent constant

ModReduce()
    = reduce the CKKS modulus-chain level

EvalAtIndex()
    = rotate/access slots used to restore sparse packing

m_FHEWtoCKKSswk
    = CKKS ciphertext containing the encoded/encrypted FHEW/LWE secret key
```

---

# 42. Final Mental Model

When reading:

```cpp
EvalFHEWtoCKKS(...)
```

think:

```text
"I have LWE ciphertexts:

    (ai, bi)

satisfying:

    bi = ai*s + Δmi + ei  (mod q)

I construct:

    A = [a0^T; a1^T; ...]

and:

    b = [b0, b1, ...]

I already have the LWE secret key s encrypted in CKKS.

So I homomorphically compute:

    A*s

using a diagonalized linear transform.

Then I compute:

    b - A*s

which is the LWE decryption expression.

Finally I use the polynomial/sine-based modular reduction,
post-scale, and post-bias to recover values in CKKS slots."
```

The central mathematical correspondence is therefore:

```text
LWE decryption:
    b - a*s

Multiple LWE ciphertexts:
    B - A*s

OpenFHE implementation:
    EvalPartialHomDecryption()
            |
            v
          A*s
            |
            v
    EvalAdd(BPlain, -AdotS)
            |
            v
        B - A*s
            |
            v
    EvalChebyshevSeries()
            |
            v
       recovered CKKS values
```

The key point is that **FHEW → CKKS is fundamentally a homomorphic
evaluation of the LWE decryption expression**, followed by an
approximate modular-reduction procedure that converts the result into a
usable CKKS representation.
