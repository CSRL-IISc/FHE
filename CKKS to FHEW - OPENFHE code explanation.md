# OpenFHE CKKS Scheme Switching — Function Call Tree and Mathematical Mapping

## 1. Overview

The CKKS → FHEW scheme-switching path can be viewed as:

```text
CKKS ciphertext
      |
      v
EvalCKKStoFHEW()
      |
      +-- EvalSlotsToCoeffsSwitch()
      |       |
      |       +-- EvalLTWithPrecomputeSwitch()
      |               |
      |               +-- rotations
      |               +-- plaintext diagonal multiplication
      |               +-- additions
      |
      +-- ModSwitch()
      |
      +-- KeySwitch()
      |
      +-- ExtractLWEpacked()
      |
      +-- ExtractLWECiphertext()
      |
      +-- RoundqQAlter()
      |
      v
LWE / FHEW ciphertexts
```

The important idea is that CKKS data is initially represented in its
**slots**, while the later extraction step needs the information in
specific **polynomial coefficient positions**.

The linear transform performed by
`EvalLTWithPrecomputeSwitch()` is the mechanism used to rearrange the
encrypted data into the required representation.

---

# 2. The Mathematical Idea: Linear Transform

Suppose the CKKS slots contain:

```text
x = [x0, x1, x2, x3]
```

and we want to apply a linear transformation:

```text
y = T · x
```

where `T` is a matrix.

For example:

```text
T =
[ t00 t01 t02 t03
  t10 t11 t12 t13
  t20 t21 t22 t23
  t30 t31 t32 t33 ]
```

Then:

```text
y0 = t00*x0 + t01*x1 + t02*x2 + t03*x3
y1 = t10*x0 + t11*x1 + t12*x2 + t13*x3
y2 = t20*x0 + t21*x1 + t22*x2 + t23*x3
y3 = t30*x0 + t31*x1 + t32*x2 + t33*x3
```

Instead of evaluating this matrix multiplication directly, CKKS uses
a diagonal decomposition.

---

# 3. Diagonal Representation

The matrix can be represented using cyclic diagonals.

For the example above, one possible diagonal representation is:

```text
D0 = [t00, t11, t22, t33]

D1 = [t01, t12, t23, t30]

D2 = [t02, t13, t20, t31]

D3 = [t03, t10, t21, t32]
```

The transform can then be written conceptually as:

```text
y = D0 ⊙ Rot0(x)
  + D1 ⊙ Rot1(x)
  + D2 ⊙ Rot2(x)
  + D3 ⊙ Rot3(x)
```

where:

```text
⊙       = slot-wise multiplication
Rotk(x) = cyclic rotation of the CKKS slots by k
```

For example:

```text
x = [x0, x1, x2, x3]

Rot0(x) = [x0, x1, x2, x3]
Rot1(x) = [x1, x2, x3, x0]
Rot2(x) = [x2, x3, x0, x1]
Rot3(x) = [x3, x0, x1, x2]
```

Therefore:

```text
D0 ⊙ Rot0(x)
= [t00*x0, t11*x1, t22*x2, t33*x3]

D1 ⊙ Rot1(x)
= [t01*x1, t12*x2, t23*x3, t30*x0]

D2 ⊙ Rot2(x)
= [t02*x2, t13*x3, t20*x0, t31*x1]

D3 ⊙ Rot3(x)
= [t03*x3, t10*x0, t21*x1, t32*x2]
```

Adding them gives:

```text
y0 = t00*x0 + t01*x1 + t02*x2 + t03*x3

y1 = t10*x0 + t11*x1 + t12*x2 + t13*x3

y2 = t20*x0 + t21*x1 + t22*x2 + t23*x3

y3 = t30*x0 + t31*x1 + t32*x2 + t33*x3
```

Therefore:

```text
y = T · x
```

This is the diagonalized matrix-vector multiplication described in the
FHE Textbook's matrix multiplication formulation.

---

# 4. Important Meaning of `A` in OpenFHE

There are two different things that can be called `A`.

## 4.1 `A` in ordinary matrix notation

In:

```text
y = A · x
```

`A` means the original dense transformation matrix.

## 4.2 `A` in `EvalLTWithPrecomputeSwitch()`

In OpenFHE:

```cpp
EvalLTWithPrecomputeSwitch(cc, ctxt, A, dim1)
```

the parameter `A` is a vector of **precomputed plaintexts**.

These plaintexts represent the diagonals of the transformation matrix.

Conceptually:

```text
Original matrix T
      |
      +-- extract shifted diagonal
      |
      v
D0, D1, D2, D3, ...
      |
      +-- MakeAuxPlaintext()
      |
      v
A[0], A[1], A[2], A[3], ...
```

Thus:

```text
A[0]  -> plaintext encoding of D0
A[1]  -> plaintext encoding of D1
A[2]  -> plaintext encoding of D2
...
```

So `A[k]` in:

```cpp
A[bStep*j + i]
```

should be interpreted as:

```text
the plaintext containing one transformation-matrix diagonal
```

not as an individual matrix element.

---

# 5. `EvalLTPrecomputeSwitch()`

## Purpose

`EvalLTPrecomputeSwitch()` prepares the plaintext diagonal representation
of the linear transformation.

Its conceptual flow is:

```text
Transformation matrix T
        |
        v
ExtractShiftedDiagonal(...)
        |
        v
D0, D1, D2, ...
        |
        v
scale / arrange coefficients
        |
        v
MakeAuxPlaintext(...)
        |
        v
A[0], A[1], A[2], ...
```

## Important operation

The source contains:

```cpp
ExtractShiftedDiagonal(newA, ji)
```

This extracts a shifted/cyclic diagonal from the transformation matrix.

Then the diagonal is converted into a CKKS plaintext through:

```cpp
FHECKKSRNS::MakeAuxPlaintext(...)
```

Therefore:

```text
EvalLTPrecomputeSwitch()
        =
construct plaintext diagonals needed by the linear transform
```

It does not yet apply the transform to the encrypted ciphertext.

---

# 6. `EvalLTWithPrecomputeSwitch()`

## Purpose

This function takes:

```text
1. the encrypted CKKS ciphertext
2. the precomputed diagonal plaintexts
```

and actually evaluates the linear transformation.

Conceptually:

```text
                 ciphertext x
                       |
          +------------+------------+
          |            |            |
          v            v            v
        Rot0          Rot1         Rot2 ...
          |            |            |
          v            v            v
       × D0          × D1         × D2
          |            |            |
          +------------+------------+
                       |
                       v
                     add
                       |
                       v
                       y
```

The mathematical form is:

```text
y = Σk Dk ⊙ Rotk(x)
```

or, for the 4-slot example:

```text
y = D0 ⊙ Rot0(x)
  + D1 ⊙ Rot1(x)
  + D2 ⊙ Rot2(x)
  + D3 ⊙ Rot3(x)
```

---

# 7. Mapping the OpenFHE Operations

## 7.1 `EvalFastRotationPrecompute()`

OpenFHE:

```cpp
digits = cc.EvalFastRotationPrecompute(ctxt);
```

Conceptually:

```text
prepare information needed for multiple ciphertext rotations
```

This is an optimization.

It does not represent a new mathematical term in the matrix multiplication.

---

# 8. `EvalFastRotationExt()`

OpenFHE:

```cpp
cc.EvalFastRotationExt(ctxt, i, digits, true)
```

Conceptually:

```text
Rot_i(x)
```

That is, it performs a homomorphic CKKS slot rotation.

So:

```text
FHE Textbook operation       OpenFHE operation
------------------------------------------------
Rot_i(x)                     EvalFastRotationExt(...)
```

---

# 9. `EvalMultExt()`

OpenFHE:

```cpp
FHECKKSRNS::EvalMultExt(
    fastRotation[i - 1],
    A[bStep*j + i]
)
```

Conceptually:

```text
Dk ⊙ Rotk(x)
```

where:

```text
Dk = corresponding plaintext diagonal
```

So:

```text
EvalFastRotationExt(...)
        |
        v
Rotk(x)

A[k]
        |
        v
Dk

EvalMultExt(...)
        |
        v
Dk ⊙ Rotk(x)
```

This is the slot-wise multiplication part of the diagonalized linear
transform.

---

# 10. `EvalAddExtInPlace()`

OpenFHE:

```cpp
FHECKKSRNS::EvalAddExtInPlace(
    inner,
    FHECKKSRNS::EvalMultExt(...)
);
```

Conceptually:

```text
accumulator = accumulator + Dk ⊙ Rotk(x)
```

After all terms are added:

```text
y = D0 ⊙ Rot0(x)
  + D1 ⊙ Rot1(x)
  + D2 ⊙ Rot2(x)
  + ...
```

Therefore:

```text
EvalAddExtInPlace()
        =
accumulate the diagonal contributions
```

---

# 11. Why BSGS Is Used

A direct implementation would require:

```text
Rot0(x)
Rot1(x)
Rot2(x)
Rot3(x)
...
RotN-1(x)
```

For a large number of CKKS slots, this is expensive.

OpenFHE therefore uses **baby-step/giant-step (BSGS)**.

The code contains:

```cpp
bStep
gStep
```

and indexes the diagonal plaintexts using:

```cpp
A[bStep*j + i]
```

The rotation index can therefore be viewed as:

```text
k = j*bStep + i
```

where:

```text
j = giant-step index
i = baby-step index
```

So:

```text
A[bStep*j + i]
```

means approximately:

```text
D_k
```

with:

```text
k = j*bStep + i
```

The mathematical operation remains:

```text
Dk ⊙ Rotk(x)
```

BSGS changes **how efficiently the rotations are generated and reused**.

---

# 12. `KeySwitchExt()`

The source also contains:

```cpp
cc.KeySwitchExt(ctxt, true)
```

This should not be interpreted as another matrix-multiplication term.

It is implementation machinery associated with the ciphertext representation
used by the extended-RNS operations.

Conceptually:

```text
linear-transform mathematics
        |
        +-- rotations
        +-- diagonal multiplication
        +-- additions
        |
        +-- key switching / RNS representation management
```

The key-switch operation supports the implementation of the transform but
does not correspond to an additional `Dk ⊙ Rotk(x)` term.

---

# 13. `KeySwitchDown()`

Similarly:

```cpp
cc.KeySwitchDown(...)
```

is associated with moving the ciphertext back from the extended
representation used during the computation.

Conceptually:

```text
extended representation
        |
        v
KeySwitchDown()
        |
        v
normal CKKS representation
```

Again, this is implementation detail rather than a new mathematical
operation in:

```text
y = Σk Dk ⊙ Rotk(x)
```

---

# 14. `EvalSlotsToCoeffsSwitch()`

This is the higher-level operation in the CKKS → FHEW path.

Its role is to perform the required linear transformation that converts
the CKKS slot-oriented representation into the representation needed for
the later coefficient extraction.

Conceptually:

```text
CKKS slot representation
        |
        v
linear transformation
        |
        v
coefficient-oriented representation
```

Internally it uses:

```text
EvalLTWithPrecomputeSwitch()
```

to apply the precomputed transformation.

Thus:

```text
EvalSlotsToCoeffsSwitch()
        |
        +-- prepare / adjust ciphertext
        |
        +-- EvalLTWithPrecomputeSwitch()
                |
                +-- rotations
                +-- diagonal multiplications
                +-- additions
```

---

# 15. Why This Transform Is Needed

CKKS logically exposes a vector of slots:

```text
[x0, x1, x2, ..., xN-1]
```

But the ciphertext itself is an RLWE polynomial object.

The values in the logical slots are not simply stored as:

```text
coefficient 0 = x0
coefficient 1 = x1
...
```

Therefore the scheme-switching process needs a linear transformation to
rearrange the encoded information.

The conceptual flow is:

```text
CKKS slots
    |
    |  linear transform
    v
desired polynomial-coefficient representation
    |
    v
extract polynomial coefficients
    |
    v
LWE ciphertext
```

This is why the linear transform appears before:

```cpp
ExtractLWEpacked(...)
```

and:

```cpp
ExtractLWECiphertext(...)
```

---

# 16. `EvalCKKStoFHEW()`

This is the main CKKS → FHEW runtime function.

The flow is:

```text
EvalCKKStoFHEW()
|
+-- EvalSlotsToCoeffsSwitch()
|      |
|      +-- EvalLTWithPrecomputeSwitch()
|
+-- ModSwitch()
|
+-- KeySwitch()
|
+-- ExtractLWEpacked()
|
+-- ExtractLWECiphertext()
|
+-- RoundqQAlter()
|
v
LWE ciphertexts
```

Each stage has a different purpose.

---

# 17. `ModSwitch()`

After the CKKS linear transform, the ciphertext is moved from the current
CKKS modulus to the modulus required by the scheme-switching operation.

Conceptually:

```text
CKKS ciphertext at Q
        |
        v
ModSwitch()
        |
        v
CKKS ciphertext at Q'
```

This changes the modulus representation.

It is not the linear transformation itself.

---

# 18. `KeySwitch()`

OpenFHE then performs:

```cpp
ccKS->KeySwitch(ctxtKS, m_CKKStoFHEWswk)
```

The purpose is to switch from the CKKS secret-key representation to the
RLWE representation corresponding to the FHEW/LWE secret key.

Conceptually:

```text
RLWE ciphertext under CKKS key
            |
            v
       KeySwitch
            |
            v
RLWE ciphertext under LWE-derived key
```

This makes the RLWE ciphertext suitable for coefficient extraction into
an LWE ciphertext.

---

# 19. `ExtractLWEpacked()`

This function takes the RLWE ciphertext and extracts its polynomial
coefficient vectors.

Conceptually:

```text
RLWE ciphertext

(A(X), B(X))

       |
       v

A = [A0, A1, A2, ...]
B = [B0, B1, B2, ...]
```

So:

```text
ExtractLWEpacked()
    =
extract coefficient representation of RLWE A and B
```

It does not yet create separate LWE ciphertexts.

---

# 20. `ExtractLWECiphertext()`

This function takes the packed coefficient representation and chooses
the coefficient positions corresponding to one LWE ciphertext.

Conceptually:

```text
packed RLWE coefficients
        |
        v
select appropriate positions
        |
        +-- a = [a0, a1, ..., an-1]
        |
        +-- b = selected B coefficient
        |
        v
LWE ciphertext (a,b)
```

The resulting LWE ciphertext has the usual form:

```text
(a, b)
```

---

# 21. LWE Equation

For an LWE ciphertext:

```text
b = a · s + Δm + e  (mod q)
```

where:

```text
a = LWE vector
s = LWE secret key
m = plaintext message
Δ = message scaling factor
e = error
q = LWE modulus
```

Decryption computes:

```text
b - a · s = Δm + e  (mod q)
```

The expression:

```text
b - a · s
```

is therefore the quantity from which the plaintext message is recovered.

---

# 22. `RoundqQAlter()`

The code may have:

```text
Q' != q
```

where:

```text
Q' = modulus used by the CKKS-side representation
q  = LWE modulus
```

Therefore the extracted LWE values have to be mapped to the LWE modulus.

The code performs this through:

```cpp
RoundqQAlter(...)
```

Conceptually:

```text
value modulo Q'
        |
        v
scale from Q' to q
        |
        v
round
        |
        v
value modulo q
```

This is the final modulus-conversion step before returning the FHEW/LWE
ciphertexts.

---

# 23. Complete CKKS → FHEW Tree

```text
EvalCKKStoFHEW()
|
+-- EvalSlotsToCoeffsSwitch()
|   |
|   +-- EvalLTWithPrecomputeSwitch()
|       |
|       +-- EvalFastRotationPrecompute()
|       |
|       +-- EvalFastRotationExt()
|       |       |
|       |       +-- Rot_k(x)
|       |
|       +-- KeySwitchExt()
|       |
|       +-- EvalMultExt()
|       |       |
|       |       +-- D_k ⊙ Rot_k(x)
|       |
|       +-- EvalAddExtInPlace()
|       |       |
|       |       +-- accumulate terms
|       |
|       +-- KeySwitchDown()
|
+-- ModSwitch()
|
+-- KeySwitch()
|
+-- ExtractLWEpacked()
|
+-- ExtractLWECiphertext()
|
+-- RoundqQAlter()
|
v
LWE ciphertexts
```

---

# 24. Mathematical View of the Whole CKKS → FHEW Process

At the highest level:

```text
CKKS slots
    |
    | y = T · x
    | where
    | y = Σk Dk ⊙ Rotk(x)
    |
    v
transformed CKKS representation
    |
    | modulus switch
    v
Q' representation
    |
    | key switch
    v
RLWE representation under
LWE-derived secret key
    |
    | coefficient extraction
    v
LWE ciphertext (a,b)
    |
    | modulus conversion Q' -> q
    v
final FHEW ciphertext
```

---

# 25. FHEW → CKKS

The reverse direction is different.

The high-level call is:

```text
EvalFHEWtoCKKS()
|
+-- construct A from LWE ciphertexts
|
+-- construct B from LWE ciphertexts
|
+-- EvalPartialHomDecryption()
|   |
|   +-- EvalLTRectPrecomputeSwitch()
|   |
|   +-- EvalLTRectWithPrecomputeSwitch()
|
v
CKKS ciphertext
```

The important mathematical operation here is the homomorphic evaluation of:

```text
b - a · s
```

inside CKKS.

For each LWE ciphertext:

```text
b = a · s + Δm + e  (mod q)
```

therefore:

```text
b - a · s = Δm + e  (mod q)
```

OpenFHE embeds the LWE secret-key information into the CKKS-side
computation so that the decryption expression can be evaluated
homomorphically rather than explicitly decrypting the LWE ciphertext.

---

# 26. Full Scheme-Switching Picture

```text
                         Scheme Switching
                                |
                +---------------+---------------+
                |                               |
                v                               v
          CKKS -> FHEW                    FHEW -> CKKS
                |                               |
                v                               v
   EvalCKKStoFHEW()                    EvalFHEWtoCKKS()
                |                               |
                v                               v
   Slots -> Coefficients              Build A and B
                |                               |
                v                               v
       Modulus switching               Homomorphic
                |                      decryption
                v                               |
         Key switching                         v
                |                         CKKS result
                v
       RLWE -> LWE extraction
                |
                v
          LWE ciphertext
```

---

# 27. Key Mental Model

When reading:

```cpp
EvalLTWithPrecomputeSwitch(...)
```

think:

```text
"I have an encrypted CKKS vector.

I have already converted the transformation matrix
into plaintext diagonal vectors.

Now I am evaluating:

    y = Σk Dk ⊙ Rotk(x)

using CKKS rotations, plaintext multiplication,
and ciphertext additions.

OpenFHE uses BSGS and extended-RNS machinery
to perform this efficiently."
```

When reading:

```cpp
EvalFHEWtoCKKS(...)
```

think:

```text
"I have LWE ciphertexts satisfying:

    b = a · s + Δm + e  (mod q)

and I want to evaluate:

    b - a · s

homomorphically in CKKS."
```
