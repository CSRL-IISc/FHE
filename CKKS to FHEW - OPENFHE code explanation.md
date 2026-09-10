# OpenFHE CKKS → FHEW Scheme Switching

This document covers **only the CKKS → FHEW direction** in `ckksrns-schemeswitching.cpp`.

## 1. Overview

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

The key idea is:

```text
CKKS slots
    |
    | linear transformation
    v
coefficient-oriented representation
    |
    | modulus switching + key switching
    v
RLWE representation
    |
    | coefficient extraction
    v
LWE ciphertexts
```

## 2. Why Is the Linear Transform Needed?

CKKS exposes a logical vector of slots:

```text
x = [x0, x1, x2, ..., xn-1]
```

However, the ciphertext itself is an RLWE polynomial. The logical slot
values are not simply stored as consecutive polynomial coefficients.

Therefore, before extracting LWE ciphertexts, OpenFHE applies a linear
transform that rearranges the encrypted information into the required
polynomial-coefficient representation.

Conceptually:

```text
CKKS slot representation
        |
        | linear transform
        v
coefficient-oriented representation
        |
        v
coefficient extraction
        |
        v
LWE ciphertexts
```

In this source, `EvalSlotsToCoeffsSwitch()` performs the higher-level
slots-to-coefficients operation, and `EvalLTWithPrecomputeSwitch()` is the
routine that evaluates the precomputed linear transform.

## 3. Mathematical Form of the Linear Transform

Suppose:

```text
x = [x0, x1, x2, x3]
```

and the desired transformation is:

```text
y = T · x
```

with:

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

The FHE Textbook shows that matrix-vector multiplication can be written
using cyclic diagonals, rotations, and slot-wise multiplication.

For a 4-slot example, define:

```text
D0 = [t00, t11, t22, t33]

D1 = [t01, t12, t23, t30]

D2 = [t02, t13, t20, t31]

D3 = [t03, t10, t21, t32]
```

Then the transform is:

```text
y = D0 ⊙ Rot0(x)
  + D1 ⊙ Rot1(x)
  + D2 ⊙ Rot2(x)
  + D3 ⊙ Rot3(x)
```

where:

```text
⊙       = slot-wise multiplication
Rotk(x) = cyclic slot rotation by k
```

For example:

```text
x = [x0, x1, x2, x3]

Rot0(x) = [x0, x1, x2, x3]
Rot1(x) = [x1, x2, x3, x0]
Rot2(x) = [x2, x3, x0, x1]
Rot3(x) = [x3, x0, x1, x2]
```

The exact rotation sign/indexing convention must match the convention used
by OpenFHE's `ExtractShiftedDiagonal()` and rotation routines. The example
above is only a concrete illustration of the diagonal method.

## 4. What Does `A` Mean in OpenFHE?

This is an important distinction.

### 4.1 `A` as a mathematical matrix

In:

```text
y = A · x
```

`A` usually means the original dense transformation matrix.

### 4.2 `A` in `EvalLTWithPrecomputeSwitch()`

The parameter:

```cpp
EvalLTWithPrecomputeSwitch(cc, ctxt, A, dim1)
```

is a vector of **precomputed plaintexts**.

These plaintexts represent the shifted/cyclic diagonals of the
transformation matrix.

Conceptually:

```text
Transformation matrix T
        |
        | ExtractShiftedDiagonal()
        v
D0, D1, D2, ...
        |
        | MakeAuxPlaintext()
        v
A[0], A[1], A[2], ...
```

So:

```text
A[0] -> plaintext encoding of D0
A[1] -> plaintext encoding of D1
A[2] -> plaintext encoding of D2
...
```

Therefore, in:

```cpp
A[bStep*j + i]
```

`A[...]` is **not an individual matrix element**. It is a plaintext
containing one of the transform's diagonals.

## 5. `EvalLTPrecomputeSwitch()`

### Purpose

`EvalLTPrecomputeSwitch()` prepares the diagonal plaintexts needed for the
linear transform.

Conceptually:

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
scaling / arrangement
        |
        v
MakeAuxPlaintext(...)
        |
        v
A[0], A[1], A[2], ...
```

A relevant source operation is:

```cpp
ExtractShiftedDiagonal(newA, ji)
```

This extracts a shifted diagonal from the transformation matrix.

The diagonal is then converted into a CKKS plaintext using:

```cpp
FHECKKSRNS::MakeAuxPlaintext(...)
```

So:

```text
EvalLTPrecomputeSwitch()
=
prepare plaintext diagonals for the linear transform
```

It prepares the data; it does not perform the ciphertext transformation itself.

## 6. `EvalLTWithPrecomputeSwitch()`

### Purpose

This function takes:

```text
1. the encrypted CKKS ciphertext
2. the precomputed diagonal plaintexts
```

and actually evaluates the linear transform.

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

The mathematical operation is:

```text
y = Σk Dk ⊙ Rotk(x)
```

For the 4-slot example:

```text
y = D0 ⊙ Rot0(x)
  + D1 ⊙ Rot1(x)
  + D2 ⊙ Rot2(x)
  + D3 ⊙ Rot3(x)
```

The function implements this diagonal decomposition using CKKS rotations,
plaintext multiplication, and additions.

## 7. `EvalFastRotationPrecompute()`

OpenFHE performs:

```cpp
digits = cc.EvalFastRotationPrecompute(ctxt);
```

Conceptually:

```text
prepare common information needed by multiple rotations
```

This is an optimization. It does not correspond to a new mathematical term
in:

```text
y = Σk Dk ⊙ Rotk(x)
```

## 8. `EvalFastRotationExt()`

OpenFHE performs operations such as:

```cpp
cc.EvalFastRotationExt(ctxt, i, digits, true)
```

Conceptually, this produces:

```text
Roti(x)
```

So the mapping is:

```text
FHE Textbook                 OpenFHE
------------------------------------------------
Roti(x)                      EvalFastRotationExt(...)
```

The operation is a homomorphic CKKS slot rotation.

## 9. `EvalMultExt()`

A core operation is:

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

The mapping is:

```text
fastRotation[i - 1]
        |
        v
Rotk(x)

A[bStep*j + i]
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

## 10. `EvalAddExtInPlace()`

OpenFHE accumulates the terms using:

```cpp
FHECKKSRNS::EvalAddExtInPlace(
    inner,
    FHECKKSRNS::EvalMultExt(...)
);
```

Conceptually:

```text
accumulator
    =
accumulator + Dk ⊙ Rotk(x)
```

After all diagonal terms are accumulated:

```text
y = D0 ⊙ Rot0(x)
  + D1 ⊙ Rot1(x)
  + D2 ⊙ Rot2(x)
  + ...
```

Thus:

```text
EvalAddExtInPlace()
=
add the diagonal contributions together
```

## 11. Why BSGS Is Used

A direct implementation would conceptually require many operations:

```text
D0 ⊙ Rot0(x)
D1 ⊙ Rot1(x)
D2 ⊙ Rot2(x)
...
Dn-1 ⊙ Rotn-1(x)
```

For a large number of slots, evaluating all rotations independently is
expensive.

OpenFHE therefore uses **baby-step/giant-step (BSGS)**.

The code contains:

```cpp
bStep
gStep
```

and accesses diagonals using:

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

Thus:

```text
A[bStep*j + i]
```

corresponds conceptually to:

```text
Dk
```

where:

```text
k = j*bStep + i
```

BSGS changes the implementation strategy used to generate and reuse
rotations. It does not change the mathematical result:

```text
y = Σk Dk ⊙ Rotk(x)
```

## 12. `KeySwitchExt()`

The source also contains:

```cpp
cc.KeySwitchExt(ctxt, true)
```

This should not be interpreted as an additional matrix operation.

It is part of the ciphertext/RNS representation machinery used by the
extended-RNS CKKS operations.

The mathematical transform remains:

```text
Dk ⊙ Rotk(x)
```

`KeySwitchExt()` supports the implementation of that operation.

## 13. `KeySwitchDown()`

Similarly:

```cpp
cc.KeySwitchDown(...)
```

is used to return from the extended representation used during the
calculation.

Conceptually:

```text
extended ciphertext representation
        |
        v
KeySwitchDown()
        |
        v
normal CKKS representation
```

Again, this is implementation machinery rather than a new mathematical
term in the linear transform.

## 14. `EvalSlotsToCoeffsSwitch()`

This is the higher-level slots-to-coefficients operation used in the
CKKS → FHEW path.

Its conceptual role is:

```text
CKKS slot representation
        |
        v
linear transform
        |
        v
coefficient-oriented representation
```

Internally, it calls:

```text
EvalLTWithPrecomputeSwitch()
```

to evaluate the precomputed linear transform.

So:

```text
EvalSlotsToCoeffsSwitch()
        |
        +-- prepare / adjust ciphertext
        |
        +-- EvalLTWithPrecomputeSwitch()
                |
                +-- rotations
                +-- diagonal plaintext multiplication
                +-- additions
```

This is the bridge between the logical CKKS slot representation and the
coefficient representation needed by the later LWE extraction.

## 15. `EvalCKKStoFHEW()`

This is the main runtime function for the CKKS → FHEW scheme switch.

The sequence is:

```text
EvalCKKStoFHEW()
|
+-- EvalSlotsToCoeffsSwitch()
|   |
|   +-- EvalLTWithPrecomputeSwitch()
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

Each stage performs a different task.

## 16. `ModSwitch()`

After the slot-to-coefficient linear transform, the ciphertext is moved
from the current CKKS modulus to the modulus required by the scheme switch.

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

This is modulus conversion, not matrix multiplication.

## 17. `KeySwitch()`

The code then performs:

```cpp
ccKS->KeySwitch(ctxtKS, m_CKKStoFHEWswk)
```

The purpose is to switch the ciphertext to the RLWE representation
associated with the LWE/FHEW secret key.

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

This makes the ciphertext suitable for coefficient extraction into LWE
ciphertexts.

## 18. `ExtractLWEpacked()`

This function extracts the polynomial coefficient vectors from the RLWE
ciphertext.

Conceptually:

```text
RLWE ciphertext

(A(X), B(X))

       |
       v

A = [A0, A1, A2, ...]
B = [B0, B1, B2, ...]
```

The function therefore converts the polynomial representation into a
packed coefficient representation.

It does not yet create an individual LWE ciphertext.

## 19. `ExtractLWECiphertext()`

This function selects the coefficient positions needed for one LWE
ciphertext.

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

The resulting object is an LWE ciphertext.

The source performs the required coefficient indexing, including the
negated/reversed indexing used by the extraction convention.

## 20. LWE Ciphertext Equation

The LWE ciphertext is written as:

```text
(a, b)
```

with:

```text
b = a · s + Δm + e   (mod q)
```

where:

```text
a = LWE vector
s = LWE secret key
m = plaintext message
Δ = message scaling factor
e = LWE error
q = LWE modulus
```

Therefore LWE decryption computes:

```text
b - a · s = Δm + e   (mod q)
```

This is the quantity from which the plaintext message is recovered.

## 21. `RoundqQAlter()`

The CKKS-side modulus and the final LWE modulus may differ:

```text
Q' != q
```

Therefore the extracted LWE values may need to be converted from the
CKKS-side modulus to the LWE modulus.

OpenFHE uses:

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

This is the final modulus conversion before the LWE ciphertexts are
returned.

## 22. Complete CKKS → FHEW Function Tree

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
|       |       +-- Rotk(x)
|       |
|       +-- KeySwitchExt()
|       |
|       +-- EvalMultExt()
|       |       |
|       |       +-- Dk ⊙ Rotk(x)
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
LWE / FHEW ciphertexts
```

## 23. Mathematical View of the CKKS → FHEW Path

The entire path can be viewed as:

```text
CKKS slots
    |
    | y = T · x
    |
    | y = D0 ⊙ Rot0(x)
    |   + D1 ⊙ Rot1(x)
    |   + D2 ⊙ Rot2(x)
    |   + ...
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
RLWE representation associated
with the LWE-derived secret key
    |
    | coefficient extraction
    v
LWE ciphertexts (a,b)
    |
    | modulus conversion Q' -> q
    v
final FHEW/LWE ciphertexts
```

## 24. One-Line Interpretation of Each Important Function

```text
EvalCKKStoFHEW()
    = perform the complete CKKS → FHEW conversion

EvalSlotsToCoeffsSwitch()
    = transform CKKS slots into the coefficient-oriented representation
      needed for extraction

EvalLTPrecomputeSwitch()
    = precompute plaintext diagonals of the linear transform

EvalLTWithPrecomputeSwitch()
    = apply the precomputed linear transform homomorphically

EvalFastRotationPrecompute()
    = prepare reusable data for fast rotations

EvalFastRotationExt()
    = homomorphically rotate CKKS slots

EvalMultExt()
    = multiply a rotated ciphertext by a plaintext diagonal

EvalAddExtInPlace()
    = accumulate the diagonal contributions

ModSwitch()
    = move the ciphertext from Q to Q'

KeySwitch()
    = switch to the RLWE representation associated with the
      LWE-derived key

ExtractLWEpacked()
    = extract polynomial coefficient vectors A and B

ExtractLWECiphertext()
    = select coefficients and construct an individual LWE ciphertext

RoundqQAlter()
    = convert/round values from the CKKS-side modulus Q' to
      the LWE modulus q
```

## 25. Source Mapping

Main implementation file:

```text
ckksrns-schemeswitching.cpp
```

Important functions in that file:

```text
switchingKeyGenRLWEcc()
    -> prepares the CKKS-to-FHEW key-switching key

EvalLTPrecomputeSwitch()
    -> constructs the plaintext diagonal representation

EvalLTWithPrecomputeSwitch()
    -> evaluates the linear transform

EvalSlotsToCoeffsSwitch()
    -> performs the slots-to-coefficients switching transform

ExtractLWEpacked()
    -> extracts packed RLWE A/B coefficient vectors

ExtractLWECiphertext()
    -> constructs an LWE ciphertext from selected coefficients

EvalCKKStoFHEW()
    -> complete CKKS-to-FHEW runtime path
```

Related CKKS functionality is implemented through `FHECKKSRNS` in:

```text
ckksrns-fhe.cpp
ckksrns-fhe.h
```

The scheme-switching parameters are defined in:

```text
scheme-swch-params.h
```

The FHEW/RLWE base context relationship involves:

```text
rns-fhe.h
```
