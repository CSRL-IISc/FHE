# OpenFHE CKKS → FHEW Scheme Switching

This note explains **only the CKKS → FHEW path** in OpenFHE's
`ckksrns-schemeswitching.cpp`, including the function-call tree,
important code operations, and their mathematical interpretation.

The main runtime function is:

```cpp
SWITCHCKKSRNS::EvalCKKStoFHEW(...)
```

---

# 1. Big-Picture Function Call Tree

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
|       |
|       +-- KeySwitchExt()
|       |
|       +-- FHECKKSRNS::EvalMultExt()
|       |
|       +-- FHECKKSRNS::EvalAddExtInPlace()
|       |
|       +-- KeySwitchDownFirstElement()
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
+-- return vector<LWECiphertext>
```

A useful conceptual view is:

```text
CKKS ciphertext
      |
      v
slots-to-coefficients linear transform
      |
      v
CKKS/RLWE ciphertext in the representation needed for extraction
      |
      v
modulus switch
      |
      v
key switch
      |
      v
extract RLWE polynomial coefficients
      |
      v
construct LWE ciphertexts
      |
      v
convert Q' -> q if necessary
      |
      v
FHEW/LWE ciphertexts
```

---

# 2. Why Is the Linear Transform Needed?

A CKKS ciphertext logically represents a vector of slot values:

```text
x = [x0, x1, x2, ..., xN-1]
```

But the ciphertext is internally an RLWE polynomial. The logical CKKS
slots are not simply stored as consecutive polynomial coefficients.

For the CKKS → FHEW switch, OpenFHE therefore first applies a linear
transformation that puts the desired information into a coefficient
representation from which LWE ciphertexts can be extracted.

Conceptually:

```text
CKKS slot representation
        |
        | linear transform
        v
coefficient-oriented representation
        |
        v
RLWE coefficient extraction
        |
        v
LWE ciphertexts
```

The higher-level function responsible for this stage is:

```cpp
EvalSlotsToCoeffsSwitch(...)
```

and the routine that actually evaluates the precomputed linear transform
is:

```cpp
EvalLTWithPrecomputeSwitch(...)
```

---

# 3. The Linear Transform as Matrix Multiplication

Suppose the encrypted slot vector is:

```text
x = [x0, x1, x2, x3]
```

and the desired linear transform is:

```text
T =
[ t00 t01 t02 t03
  t10 t11 t12 t13
  t20 t21 t22 t23
  t30 t31 t32 t33 ]
```

Then:

```text
y = T · x
```

means:

```text
y0 = t00*x0 + t01*x1 + t02*x2 + t03*x3

y1 = t10*x0 + t11*x1 + t12*x2 + t13*x3

y2 = t20*x0 + t21*x1 + t22*x2 + t23*x3

y3 = t30*x0 + t31*x1 + t32*x2 + t33*x3
```

The FHE Textbook's matrix-multiplication formulation rewrites this using
cyclic diagonals and rotations.

For this small example, define:

```text
D0 = [t00, t11, t22, t33]
D1 = [t01, t12, t23, t30]
D2 = [t02, t13, t20, t31]
D3 = [t03, t10, t21, t32]
```

Then, for one consistent rotation convention:

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

Then:

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

Adding them gives the matrix-vector product.

**Important:** the exact sign/direction of `Rotk` and the corresponding
shifted diagonal depends on the convention used by OpenFHE's
`ExtractShiftedDiagonal()` and rotation routines. The example above is a
mathematical illustration of the diagonal method, not a claim about the
literal sign convention in every OpenFHE helper.

---

# 4. What Does `A` Mean?

There are two meanings of `A` that should not be confused.

## 4.1 `A` as a mathematical matrix

In ordinary matrix notation:

```text
 y = A · x
```

`A` is the original dense transformation matrix.

## 4.2 `A` in `EvalLTWithPrecomputeSwitch()`

In:

```cpp
EvalLTWithPrecomputeSwitch(cc, ctxt, A, dim1)
```

`A` is a vector of **precomputed CKKS plaintexts**, not the original dense
matrix.

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

Thus:

```text
A[0] -> plaintext encoding of D0
A[1] -> plaintext encoding of D1
A[2] -> plaintext encoding of D2
...
```

Therefore, when the code contains:

```cpp
A[bStep*j + i]
```

read it as:

```text
one precomputed plaintext diagonal of the transform
```

not:

```text
one scalar element of the original matrix
```

---

# 5. Precomputation Path

The runtime call depends on linear-transform precomputation performed by:

```text
EvalCKKStoFHEWPrecompute()
        |
        +-- construct transform matrices / U matrices
        |
        +-- EvalLTPrecomputeSwitch()
                |
                +-- ExtractShiftedDiagonal()
                +-- scaling / coefficient arrangement
                +-- MakeAuxPlaintext()
                |
                v
             m_U0Pre / related plaintext vectors
```

The exact transform used by the scheme-switching code is constructed from
its internal `U0`/`U1` data and related parameters. The generic mathematical
role is still:

```text
transformation matrix
        -> cyclic diagonals
        -> CKKS plaintext diagonals
```

---

# 6. `EvalLTPrecomputeSwitch()`

## Purpose

`EvalLTPrecomputeSwitch()` prepares the plaintext representation of a
linear transformation so it can later be applied efficiently to a CKKS
ciphertext.

There are overloads for square and rectangular transformations.

The square-matrix version conceptually does:

```text
matrix T
  |
  +-- choose number of dimensions / slots
  |
  +-- determine modulus/tower representation
  |
  +-- extract shifted diagonals
  |
  +-- convert each diagonal to CKKS plaintext
  |
  v
plaintext diagonal vector
```

## Important source operations

The code performs operations of the form:

```cpp
ExtractShiftedDiagonal(newA, ji)
```

and then constructs a plaintext using:

```cpp
FHECKKSRNS::MakeAuxPlaintext(
    cc,
    elementParamsPtr,
    ...,
    1,
    towersToDrop,
    M4
)
```

The conceptual mapping is:

```text
ExtractShiftedDiagonal()
    -> obtain Dk

MakeAuxPlaintext()
    -> encode Dk as a CKKS plaintext
```

So this function is **preparation**, not ciphertext evaluation.

---

# 7. Why `MakeAuxPlaintext()`?

The diagonal vectors are ordinary vectors of values. CKKS ciphertext
operations need a CKKS plaintext object to multiply by a ciphertext.

Therefore:

```text
Dk vector
   |
   v
MakeAuxPlaintext()
   |
   v
CKKS plaintext containing Dk
```

That plaintext can then participate in:

```text
ciphertext × plaintext
```

inside `EvalLTWithPrecomputeSwitch()`.

---

# 8. `EvalLTWithPrecomputeSwitch()`

## Purpose

This is the **execution stage** of the linear transform.

Input:

```text
ctxt = encrypted CKKS vector
A    = precomputed plaintext diagonals
```

Output:

```text
encrypted transformed vector
```

Conceptually:

```text
                     ctxt = x
                         |
             +-----------+-----------+
             |           |           |
             v           v           v
           Rot0         Rot1        Rot2 ...
             |           |           |
             v           v           v
           × D0        × D1         × D2
             |           |           |
             +-----------+-----------+
                         |
                         v
                       sum
                         |
                         v
                         y
```

Mathematically:

```text
y = Σk Dk ⊙ Rotk(x)
```

The important distinction is:

```text
EvalLTPrecomputeSwitch()
    = build D0, D1, D2, ... plaintexts

EvalLTWithPrecomputeSwitch()
    = apply those plaintexts to the encrypted ciphertext
```

---

# 9. The First Important Code Operation: Fast-Rotation Precomputation

Inside `EvalLTWithPrecomputeSwitch()` the code first obtains data needed
for multiple rotations:

```cpp
auto digits = cc.EvalFastRotationPrecompute(ctxt);
```

Conceptually:

```text
ciphertext x
    |
    v
precompute common rotation/key-switch information
    |
    v
digits
```

This is an implementation optimization.

It does **not** correspond to an extra term in the mathematical equation:

```text
y = Σk Dk ⊙ Rotk(x)
```

---

# 10. `EvalFastRotationExt()`

The source then evaluates fast rotations such as:

```cpp
cc.EvalFastRotationExt(ctxt, j, digits, true)
```

Conceptually:

```text
EvalFastRotationExt(ctxt, j, ...)
             |
             v
         Rotj(x)
```

The mapping is:

```text
FHE Textbook                 OpenFHE
------------------------------------------------
Rotj(x)                      EvalFastRotationExt(...)
```

This is a **homomorphic slot rotation**. The plaintext is never
recovered from the ciphertext to perform this operation.

---

# 11. `KeySwitchExt()` Inside the Linear Transform

The code contains operations such as:

```cpp
cc.KeySwitchExt(ctxt, true)
```

This is not another mathematical matrix term.

It is part of the extended-RNS / ciphertext representation machinery
needed for OpenFHE's CKKS operations.

At the mathematical level, the operation we care about remains:

```text
Dk ⊙ Rotk(x)
```

`KeySwitchExt()` helps put the ciphertext into the representation required
by the subsequent extended multiplication/rotation operations.

---

# 12. `EvalMultExt()`

A core operation is of the form:

```cpp
inner = FHECKKSRNS::EvalMultExt(
    cc.KeySwitchExt(ctxt, true),
    A[bStep*j]
);
```

or, for baby-step terms:

```cpp
FHECKKSRNS::EvalMultExt(
    fastRotation[i - 1],
    A[bStep*j + i]
)
```

Conceptually:

```text
fastRotation[i-1]
        |
        v
     Rotk(x)

A[bStep*j+i]
        |
        v
       Dk

EvalMultExt()
        |
        v
Dk ⊙ Rotk(x)
```

This is the slot-wise multiplication part of the diagonalized linear
transform.

---

# 13. `EvalAddExtInPlace()`

The resulting diagonal contributions are accumulated using:

```cpp
FHECKKSRNS::EvalAddExtInPlace(
    inner,
    FHECKKSRNS::EvalMultExt(...)
);
```

Conceptually:

```text
inner = inner + Dk ⊙ Rotk(x)
```

After enough iterations:

```text
y = D0 ⊙ Rot0(x)
  + D1 ⊙ Rot1(x)
  + D2 ⊙ Rot2(x)
  + ...
```

Thus:

```text
EvalAddExtInPlace()
    = accumulate the diagonal contributions
```

---

# 14. How the Code Organizes the Terms: BSGS

A direct implementation would conceptually evaluate:

```text
D0 ⊙ Rot0(x)
D1 ⊙ Rot1(x)
D2 ⊙ Rot2(x)
...
DN-1 ⊙ RotN-1(x)
```

That can require many expensive rotations.

OpenFHE uses a **baby-step/giant-step (BSGS)** organization.

The function computes quantities such as:

```cpp
uint32_t bStep = dim1;
uint32_t gStep = ceil(slots / bStep);
```

and uses indices such as:

```cpp
A[bStep*j + i]
```

The conceptual rotation index is:

```text
k = j*bStep + i
```

where:

```text
j = giant-step index
i = baby-step index
```

Therefore:

```text
A[bStep*j + i]
```

means approximately:

```text
Dk
```

for:

```text
k = j*bStep + i
```

BSGS changes **how the terms are evaluated efficiently**. It does not
change the mathematical linear transformation.

---

# 15. What BSGS Looks Like in a Small Example

Suppose:

```text
number of diagonal terms = 8
bStep = 2
```

Then:

```text
k = j*2 + i
```

The terms are grouped as:

```text
j = 0:  k=0, k=1
j = 1:  k=2, k=3
j = 2:  k=4, k=5
j = 3:  k=6, k=7
```

So the plaintext-diagonal indices look like:

```text
A[0], A[1]
A[2], A[3]
A[4], A[5]
A[6], A[7]
```

The implementation can reuse common rotation information instead of
independently performing every possible rotation from scratch.

---

# 16. `KeySwitchDownFirstElement()` and `KeySwitchDown()`

The extended operations used inside the linear transform need to be
brought back to the normal representation at appropriate points.

The code uses operations such as:

```cpp
cc.KeySwitchDownFirstElement(...)
```

and:

```cpp
cc.KeySwitchDown(...)
```

Conceptually:

```text
extended-RNS representation
        |
        v
KeySwitchDown / KeySwitchDownFirstElement
        |
        v
normal CKKS representation
```

These operations are part of the implementation-level representation and
key-management machinery. They should not be interpreted as additional
terms in:

```text
y = Σk Dk ⊙ Rotk(x)
```

---

# 17. `EvalSlotsToCoeffsSwitch()`

`EvalSlotsToCoeffsSwitch()` is the higher-level wrapper used by the
scheme-switching path.

Conceptually:

```text
CKKS slots
    |
    v
slots-to-coefficients linear transform
    |
    v
coefficient-oriented CKKS representation
```

The function also contains ciphertext preparation related to the CKKS
parameters. For example, the source checks whether the required
precomputation exists and handles modulus/scaling details for flexible
scaling modes.

The critical nested call is:

```cpp
EvalLTWithPrecomputeSwitch(...)
```

So the relationship is:

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

---

# 18. `m_U0Pre`

`EvalSlotsToCoeffsSwitch()` uses precomputed transform data such as:

```text
m_U0Pre
```

The important point is that this is not the ciphertext itself.

It is precomputed information representing part of the linear transform.

Conceptually:

```text
U0 transformation data
        |
        | EvalLTPrecomputeSwitch()
        v
m_U0Pre
        |
        | EvalLTWithPrecomputeSwitch()
        v
transformed ciphertext
```

If the required precomputation is missing, the switch function cannot
perform the linear transform and reports a precomputation error.

---

# 19. Why This Is Related to SlotToCoeff

The FHE Textbook describes `SlotToCoeff` as a linear transformation that
moves information from CKKS slots into polynomial-coefficient positions.

The OpenFHE scheme-switching helper is implementing the same **kind of
linear-transform mechanism**, using matrix diagonals, rotations,
plaintext multiplication, and addition.

Thus the conceptual mapping is:

```text
FHE Textbook
------------------------------------------------
Slot vector
    |
SlotToCoeff linear transform
    |
coefficient-oriented representation

OpenFHE
------------------------------------------------
CKKS slots
    |
EvalSlotsToCoeffsSwitch()
    |
EvalLTWithPrecomputeSwitch()
    |
coefficient-oriented representation
```

The exact matrices (`U0`, `U1`, etc.) are OpenFHE implementation details;
the textbook supplies the mathematical framework for evaluating a linear
transform homomorphically.

---

# 20. Main Runtime Function: `EvalCKKStoFHEW()`

The main function is conceptually:

```text
EvalCKKStoFHEW()
|
+-- slots -> coefficients
|
+-- modulus switch
|
+-- key switch
|
+-- extract RLWE coefficients
|
+-- form LWE ciphertexts
|
+-- convert modulus Q' -> q if required
|
+v
LWE/FHEW ciphertexts
```

The source first limits the number of requested ciphertexts to the CKKS
slot capacity, then performs the linear transform and extraction steps.

A simplified representation of the code is:

```cpp
ctxtDecoded = EvalSlotsToCoeffsSwitch(*ccCKKS, ciphertext);

ModReduceInternalInPlace(ctxtDecoded, 1);

ctxtKS = m_ctxtKS->Clone();
ModSwitch(ctxtDecoded, ctxtKS, m_modulus_CKKS_from);

ctSwitched = ccKS->KeySwitch(ctxtKS, m_CKKStoFHEWswk);

AandB = ExtractLWEpacked(ctSwitched);

for (...) {
    auto lwe = ExtractLWECiphertext(
        AandB,
        m_modulus_CKKS_from,
        n,
        index
    );

    // if needed, convert Q' -> q using RoundqQAlter()
}
```

The actual source contains additional bookkeeping and parameter handling,
but these are the essential operations.

---

# 21. `ModReduceInternalInPlace()`

Before the explicit scheme-switch modulus conversion, the code performs:

```cpp
ModReduceInternalInPlace(ctxtDecoded, 1);
```

Conceptually, this removes one level/tower according to the CKKS modulus
chain and adjusts the ciphertext representation accordingly.

This is distinct from the later explicit mapping to
`m_modulus_CKKS_from`.

The relevant distinction is:

```text
ModReduceInternalInPlace()
    = CKKS modulus-chain level reduction

ModSwitch(..., m_modulus_CKKS_from)
    = move the representation to the specific modulus
      required by scheme switching
```

---

# 22. `ModSwitch()`

The code then performs the scheme-switch modulus conversion:

```cpp
ModSwitch(ctxtDecoded, ctxtKS, m_modulus_CKKS_from);
```

Conceptually:

```text
CKKS ciphertext at Q
        |
        v
     ModSwitch
        |
        v
CKKS ciphertext at Q'
```

Here:

```text
Q  = current CKKS-side modulus
Q' = modulus chosen for the CKKS -> FHEW switch
```

This operation changes the modulus representation. It is not the matrix
linear transform.

---

# 23. `switchingKeyGenRLWEcc()`

Before the runtime `KeySwitch()` can happen, the appropriate switching
key has to be generated.

The relevant helper is:

```cpp
switchingKeyGenRLWEcc(
    ckksSKto,
    ckksSKfrom,
    LWEsk
)
```

Conceptually:

```text
CKKS secret key
       |
       | construct RLWE representation
       | associated with LWE secret key
       v
transformed secret-key representation
       |
       v
KeySwitchGen(...)
       |
       v
CKKS -> FHEW key-switching key
```

The source maps the secret-key coefficients into a representation that
corresponds to the LWE secret key and then invokes the CKKS key-switch-key
generation mechanism.

This is important because the post-transform RLWE ciphertext must be
associated with the secret-key representation from which the LWE
ciphertext will eventually be extracted.

---

# 24. `KeySwitch()`

At runtime:

```cpp
ccKS->KeySwitch(ctxtKS, m_CKKStoFHEWswk)
```

is performed.

Conceptually:

```text
RLWE ciphertext under CKKS key
            |
            v
        KeySwitch
            |
            v
RLWE ciphertext under the RLWE representation
corresponding to the LWE/FHEW secret key
```

This is a **key representation transformation**.

It is not the LWE decryption operation and not a plaintext extraction.

---

# 25. `ExtractLWEpacked()`

After key switching, the result is still an RLWE ciphertext.

The helper:

```cpp
ExtractLWEpacked(ctSwitched)
```

takes the polynomial components and puts them into explicit coefficient
vectors.

Conceptually, an RLWE ciphertext:

```text
(A(X), B(X))
```

becomes:

```text
A = [A0, A1, A2, ..., A(N-1)]
B = [B0, B1, B2, ..., B(N-1)]
```

In the source, the function accesses the first relevant element, switches
it to coefficient format, and obtains the underlying values.

The conceptual mapping is:

```text
RLWE polynomial A(X)  -> coefficient vector A
RLWE polynomial B(X)  -> coefficient vector B
```

---

# 26. `ExtractLWECiphertext()`

`ExtractLWEpacked()` gives a packed set of coefficients. The next helper
chooses the positions belonging to one LWE ciphertext:

```cpp
ExtractLWECiphertext(
    AandB,
    modulus,
    n,
    index
)
```

The resulting LWE ciphertext is:

```text
(a, b)
```

where:

```text
a = [a0, a1, ..., a(n-1)]
```

and `b` is the selected coefficient from the packed B vector.

The source performs a specific reversed/negated coefficient indexing when
building `a`. This is part of the RLWE-to-LWE extraction convention used
by the implementation.

Conceptually:

```text
packed RLWE coefficients
        |
        v
select coefficient positions
        |
        +----> a = [a0, ..., a(n-1)]
        |
        +----> b = selected B[index]
        |
        v
LWE ciphertext (a,b)
```

---

# 27. Why Can an RLWE Ciphertext Become an LWE Ciphertext by Extraction?

The coefficient-extraction step works because, after the preceding linear
transform and key-switching operations, the relevant polynomial
coefficients encode the required LWE samples.

So the sequence is not:

```text
arbitrary RLWE
    -> arbitrary coefficient
    -> LWE
```

Instead it is:

```text
CKKS slots
    -> carefully chosen linear transform
    -> required coefficient positions
    -> key switch to the LWE-derived key
    -> extract those coefficients
    -> LWE ciphertexts
```

This is why the slots-to-coefficients stage is essential to the complete
CKKS → FHEW conversion.

---

# 28. LWE Ciphertext Equation

Once a ciphertext has been extracted, it has the usual LWE form:

```text
(a, b)
```

with:

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

Therefore the LWE decryption expression is:

```text
b - a · s = Δm + e   (mod q)
```

The expression `b - a · s` should therefore be read as the LWE
**decryption expression**, not as a step performed during
`EvalCKKStoFHEW()` itself. The CKKS → FHEW routine constructs the LWE
ciphertext; FHEW operations such as sign evaluation then operate on that
LWE ciphertext.

---

# 29. `RoundqQAlter()`

The modulus used on the CKKS-side of the switch may differ from the final
LWE modulus.

Let:

```text
Q' = CKKS-side modulus used before extraction
q  = final LWE modulus
```

If:

```text
Q' != q
```

the extracted values are converted using `RoundqQAlter()`.

Conceptually:

```text
value represented modulo Q'
        |
        v
scale from Q' to q
        |
        v
round
        |
        v
value represented modulo q
```

The code performs this conversion for both the `a` entries and the `b`
value when the moduli differ.

---

# 30. Complete Detailed Call Tree

```text
SWITCHCKKSRNS::EvalCKKStoFHEW()
|
+-- determine numCtxts
|
+-- ccCKKS = ciphertext->GetCryptoContext()
|
+-- EvalSlotsToCoeffsSwitch(*ccCKKS, ciphertext)
|   |
|   +-- check m_U0Pre / precomputation
|   |
|   +-- prepare ciphertext / scaling / modulus representation
|   |
|   +-- EvalLTWithPrecomputeSwitch(...)
|       |
|       +-- determine slots, bStep, gStep
|       |
|       +-- EvalFastRotationPrecompute(ctxt)
|       |
|       +-- EvalFastRotationExt(...)
|       |       -> homomorphic slot rotation
|       |
|       +-- KeySwitchExt(...)
|       |       -> extended ciphertext representation
|       |
|       +-- EvalMultExt(...)
|       |       -> Dk ⊙ Rotk(x)
|       |
|       +-- EvalAddExtInPlace(...)
|       |       -> accumulate diagonal terms
|       |
|       +-- KeySwitchDownFirstElement(...)
|       |
|       +-- KeySwitchDown(...)
|               -> normal CKKS representation
|
+-- ModReduceInternalInPlace(ctxtDecoded, 1)
|       -> reduce CKKS modulus-chain level
|
+-- Clone m_ctxtKS
|
+-- ModSwitch(ctxtDecoded, ctxtKS, m_modulus_CKKS_from)
|       -> convert to scheme-switch modulus Q'
|
+-- ccKS->KeySwitch(ctxtKS, m_CKKStoFHEWswk)
|       -> switch RLWE key representation
|
+-- ExtractLWEpacked(ctSwitched)
|       |
|       +-- extract coefficient vector A
|       +-- extract coefficient vector B
|
+-- for each requested output index
|   |
|   +-- ExtractLWECiphertext(AandB, Q', n, index)
|   |       |
|   |       +-- construct a vector
|   |       +-- select b coefficient
|   |       +-- return (a,b)
|   |
|   +-- if Q' != q
|           |
|           +-- RoundqQAlter(a_j, q, Q')
|           +-- RoundqQAlter(b, q, Q')
|           +-- construct corrected LWE ciphertext
|
+-- return vector<LWECiphertext>
```

---

# 31. Mathematical View of the Same Tree

```text
              CKKS slots
                  |
                  | x
                  v
        linear transformation T
                  |
                  | y = T · x
                  |
                  | y = Σk Dk ⊙ Rotk(x)
                  v
       coefficient-oriented CKKS
                  |
                  | modulus reduction/switch
                  v
                Q'
                  |
                  | key switch
                  v
       RLWE under LWE-derived key
                  |
                  | coefficient extraction
                  v
             (a, b) LWE
                  |
                  | Q' -> q if necessary
                  v
           FHEW/LWE ciphertext
```

---

# 32. Exact Conceptual Mapping of the Important Operations

```text
Mathematical operation                 OpenFHE operation
================================================================
Extract a matrix diagonal              ExtractShiftedDiagonal()
Encode diagonal as CKKS plaintext      MakeAuxPlaintext()
Prepare repeated rotation data         EvalFastRotationPrecompute()
Homomorphic slot rotation              EvalFastRotationExt()
Ciphertext representation extension    KeySwitchExt()
Dk ⊙ Rotk(x)                            EvalMultExt()
Sum diagonal contributions             EvalAddExtInPlace()
Return from extended representation    KeySwitchDown*()
Apply the complete linear transform    EvalLTWithPrecomputeSwitch()
Slots -> coefficient-oriented data     EvalSlotsToCoeffsSwitch()
CKKS modulus-chain reduction           ModReduceInternalInPlace()
Convert to scheme-switch modulus       ModSwitch()
RLWE key representation change         KeySwitch()
Extract RLWE coefficients              ExtractLWEpacked()
Create individual LWE ciphertext       ExtractLWECiphertext()
Q' -> q conversion                     RoundqQAlter()
```

---

# 33. The Most Important Distinction: Precompute vs Evaluation

This distinction is useful when reading the source.

```text
PRECOMPUTATION

EvalLTPrecomputeSwitch()
        |
        +-- ExtractShiftedDiagonal()
        +-- MakeAuxPlaintext()
        |
        v
A[0], A[1], A[2], ...
```

versus:

```text
EVALUATION

EvalLTWithPrecomputeSwitch()
        |
        +-- rotate ciphertext
        +-- multiply by A[k]
        +-- add results
        |
        v
transformed ciphertext
```

So if you see:

```cpp
A[bStep*j + i]
```

you should think:

```text
precomputed diagonal plaintext
```

If you see:

```cpp
EvalFastRotationExt(...)
```

you should think:

```text
homomorphic slot rotation
```

If you see:

```cpp
EvalMultExt(...)
```

you should think:

```text
plaintext diagonal × rotated ciphertext
```

If you see:

```cpp
EvalAddExtInPlace(...)
```

you should think:

```text
accumulate the matrix-transform terms
```

---

# 34. The Core Equation to Keep in Mind

For the linear-transform portion, the clean mathematical model is:

```text
y = Σk Dk ⊙ Rotk(x)
```

where:

```text
x   = encrypted CKKS slot vector
Dk  = plaintext encoding of the k-th shifted diagonal
Rotk(x) = homomorphic CKKS slot rotation
⊙   = slot-wise multiplication
```

This realizes:

```text
y = T · x
```

for the transformation matrix `T` represented by those diagonals.

The OpenFHE implementation then wraps this mathematical operation with
BSGS optimization and extended-RNS/key-switching machinery.

---

# 35. One-Line Explanation of Each Function

```text
switchingKeyGenRLWEcc()
    = generate the CKKS -> FHEW key-switching key

EvalLTPrecomputeSwitch()
    = convert a linear-transform matrix into plaintext diagonals

EvalLTWithPrecomputeSwitch()
    = evaluate the linear transform on an encrypted CKKS ciphertext

EvalSlotsToCoeffsSwitch()
    = perform the slots-to-coefficients switching transform

EvalFastRotationPrecompute()
    = precompute data reused by multiple fast rotations

EvalFastRotationExt()
    = homomorphically rotate CKKS slots

KeySwitchExt()
    = prepare/use an extended ciphertext representation

EvalMultExt()
    = multiply the rotated ciphertext by a plaintext diagonal

EvalAddExtInPlace()
    = add the current diagonal contribution to the accumulator

KeySwitchDownFirstElement()
    = return the first extended result to the normal representation

KeySwitchDown()
    = return an extended result to the normal CKKS representation

ModReduceInternalInPlace()
    = reduce the CKKS modulus-chain level

ModSwitch()
    = move to the modulus selected for scheme switching

KeySwitch()
    = switch the RLWE ciphertext to the LWE-derived key representation

ExtractLWEpacked()
    = extract RLWE polynomial coefficients into packed A/B vectors

ExtractLWECiphertext()
    = select the coefficient positions forming one LWE ciphertext

RoundqQAlter()
    = map values from the CKKS-side modulus Q' to the LWE modulus q

EvalCKKStoFHEW()
    = execute the complete CKKS -> FHEW conversion
```

---

# 36. Final Mental Model

When reading the code, keep this picture in mind:

```text
                CKKS ciphertext
                       |
                       v
             EvalSlotsToCoeffsSwitch()
                       |
                       v
            EvalLTWithPrecomputeSwitch()
                       |
             +---------+---------+
             |         |         |
           rotate    rotate    rotate
             |         |         |
            ×D0       ×D1       ×D2 ...
             |         |         |
             +---------+---------+
                       |
                       v
                      sum
                       |
                       v
            transformed CKKS/RLWE
                       |
                       v
                 modulus switch
                       |
                       v
                   key switch
                       |
                       v
             RLWE coefficient data
                       |
                       v
              ExtractLWEpacked()
                       |
                       v
           ExtractLWECiphertext()
                       |
                       v
                  Q' -> q
                       |
                       v
                  LWE / FHEW
```

The core reason for the seemingly complicated code is therefore:

```text
1. Re-express the CKKS slot transformation as diagonal operations.
2. Apply those operations homomorphically using rotations,
   plaintext multiplication, and addition.
3. Put the result at the coefficient positions needed for extraction.
4. Switch the RLWE key representation.
5. Extract those coefficients as LWE ciphertexts.
6. Convert the modulus to the final LWE modulus when required.
```
