# OpenFHE CKKS ↔ TFHE Scheme Switching

## 1. Scope

This note explains the function-call flow around `ckksrns-schemeswitching.cpp`, with emphasis on:

- `EvalCKKStoFHEW()`
- `EvalFHEWtoCKKS()`
- the CKKS linear transforms used during switching
- LWE extraction from the RLWE ciphertext
- the caller `EvalMinSchemeSwitching()`

The descriptions below are based on the supplied OpenFHE source files and follow the terminology used there.

---

# 2. Main Runtime Call Tree: CKKS → FHEW

```text
SWITCHCKKSRNS::EvalCKKStoFHEW()
│
├── 1. EvalSlotsToCoeffsSwitch()
│      │
│      └── EvalLTWithPrecomputeSwitch()
│             │
│             ├── cc.EvalFastRotationPrecompute()
│             ├── cc.EvalFastRotationExt()
│             ├── cc.KeySwitchExt()
│             ├── FHECKKSRNS::EvalMultExt()
│             ├── FHECKKSRNS::EvalAddExtInPlace()
│             ├── cc.KeySwitchDownFirstElement()
│             ├── cc.KeySwitchDown()
│             └── cc.EvalFastRotationExt()
│
├── 2. ModReduceInternalInPlace()
│
├── 3. ModSwitch()
│
├── 4. ccKS->KeySwitch()
│
├── 5. ExtractLWEpacked()
│      ├── extract RLWE A polynomial coefficients
│      └── extract RLWE B polynomial coefficients
│
├── 6. ExtractLWECiphertext()
│      └── construct LWECiphertextImpl(a, b)
│
└── 7. RoundqQAlter()   [only when Q' != q]
       ├── round each LWE a coefficient
       └── round LWE b
```

The conceptual flow is:

```text
CKKS ciphertext
      │
      ▼
Slots → coefficients
      │
      ▼
Modulus reduce / switch
      │
      ▼
Key-switch CKKS secret-key domain
into an RLWE representation of the LWE/FHEW key
      │
      ▼
Packed RLWE polynomial
      │
      ▼
Extract coefficient slices
      │
      ▼
LWE ciphertexts
      │
      ▼
Optional modulus conversion Q' → q
      │
      ▼
FHEW / BinFHE ciphertexts
```

---

# 3. Step-by-Step: `EvalCKKStoFHEW()`

## 3.1 `SWITCHCKKSRNS::EvalCKKStoFHEW()`

**Role:** This is the top-level runtime function that converts selected CKKS slots into a vector of FHEW/LWE ciphertexts.

The implementation first limits the requested number of ciphertexts to the CKKS slot count. It obtains the CKKS crypto context and then performs a slots-to-coefficients transform.

Relevant source behavior:

```cpp
ctxtDecoded = EvalSlotsToCoeffsSwitch(*ccCKKS, ciphertext);
ModReduceInternalInPlace(ctxtDecoded, 1);
ctxtKS = m_ctxtKS->Clone();
ModSwitch(ctxtDecoded, ctxtKS, m_modulus_CKKS_from);
ctSwitched = ccKS->KeySwitch(ctxtKS, m_CKKStoFHEWswk);
```

So this function is mainly an **orchestrator**. It does not implement the linear transform, RLWE coefficient extraction, or LWE construction itself; it calls specialized helpers for those operations.

**File:** `ckksrns-schemeswitching.cpp`

---

# 4. Step 1 — `EvalSlotsToCoeffsSwitch()`

## 4.1 What it does

```text
EvalCKKStoFHEW()
      │
      └── EvalSlotsToCoeffsSwitch()
```

CKKS naturally represents packed values in polynomial slots. For scheme switching, OpenFHE needs a representation whose relevant information can be arranged into polynomial coefficients suitable for extraction into LWE ciphertexts.

`EvalSlotsToCoeffsSwitch()` therefore performs the CKKS linear transformation:

```text
slot representation
        ↓
coefficient representation
```

It uses precomputed transform data such as `m_U0Pre`.

The function checks that the required precomputation exists. It also determines the cyclotomic order and checks whether the ring is sparse relative to the expected transform size. For flexible scaling techniques, it can compress the ciphertext to fewer towers and adjust the scaling factor before continuing.

---

# 5. `EvalSlotsToCoeffsSwitch()` → `EvalLTWithPrecomputeSwitch()`

The actual linear transform is evaluated through the precomputed diagonal representation.

```text
EvalSlotsToCoeffsSwitch()
        │
        └── EvalLTWithPrecomputeSwitch()
```

The mathematical structure is approximately:


y = Σᵢ D_i Rotᵢ(x)


where:

- `x` = input CKKS ciphertext
- `Rot_i(x)` = rotated ciphertext
- `D_i` = plaintext diagonal corresponding to one rotation
- `y` = transformed CKKS ciphertext

The implementation uses baby-step/giant-step evaluation so that the number of expensive rotations/key-switch operations is reduced.

---

# 6. `EvalLTWithPrecomputeSwitch()`

## 6.1 What it does

This function evaluates a precomputed CKKS linear transform.

Its major stages are:

```text
EvalLTWithPrecomputeSwitch()
│
├── determine baby-step / giant-step dimensions
│
├── EvalFastRotationPrecompute()
│
├── generate required rotated ciphertexts
│
├── multiply rotations by precomputed plaintext diagonals
│
├── add partial results
│
└── key-switch / combine results
```

The implementation computes quantities such as:

```cpp
slots  = ...
bStep  = dim1
gStep  = ceil(slots / bStep)
digits = cc.EvalFastRotationPrecompute(ctxt);
```

The important idea is that the transform is decomposed into **giant steps** and **baby steps**.

---

## 6.2 `EvalFastRotationPrecompute()`

**Purpose:** Prepare the information needed to evaluate many rotations of the same CKKS ciphertext efficiently.

Instead of independently performing the full work for every rotation, OpenFHE performs a common precomputation once and reuses it.

Conceptually:

```text
ciphertext
   │
   └── hoisted / shared rotation data
             │
             ├── rotation 1
             ├── rotation 2
             ├── rotation 3
             └── ...
```

This is particularly useful inside baby-step/giant-step linear transforms, where many related rotations are needed.

---

## 6.3 `EvalFastRotationExt()`

**Purpose:** Produce a rotated version of the CKKS ciphertext using the precomputed rotation information.

The `Ext` form operates with the extended ciphertext representation used by the RNS CKKS implementation.

In the linear transform it is used to obtain the rotations corresponding to the baby steps and giant steps.

Conceptually:

```text
precomputed rotation data
          + rotation index
                 │
                 ▼
        rotated CKKS ciphertext
```

---

## 6.4 `KeySwitchExt()`

**Purpose:** Convert the extended representation into the form required for the extended CKKS multiplication/evaluation path.

In this code it appears in expressions such as:

```cpp
FHECKKSRNS::EvalMultExt(cc.KeySwitchExt(ctxt, true), A[...])
```

So its job is part of the internal RNS/extended ciphertext path used before multiplying by a plaintext diagonal.

---

## 6.5 `FHECKKSRNS::EvalMultExt()`

**Purpose:** Multiply an extended CKKS ciphertext by a plaintext in the extended-RNS evaluation path.

In the linear transform:

```cpp
FHECKKSRNS::EvalMultExt(
    cc.KeySwitchExt(ctxt, true),
    A[index]
)
```

Here:

- the first argument is the rotated/extended ciphertext
- the second argument is a precomputed plaintext diagonal

Mathematically:


D_i · Rotᵢ(x)


The result is one contribution to the final linear transform.

---

## 6.6 `FHECKKSRNS::EvalAddExtInPlace()`

**Purpose:** Add another extended ciphertext contribution into the current accumulator.

The code uses it in the form:

```cpp
FHECKKSRNS::EvalAddExtInPlace(
    inner,
    FHECKKSRNS::EvalMultExt(fastRotation[i - 1], A[...])
);
```

Conceptually:

```text
inner
  + D_i · Rot_i(x)
  ─────────────────
       inner
```

Repeated additions build the output of the diagonal linear transform.

---

## 6.7 `KeySwitchDownFirstElement()`

**Purpose:** Bring the first partial result down from the extended representation into the normal CKKS ciphertext representation.

The implementation treats the first giant-step result specially because it can be used as the initial accumulator.

Conceptually:

```text
extended partial result
         │
         ▼
normal CKKS ciphertext
```

---

## 6.8 `KeySwitchDown()`

**Purpose:** Convert an extended CKKS ciphertext back into the standard CKKS ciphertext form after the extended operations are finished.

It is used both while combining giant-step results and at the end of the transform.

So the extended path looks approximately like:

```text
normal CKKS
    ↓
extended representation
    ↓
rotations + plaintext multiplies + adds
    ↓
KeySwitchDown()
    ↓
normal CKKS
```

---

# 7. Step 2 — `ModReduceInternalInPlace()`

After the slot-to-coefficient transform, the code performs:

```cpp
ModReduceInternalInPlace(ctxtDecoded, 1);
```

**Purpose:** Reduce the modulus level of the CKKS ciphertext by removing one level according to the CKKS RNS representation.

This is part of preparing the ciphertext modulus for the subsequent scheme-switching steps.

It is important to distinguish this from the later explicit conversion of the CKKS modulus to the FHEW-compatible modulus. This operation changes the current CKKS modulus level; it is not yet the final CKKS-to-LWE modulus conversion.

---

# 8. Step 3 — `ModSwitch()`

The code then creates a clone of a precomputed ciphertext and performs:

```cpp
ctxtKS = m_ctxtKS->Clone();
ModSwitch(ctxtDecoded, ctxtKS, m_modulus_CKKS_from);
```

**Purpose:** Convert the ciphertext representation so that its modulus matches the modulus expected by the key-switching/extraction stage.

At this point it is useful to name the moduli:

```text
Q  = original / current CKKS modulus
Q' = modulus used as the CKKS → FHEW intermediate modulus
q  = final LWE/FHEW modulus
```

The exact implementation uses `m_modulus_CKKS_from` for the relevant CKKS-side modulus value.

The later `RoundqQAlter()` operation handles the final `Q' → q` conversion when these moduli differ.

---

# 9. Step 4 — `ccKS->KeySwitch()`

The key switching step is:

```cpp
ccKS->KeySwitch(ctxtKS, m_CKKStoFHEWswk);
```

**Purpose:** Change the ciphertext from the original CKKS secret-key domain to the key domain that corresponds to the LWE/FHEW secret key.

This is the crucial cryptographic bridge between the CKKS ciphertext and the LWE representation.

The switching key `m_CKKStoFHEWswk` is generated so that the resulting RLWE ciphertext is related to the FHEW/LWE secret key rather than only to the original CKKS secret key.

Conceptually:

```text
CKKS ciphertext under s_CKKS
             │
             │ key switch
             ▼
RLWE ciphertext under RLWE representation of s_LWE
```

---

# 10. How `m_CKKStoFHEWswk` Is Generated

The helper involved is:

```text
switchingKeyGenRLWEcc()
```

Call structure:

```text
Key generation
│
└── switchingKeyGenRLWEcc()
       │
       ├── read CKKS secret key
       ├── read LWE secret key
       ├── construct transformed secret-key representations
       └── ccCKKSto->KeySwitchGen(oldTransformedSK, RLWELWEsk)
```

### `switchingKeyGenRLWEcc()`

This helper maps the LWE secret-key information into an RLWE-compatible polynomial representation and generates the CKKS key-switching key that maps from the CKKS secret-key representation to it.

The implementation converts the relevant coefficients into values represented by `0`, `1`, or `modulus - 1`, depending on the source CKKS secret-key coefficient and the corresponding LWE secret-key element.

The final call is effectively:

```cpp
ccCKKSto->KeySwitchGen(oldTranformedSK, RLWELWEsk);
```

The result is the switching key later used by:

```cpp
ccKS->KeySwitch(ctxtKS, m_CKKStoFHEWswk);
```

---

# 11. Step 5 — `ExtractLWEpacked()`

After key switching, the ciphertext is still an RLWE-style ciphertext.

The next helper is:

```text
ExtractLWEpacked()
```

Call structure:

```text
RLWE ciphertext
     │
     └── ExtractLWEpacked()
            ├── extract A polynomial coefficients
            └── extract B polynomial coefficients
```

The source takes the first tower/element, converts it to coefficient format, and copies the coefficients of the two RLWE components into vectors.

The conceptual transformation is:

```text
RLWE ciphertext
   = (A(X), B(X))
          │
          ▼
A coefficient vector + B coefficient vector
```

The function does **not** yet create individual LWE ciphertext objects. It only exposes the packed polynomial coefficients.

---

# 12. Step 6 — `ExtractLWECiphertext()`

Next, `EvalCKKStoFHEW()` repeatedly calls:

```text
ExtractLWECiphertext()
```

Call structure:

```text
packed A/B coefficient vectors
          │
          ▼
ExtractLWECiphertext(index)
          │
          ├── select n coefficients from A
          ├── rearrange / negate according to the extraction convention
          └── use B[index] as the LWE b value
                  │
                  ▼
           LWECiphertextImpl(a, b)
```

The function allocates:

```cpp
NativeVector a(n, modulus);
```

and fills it from selected positions in the packed RLWE `A` polynomial.

It then constructs:

```cpp
LWECiphertextImpl(a, aANDb[0][index])
```

Thus the packed RLWE result contains enough coefficient information to recover multiple LWE ciphertexts.

---

# 13. Why `gap` Exists

The code computes:

```cpp
gap = ccKS->GetRingDimension() / (2 * m_numSlotsCKKS);
```

and advances the extraction index by this gap.

Conceptually:

```text
RLWE polynomial coefficients

|--- slot 0 ---| gap |--- slot 1 ---| gap |--- slot 2 ---| ...
                       ↑
                 extraction spacing
```

So the LWE ciphertexts are not simply formed from consecutive polynomial coefficients. They are packed at regularly spaced positions determined by the ring dimension and number of CKKS slots.

---

# 14. Step 7 — `RoundqQAlter()`

The final conversion is conditional:

```cpp
if (m_modulus_LWE != m_modulus_CKKS_from)
```

Then the code applies:

```cpp
RoundqQAlter(original_a[j], m_modulus_LWE, m_modulus_CKKS_from);
RoundqQAlter(original_b,   m_modulus_LWE, m_modulus_CKKS_from);
```

**Purpose:** Convert the extracted LWE values from the intermediate modulus `Q'` to the final LWE modulus `q` by scaling and rounding.

The intended operation is described in the source as multiplying by the ratio of the target and source moduli and rounding.

Conceptually:

```text
value modulo Q'
      │
      │ × q / Q'
      │ + rounding
      ▼
value modulo q
```

The result is then stored in a new `LWECiphertextImpl`.

---

# 15. CKKS → FHEW: Complete Detailed Tree

```text
SWITCHCKKSRNS::EvalCKKStoFHEW()
│
├── EvalSlotsToCoeffsSwitch()
│   │
│   └── EvalLTWithPrecomputeSwitch()
│       │
│       ├── EvalFastRotationPrecompute()
│       │
│       ├── EvalFastRotationExt()
│       │
│       ├── KeySwitchExt()
│       │
│       ├── FHECKKSRNS::EvalMultExt()
│       │   └── multiply rotated ciphertext by plaintext diagonal
│       │
│       ├── FHECKKSRNS::EvalAddExtInPlace()
│       │   └── accumulate diagonal-transform terms
│       │
│       ├── KeySwitchDownFirstElement()
│       │
│       └── KeySwitchDown()
│
├── ModReduceInternalInPlace()
│   └── reduce CKKS modulus level
│
├── ModSwitch()
│   └── prepare modulus Q'
│
├── ccKS->KeySwitch()
│   └── switch CKKS key → RLWE representation of LWE key
│
├── ExtractLWEpacked()
│   ├── RLWE A → coefficient vector
│   └── RLWE B → coefficient vector
│
├── ExtractLWECiphertext()
│   ├── choose coefficient slice for a
│   ├── use selected B coefficient
│   └── construct LWECiphertextImpl
│
└── RoundqQAlter()
    └── Q' → q when required
```

---

# 16. Precomputation Tree

The runtime transform depends on precomputed plaintext diagonals.

```text
EvalCKKStoFHEWPrecompute()
│
├── construct transform matrices / U0 and U1
│
└── EvalLTPrecomputeSwitch()
    │
    ├── determine CKKS modulus / RNS parameters
    ├── determine BSGS step
    ├── ExtractShiftedDiagonal()
    ├── scale diagonal coefficients
    ├── Fill()
    ├── Rotate()
    └── FHECKKSRNS::MakeAuxPlaintext()
           │
           └── produces plaintext representation used later
               by EvalLTWithPrecomputeSwitch()
```

The key relationship is:

```text
Precomputation
     │
     ▼
m_U0Pre / related diagonal plaintexts
     │
     ▼
EvalSlotsToCoeffsSwitch()
     │
     ▼
EvalLTWithPrecomputeSwitch()
```

---

# 17. `EvalLTPrecomputeSwitch()`

There are overloads for different matrix forms.

## 17.1 Matrix A + B form

```text
EvalLTPrecomputeSwitch(cc, A, B, dim1, L, scale)
```

Main work:

1. Obtain CKKS cryptographic and element parameters.
2. Determine how many RNS towers need to be dropped.
3. Collect the relevant moduli and roots.
4. Determine the baby-step size.
5. Concatenate matrix portions as needed.
6. Extract shifted diagonals.
7. Rotate/arrange the diagonal coefficients.
8. Create auxiliary plaintexts using `MakeAuxPlaintext()`.

The important output is a vector of plaintexts representing the matrix diagonals in the form expected by the linear-transform evaluator.

## 17.2 Square A form

```text
EvalLTPrecomputeSwitch(cc, A, dim1, L, scale)
```

This performs the analogous work for a square matrix.

---

# 18. `ExtractShiftedDiagonal()`

**Purpose:** Convert a matrix representation into the diagonal associated with a particular cyclic shift.

This is necessary because the linear-transform evaluator is using a diagonal method rather than directly multiplying a ciphertext by a full matrix.

Conceptually:

```text
Matrix M

m00 m01 m02 ...
m10 m11 m12 ...
m20 m21 m22 ...
...

      ↓ shifted diagonals

D0, D1, D2, ...
```

Then the transform can be represented as:


M x = Σᵢ D_i Rotᵢ(x).


---

# 19. `FHECKKSRNS::MakeAuxPlaintext()`

**Purpose:** Build the plaintext object used by the RNS CKKS linear-transform machinery from the prepared diagonal coefficient vector.

The precomputation code supplies:

- CKKS crypto context
- element parameters
- transformed coefficient vector
- scaling information
- towers to drop
- transform modulus information

The result is stored in precomputed structures such as `m_U0Pre` and later consumed by `EvalLTWithPrecomputeSwitch()`.

---

# 20. Reverse Direction: FHEW → CKKS

The reverse switch has a different structure.

```text
SWITCHCKKSRNS::EvalFHEWtoCKKS()
│
├── collect LWE a vectors into matrix A
├── collect LWE b values into vector B
│
├── compute prescale
│
└── EvalPartialHomDecryption()
       │
       ├── prepare / pad matrix
       ├── EvalLTRectPrecomputeSwitch()
       │      ├── prepare rectangular transform
       │      ├── ExtractShiftedDiagonal()
       │      └── MakeAuxPlaintext()
       │
       └── EvalLTRectWithPrecomputeSwitch()
              ├── EvalFastRotationPrecompute()
              ├── EvalFastRotationExt()
              ├── KeySwitchExt()
              ├── EvalMultExt()
              ├── EvalAddExtInPlace()
              └── KeySwitchDown()

                ↓

        CKKS ciphertext
```

The important conceptual difference is that FHEW → CKKS is **not** ordinary decryption followed by re-encryption.

Instead, the LWE decryption expression


b - \langle a,s \rangle


is evaluated homomorphically in CKKS using an encoded representation of the LWE secret key.

---

# 21. `EvalFHEWtoCKKS()`

**Role:** Convert a vector of LWE ciphertexts back into a CKKS ciphertext.

The function builds:

```text
A = matrix of LWE a-vectors
B = vector of LWE b-values
```

Then it applies a prescaling factor involving the LWE modulus and the switching scale and invokes:

```cpp
EvalPartialHomDecryption(
    *ccCKKS,
    A,
    m_FHEWtoCKKSswk,
    dim1,
    prescale,
    0
)
```

The result represents the homomorphic evaluation of the LWE decryption expression.

Conceptually:

```text
LWE ciphertexts
   │
   ├── A matrix
   └── B vector
        │
        ▼
CKKS homomorphic computation of B - A·s
        │
        ▼
CKKS ciphertext
```

---

# 22. `EvalPartialHomDecryption()`

**Purpose:** Evaluate the inner-product part of LWE decryption homomorphically.

The code prepares the matrix for a rectangular linear transform, including padding to a suitable power-of-two size when necessary.

Then:

```text
EvalPartialHomDecryption()
│
├── pad matrix if necessary
│
├── EvalLTRectPrecomputeSwitch()
│
└── EvalLTRectWithPrecomputeSwitch()
```

The mathematical target is the vector of inner products:


\langle a_i, s\rangle.


This produces the quantity needed to combine with `b_i` and recover the encoded plaintext value under the CKKS computation.

---

# 23. `EvalLTRectPrecomputeSwitch()`

This is the rectangular-matrix equivalent of the ordinary linear-transform precomputation.

It prepares the matrix diagonals for a rectangular transform and converts them into plaintext objects suitable for the CKKS evaluator.

The important operations are:

```text
matrix
  ↓
shifted diagonals
  ↓
scaling / padding
  ↓
auxiliary plaintexts
  ↓
precomputed rectangular transform
```

---

# 24. `EvalLTRectWithPrecomputeSwitch()`

This evaluates the rectangular transform using the same general baby-step/giant-step strategy.

Its call pattern is approximately:

```text
EvalLTRectWithPrecomputeSwitch()
│
├── EvalFastRotationPrecompute()
├── EvalFastRotationExt()
├── KeySwitchExt()
├── EvalMultExt()
├── EvalAddExtInPlace()
└── KeySwitchDown()
```

The difference from the square transform is the shape of the matrix and the number of output rows/columns that must be handled.

---

# 25. Application Example: `EvalMinSchemeSwitching()`

One of the important callers is the minimum-selection flow.

```text
EvalMinSchemeSwitching()
│
├── EvalSub()
│     └── compute pairwise CKKS differences
│
├── EvalAtIndex()
│     └── obtain rotated comparison values
│
├── EvalCKKStoFHEW()
│     └── convert comparison values to FHEW/LWE
│
├── m_ccLWE->EvalSign()
│     └── perform sign evaluation in FHEW
│
├── EvalFHEWtoCKKS()
│     └── return comparison / selection result to CKKS
│
└── use result to select the minimum-related value
```

The high-level computation is therefore:

```text
CKKS values
   │
   ▼
construct differences
   │
   ▼
CKKS → FHEW
   │
   ▼
FHEW EvalSign()
   │
   ▼
FHEW → CKKS
   │
   ▼
selection / minimum computation
```

---

# 26. Where the Functions Live

| Function | Main source / role |
|---|---|
| `SWITCHCKKSRNS::EvalCKKStoFHEW` | `ckksrns-schemeswitching.cpp` — top-level CKKS → FHEW switch |
| `SWITCHCKKSRNS::EvalSlotsToCoeffsSwitch` | `ckksrns-schemeswitching.cpp` — slots → coefficients |
| `SWITCHCKKSRNS::EvalLTPrecomputeSwitch` | `ckksrns-schemeswitching.cpp` — linear-transform precomputation |
| `SWITCHCKKSRNS::EvalLTWithPrecomputeSwitch` | `ckksrns-schemeswitching.cpp` — evaluates precomputed transform |
| `SWITCHCKKSRNS::ExtractLWEpacked` | `ckksrns-schemeswitching.cpp` — extracts A/B coefficient vectors |
| `SWITCHCKKSRNS::ExtractLWECiphertext` | `ckksrns-schemeswitching.cpp` — constructs an LWE ciphertext |
| `switchingKeyGenRLWEcc` | `ckksrns-schemeswitching.cpp` — generates CKKS → RLWE(LWE-key) switching key |
| `SWITCHCKKSRNS::EvalFHEWtoCKKS` | `ckksrns-schemeswitching.cpp` — top-level FHEW → CKKS switch |
| `EvalPartialHomDecryption` | `ckksrns-schemeswitching.cpp` — homomorphic LWE decryption component |
| `EvalLTRectPrecomputeSwitch` | `ckksrns-schemeswitching.cpp` — rectangular transform precompute |
| `EvalLTRectWithPrecomputeSwitch` | `ckksrns-schemeswitching.cpp` — rectangular transform evaluation |
| `FHECKKSRNS::MakeAuxPlaintext` | `ckksrns-fhe.cpp/.h` — RNS CKKS plaintext preparation |
| `FHECKKSRNS::EvalMultExt` | `ckksrns-fhe.cpp/.h` — extended CKKS multiplication |
| `FHECKKSRNS::EvalAddExtInPlace` | `ckksrns-fhe.cpp/.h` — extended CKKS addition |
| `FHECKKSRNS::KeySwitchDown` | CKKS/RNS implementation — returns from extended representation |
| `BinFHEContext::EvalSign` | BinFHE/FHEW implementation — sign evaluation |
| `FHERNS` | `rns-fhe.h` — RNS FHE base/type layer |
| `SchSwchParams` | `scheme-swch-params.h` — scheme-switching configuration |

---

# 27. One-Page Mental Model

When reading `ckksrns-schemeswitching.cpp`, keep this picture in mind:

```text
                    CKKS WORLD
                       │
                       │
             packed slot values
                       │
                       ▼
          EvalSlotsToCoeffsSwitch
                       │
                       ▼
            diagonal linear transform
                       │
                       ▼
          coefficient-oriented CKKS
                       │
              modulus preparation
                       │
                       ▼
                  KeySwitch
                       │
                       ▼
            RLWE under LWE-key form
                       │
               coefficient extraction
                       │
                       ▼
                 packed A/B data
                       │
                 LWE extraction
                       │
                       ▼
                LWE ciphertexts
                       │
                       ▼
                    FHEW WORLD
                       │
                  EvalSign(), etc.
                       │
                       ▼
                LWE ciphertexts
                       │
                       ▼
             homomorphic decryption
                       │
                       ▼
            rectangular CKKS transform
                       │
                       ▼
                 CKKS ciphertext
                       │
                       ▼
                    CKKS WORLD
```

---

# 28. The Most Important Distinctions

## `EvalSlotsToCoeffsSwitch()` vs `EvalLTWithPrecomputeSwitch()`

`EvalSlotsToCoeffsSwitch()` is the **specific scheme-switching transform wrapper**.

`EvalLTWithPrecomputeSwitch()` is the **generic engine that evaluates a precomputed linear transform**.

---

## `EvalLTPrecomputeSwitch()` vs `EvalLTWithPrecomputeSwitch()`

```text
EvalLTPrecomputeSwitch()
        ↓
prepare plaintext diagonals
        ↓
EvalLTWithPrecomputeSwitch()
        ↓
actually evaluate the transform
```

So one is the **setup/precomputation stage**, while the other is the **runtime evaluation stage**.

---

## `KeySwitch()` vs `KeySwitchExt()` / `KeySwitchDown()`

`KeySwitch()` in the top-level switching flow changes the ciphertext's key domain.

The `Ext` and `KeySwitchDown` operations are internal pieces of the RNS CKKS extended evaluation path used during linear-transform computation.

They should not be mentally treated as the same operation even though all involve key switching.

---

## `ExtractLWEpacked()` vs `ExtractLWECiphertext()`

```text
ExtractLWEpacked()
    = unpack RLWE polynomials into A/B coefficient arrays

ExtractLWECiphertext()
    = take one slice/index from those arrays and build one LWE ciphertext
```

---

# 29. Minimal Execution Sequence to Memorize

For CKKS → FHEW:

```text
EvalCKKStoFHEW
    ↓
EvalSlotsToCoeffsSwitch
    ↓
EvalLTWithPrecomputeSwitch
    ↓
ModReduce / ModSwitch
    ↓
KeySwitch
    ↓
ExtractLWEpacked
    ↓
ExtractLWECiphertext
    ↓
RoundqQAlter
    ↓
FHEW/LWE ciphertexts
```

For FHEW → CKKS:

```text
EvalFHEWtoCKKS
    ↓
EvalPartialHomDecryption
    ↓
EvalLTRectWithPrecomputeSwitch
    ↓
CKKS ciphertext
```

For a min/sign application:

```text
EvalMinSchemeSwitching
    ↓
EvalCKKStoFHEW
    ↓
EvalSign
    ↓
EvalFHEWtoCKKS
```

---

# 30. Source Notes

The source files relevant to this explanation are:

- `ckksrns-schemeswitching.cpp`
- `ckksrns-fhe.cpp`
- `ckksrns-fhe.h`
- `scheme-swch-params.h`
- `rns-fhe.h`

The key implementation facts used here include the direct calls visible in the supplied scheme-switching implementation: `EvalSlotsToCoeffsSwitch()`, `ModReduceInternalInPlace()`, `ModSwitch()`, `KeySwitch()`, `ExtractLWEpacked()`, `ExtractLWECiphertext()`, and the conditional `RoundqQAlter()` path, together with the precomputed linear-transform machinery and the reverse `EvalFHEWtoCKKS()` path.
