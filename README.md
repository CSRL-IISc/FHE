# Fully Homomorphic Encryption: FHE, CKKS, TFHE, and Scheme Switching in OpenFHE

## 1. What is Fully Homomorphic Encryption (FHE)?

Fully Homomorphic Encryption lets you compute on **encrypted data** without ever decrypting it. A server can take ciphertexts `Enc(a)` and `Enc(b)`, run some computation, and produce `Enc(f(a,b))` — all without learning `a`, `b`, or the result. Only the key holder can decrypt.

Almost all practical FHE schemes are built on variants of the **Learning With Errors (LWE)** problem or its ring/module analogues (**RLWE**, **MLWE**). Security comes from the hardness of recovering a secret vector from noisy linear equations — the "error"/"noise" is what hides the plaintext, but it also *grows* with every homomorphic operation. Managing that noise growth is the central engineering problem in FHE, and it's why different schemes exist: each makes different trade-offs about *what kind of data* you encrypt and *what kind of noise growth* you're willing to tolerate.

Three main scheme families dominate practical FHE today:

| Scheme | Plaintext type | Best at | Weak at |
|---|---|---|---|
| **BGV / BFV** | Integers (modular arithmetic) | Exact integer arithmetic, SIMD batching | Comparisons, non-polynomial functions |
| **CKKS** | Real / complex numbers (approximate) | Deep arithmetic circuits, ML inference, SIMD | Comparisons, exact results, bootstrapping cost |
| **TFHE / FHEW** | Single bits or small integers | Arbitrary boolean/lookup-table gates, cheap bootstrapping per gate | Very slow for large arithmetic circuits (bit-by-bit) |

CKKS↔TFHE switching and RLWE-to-LWE conversion sit exactly at the seam between the second and third rows — which is the focus of this document.

---

## 2. CKKS: Approximate Arithmetic on Encrypted Reals

CKKS (Cheon–Kim–Kim–Song) encrypts vectors of real/complex numbers and supports approximate addition and multiplication, similar to floating-point arithmetic but homomorphic.

**Core ideas:**
- Plaintext is encoded into a polynomial in `Z[X]/(X^N + 1)` via a **canonical embedding**, which packs up to `N/2` real (or complex) values into a single ciphertext — enabling SIMD-style batched arithmetic.
- A ciphertext is a pair `(c0, c1)` in `R_Q^2` (R = ring, Q = ciphertext modulus), and decryption recovers `m + e` for small noise `e` — the result is *approximate*, not exact, by design.
- **Rescaling**: after a multiplication, the ciphertext modulus is divided down (like floating-point exponent management) to control noise growth and keep the scale consistent. This consumes "levels" — each multiplication uses up one level of the modulus chain `Q = q_0 · q_1 · ... · q_L`.
- Once levels run out, **bootstrapping** refreshes the ciphertext back to a high modulus, at high computational cost — this is usually the single most expensive operation in a CKKS pipeline.

CKKS is the workhorse for encrypted machine learning (linear layers, convolutions, polynomial approximations of activation functions) because multiplications and additions are relatively cheap and SIMD packing amortizes cost across thousands of slots.

**What CKKS is bad at:** comparisons, max/min, ReLU, sign function — anything non-polynomial. You can approximate these with high-degree polynomials, but it's expensive and imprecise. This is exactly the gap TFHE fills.

---

## 3. TFHE / FHEW: Fast Bit-Level Bootstrapping

TFHE (and its close relative FHEW) take the opposite approach: encrypt single bits (or small integers) using plain **LWE** (not ring-LWE) ciphertexts, and make **bootstrapping itself cheap enough to run after every gate**.

**Core ideas:**
- An LWE ciphertext is `(a, b) = (a, ⟨a,s⟩ + m·Δ + e)` — a vector `a` of dimension `n`, plus one scalar `b`. Much smaller and simpler algebraic structure than an RLWE/CKKS ciphertext.
- **Programmable Bootstrapping (PBS)**: TFHE's signature contribution. Instead of bootstrapping being just noise-refresh, it's engineered so that the bootstrapping step *itself* evaluates an arbitrary **lookup table** (any function of the input bit/small-integer), while simultaneously refreshing the noise. This is done via:
  - **Blind rotation**: homomorphically rotating a test polynomial (encoding the LUT) by an encrypted amount, using a sequence of **external products** (GGSW × RLWE) driven by the bits of the LWE secret key, via the **bootstrapping key** (a set of RGSW/GGSW encryptions of the secret key bits).
  - **Sample extraction**: pulling a single LWE ciphertext (one coefficient) back out of the resulting RLWE ciphertext.
  - **Key switching**: converting the extracted LWE ciphertext (which is under a large-dimension RLWE-derived secret) back down to the original small-dimension LWE secret key, to keep sizes bounded for the next round.
- Because bootstrapping is fast (milliseconds) and does useful work (the LUT evaluation), TFHE is ideal for **boolean circuits, comparisons, and small lookup-table-based functions** — the exact operations where CKKS struggles.

**What TFHE is bad at:** large arithmetic circuits over big numbers. Doing a 32-bit multiplication bit-by-bit through gate bootstrapping is far slower than CKKS doing it as one SIMD-packed polynomial multiplication.

---

## 4. Why Scheme Switching? (CKKS ⇄ FHEW/TFHE)

Real workloads (especially ML inference) need **both** kinds of operations:
- Linear layers, convolutions, dot products → CKKS is fast (SIMD, polynomial arithmetic).
- Activation functions like ReLU/max/sign, comparisons, argmax → TFHE is fast (exact via LUT bootstrapping).

Rather than picking one scheme and suffering on the other's operations, **scheme switching** lets you convert ciphertexts back and forth mid-computation: run the linear parts in CKKS, hop over to FHEW/TFHE for the nonlinearity, then hop back to CKKS to continue. This is the "best of both worlds" approach, and it's the feature OpenFHE exposes.

This CKKS↔TFHE + RLWE-to-LWE conversion pipeline is a natural target for hardware acceleration — OpenFHE's software implementation effectively serves as the reference algorithm for mapping onto FPGA/hardware datapaths (NTT, base conversion, blind rotation scheduling, etc.).

### 4.1 A note on naming: "FHEW" vs "TFHE" in OpenFHE's switching API

OpenFHE's switching functions are named `EvalCKKStoFHEW` / `EvalFHEWtoCKKS` — but for work specifically on **CKKS ↔ TFHE**, these are still the correct functions to use. Here's why:

- **FHEW** and **TFHE** are distinct schemes with different bootstrapping algorithms (FHEW's original **AP** accumulator update vs. TFHE's **GINX/CGGI-style** blind rotation), but both operate on **the same underlying object**: a plain LWE ciphertext `(a, b)` under a binary/ternary secret key.
- OpenFHE unifies both under one module — `src/binfhe/` — and one context class, `BinFHEContext`, where you select the bootstrapping method as a parameter (`AP` for FHEW-style, `GINX` for TFHE-style) at setup time. Everything downstream — the ciphertext format, key structures, and crucially the **CKKS-switching code** — is agnostic to which method you picked.
- In other words, `EvalCKKStoFHEW` / `EvalFHEWtoCKKS` don't mean "only compatible with FHEW's bootstrapping." The name reflects the shared LWE/`binfhe` API layer (which was originally built around the FHEW scheme and kept its name), not a restriction to FHEW specifically. **Configuring `BinFHEContext` with `GINX` makes the pipeline CKKS↔TFHE switching, using these exact same functions.**
- So concretely: the extraction/packing math (`PackLWEs`, field trace, modulus switching) is identical regardless of which bootstrapping style sits on the FHEW/TFHE side — the only thing that changes is what happens *inside* the bootstrap once on that side (AP's vs. GINX's blind-rotation accumulator update). For hardware design purposes, this means the CKKS-side extraction/packing datapath is reusable across both, and the only scheme-specific hardware is the blind-rotation core itself.

---

## 5. How OpenFHE Implements Scheme Switching

OpenFHE (the merged successor of PALISADE + HElib-style tooling) implements CKKS↔FHEW/TFHE switching primarily in:

```
src/pke/include/scheme/ckksrns/ckksrns-schemeswitching.h
src/pke/lib/scheme/ckksrns/ckksrns-schemeswitching.cpp
```

with supporting FHEW-side logic in `src/binfhe/`. The example driver program is typically:

```
src/pke/examples/scheme-switching.cpp
```

### 5.1 CKKS → FHEW (the "extraction" direction)

Conceptually, this direction needs to pull **individual real-valued slots** out of a SIMD-packed CKKS ciphertext and re-encrypt each one as a small, single-value LWE ciphertext under FHEW's parameters. Key steps in OpenFHE's pipeline:

1. **`EvalCKKStoFHEWSetup(...)`** — generates the necessary switching keys: a key-switching key from the CKKS RLWE secret to an intermediate representation, and the FHEW bootstrapping/refreshing keys needed on the other side.
2. **Homomorphic slot extraction**: since CKKS packs many values per ciphertext, the library uses techniques closely related to **`PackLWEs`/`EvalTr` (partial/field trace)** operations to isolate a single slot's contribution.
3. **`EvalCKKStoFHEW(...)`** — performs the actual conversion, producing a vector of LWE ciphertexts (one FHEW ciphertext per extracted CKKS slot), rescaled to FHEW's plaintext modulus.
4. Internally this relies on a **modulus switch** down to FHEW's much smaller ciphertext modulus, since CKKS operates over a large RNS modulus chain and FHEW/TFHE ciphertexts must stay small for cheap bootstrapping.

### 5.2 FHEW → CKKS (the "packing" direction)

The reverse direction takes many small FHEW/LWE ciphertexts (e.g., after evaluating a nonlinear function bit-by-bit or value-by-value) and **repacks** them into a single SIMD CKKS ciphertext:

1. **`EvalFHEWtoCKKSSetup(...)`** — generates the packing keys.
2. **`EvalFHEWtoCKKS(...)`** — performs the RLWE packing: each LWE ciphertext is homomorphically embedded into a coefficient/slot of a fresh RLWE (CKKS) ciphertext, then combined via automorphisms (Galois key rotations) and additions — this is the **PackLWEs** algorithm, which is the LWE→RLWE analogue of the extraction step above, run in reverse.
3. Because FHEW ciphertexts individually carry noise from bootstrapping, OpenFHE also handles rescaling back into CKKS's modulus chain so the packed result behaves like a normal fresh(ish) CKKS ciphertext, ready for further homomorphic arithmetic.

### 5.3 The glue: a shared, compatible parameter world

For switching to be efficient at all, both schemes need **compatible ring dimensions and moduli** — OpenFHE's setup functions (`EvalSchemeSwitchingSetup`, `EvalCompareSwitchPrecompute`, etc.) negotiate a shared parameter set so that:
- The RLWE ring used for the CKKS side and the RLWE ring used internally by FHEW's blind rotation (Ring-GSW bootstrapping) can interoperate.
- Automorphism/Galois keys used for packing are consistent with the CKKS rotation keys already in use.
- The whole pipeline can target **128-bit security** end-to-end — but OpenFHE's default FHEW/CKKS combinations aren't automatically jointly 128-bit secure; parameters need to be chosen carefully for both schemes together.




## 7. Suggested Reading Order

1. CKKS original paper — Cheon, Kim, Kim, Song, *"Homomorphic Encryption for Arithmetic of Approximate Numbers"* (ASIACRYPT 2017)
2. TFHE original paper — Chillotti et al., *"TFHE: Fast Fully Homomorphic Encryption over the Torus"*
3. CKKS↔FHEW switching — Lu et al. / Bossuat, Cong, et al., *"Efficiently Switching between Fully Homomorphic Encryption Schemes"* (the paper OpenFHE's implementation is largely based on)
4. `PackLWEs` / trace-based packing — Micciancio–Polyakov style RLWE packing algorithms
5. OpenFHE documentation: https://openfhe-development.readthedocs.io

---

