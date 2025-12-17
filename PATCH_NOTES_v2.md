# Patch notes for short_crap.pdf -> Protocol v2.0

These are the *minimum* edits to bring the PDF's protocol into a reproducible v2 state.

## 1) Recipe 1.1 (Finite Part fit): define sigma0 explicitly
- Current text: "fit AB(t) log(σ − σ0) + CB,N(t)"
- Patch: add a sentence defining σ0 as the actual pole location of the analytically continued test function.
- Also: forbid hardcoding σ0=0.5 unless the function genuinely has a pole at 0.5.

## 2) Recipe 3.2 (GUE Wigner surmise): fix the exponent typo
- Current: exp(−4 s^2 / π^2)
- Patch: exp(−4 s^2 / π)

Correct form:
  P_GUE(s) = (32/π^2) s^2 exp(−4 s^2 / π)

## 3) "Pole at 0.5" example: correct it
- Bad example: ζ(2s − 0.5) (pole at s=0.75)
- Good examples:
    ζ(2s)      (pole at s=0.5)
    ζ(s + 0.5) (pole at s=0.5)

## 4) Keep the cyclic-truncation warning as a hard constraint
Any banding or roots-of-unity spectrum arising from cyclic shifts is an artifact; require cross-block interference + controls.
