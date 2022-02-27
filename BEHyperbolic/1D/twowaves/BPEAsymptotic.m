function cgreturn = BPEAsymptotic(p, q, k, a1, a2, a12)

  cgreturn = -0.2e1 * a1 * a12 * k ^ 2 / (a1 ^ 2 * exp(-q) + 0.2e1 * a1 * a12 + a12 ^ 2 * exp(q)) ...
             -0.2e1 * a2 * a12 * k ^ 2 / (a2 ^ 2 * exp(-p) + 0.2e1 * a12 * a2 + a12 ^ 2 * exp(p));
end