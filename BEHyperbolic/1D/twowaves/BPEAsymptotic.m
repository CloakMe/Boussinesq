function cgreturn = BPEAsymptotic(q,k)

  cgreturn = -0.2e1 * k ^ 2 / (exp(-0.1e1 * q) + 0.2e1 + exp(q));
end