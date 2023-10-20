classdef (ConstructOnLoad) IC_2Waves
    
methods( Access = public )
  function [ u, dudt ] = GetInitialCondition(this, x, t, k, b, a1, a2, a12)
    if(nargin == 3)
          k = sqrt(1/8);
          b = 0;
          a1 = 1; 
          a2 = 1; 
          a12 = 1;
    end
    u = 0.2e1 * (a1 * (k ^ 2) * exp((k * x) + sqrt((k ^ 2 * (1 - k ^ 2))) * t + b) + a2 * (k ^ 2) * exp(-(k * x) +...
        sqrt((k ^ 2 * (1 - k ^ 2))) * t + b)) ./...
        ((0.4e1 * a1 * a2 * (k ^ 2) - a1 * a2) / (-a12 + a12 * k ^ 2) + a1 * exp((k * x) + sqrt((k ^ 2 * (1 - k ^ 2))) * t + b) +...
        a2 * exp(-(k * x) + sqrt((k ^ 2 * (1 - k ^ 2))) * t + b) + a12 * exp(0.2e1 * sqrt((k ^ 2 * (1 - k ^ 2))) * t + 0.2e1 * b)) -...
        0.2e1 * (a1 * k * exp((k * x) + sqrt((k ^ 2 * (1 - k ^ 2))) * t + b) - a2 * k * exp(-(k * x) + sqrt((k ^ 2 * (1 - k ^ 2))) * t + b)) .^ 2 ./...
        ((0.4e1 * a1 * a2 * (k ^ 2) - a1 * a2) / (-a12 + a12 * k ^ 2) + a1 * exp((k * x) + sqrt((k ^ 2 * (1 - k ^ 2))) * t + b) + a2 * exp(-(k * x) + sqrt((k ^ 2 * (1 - k ^ 2))) * t + b) + a12 * exp(0.2e1 * sqrt((k ^ 2 * (1 - k ^ 2))) * t + 0.2e1 * b)) .^ 2;
    dudt = -0.2e1 * (k ^ 2) * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * a12 * (-1 + k ^ 2) * (-0.16e2 * a1 ^ 3 * exp((k * x) + sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + b) * a2 ^ 2 * (k ^ 4) + 0.8e1 * a1 ^ 3 * exp((k * x) + sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + b) * a2 ^ 2 * (k ^ 2) - 0.16e2 * a2 ^ 3 * exp(-(k * x) + sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + b) * a1 ^ 2 * (k ^ 4) + 0.8e1 * a2 ^ 3 * exp(-(k * x) + sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + b) * a1 ^ 2 * (k ^ 2) + 0.2e1 * a1 ^ 2 * exp((2 * k * x) + 0.4e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.4e1 * b) * a12 ^ 3 * (k ^ 2) + 0.2e1 * a2 ^ 2 * exp(-(2 * k * x) + 0.4e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.4e1 * b) * a12 ^ 3 * (k ^ 2) - 0.6e1 * a1 ^ 2 * exp(0.2e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.2e1 * b) * a2 ^ 2 * a12 - 0.2e1 * a1 * exp((k * x) + 0.5e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.5e1 * b) * a12 ^ 4 * (k ^ 2) - 0.2e1 * a2 * exp(-(k * x) + 0.5e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.5e1 * b) * a12 ^ 4 * (k ^ 2) + 0.6e1 * a1 * exp(0.4e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.4e1 * b) * a2 * a12 ^ 3 + a2 ^ 3 * exp(-(2 * k * x) + 0.2e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.2e1 * b) * a1 * a12 + a1 ^ 3 * exp((2 * k * x) + 0.2e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.2e1 * b) * a2 * a12 - 0.12e2 * a1 * exp(0.4e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.4e1 * b) * a2 * a12 ^ 3 * (k ^ 2) + 0.6e1 * a1 * exp(0.4e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.4e1 * b) * a2 * a12 ^ 3 * (k ^ 4) - 0.5e1 * a1 ^ 3 * exp((2 * k * x) + 0.2e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.2e1 * b) * a2 * (k ^ 2) * a12 + 0.4e1 * a1 ^ 3 * exp((2 * k * x) + 0.2e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.2e1 * b) * a2 * (k ^ 4) * a12 + 0.4e1 * a2 ^ 3 * exp(-(2 * k * x) + 0.2e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.2e1 * b) * a1 * (k ^ 4) * a12 - 0.24e2 * a1 ^ 2 * exp(0.2e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.2e1 * b) * a2 ^ 2 * (k ^ 4) * a12 - 0.5e1 * a2 ^ 3 * exp(-(2 * k * x) + 0.2e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.2e1 * b) * a1 * (k ^ 2) * a12 - a2 ^ 2 * exp(-(2 * k * x) + 0.4e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.4e1 * b) * a12 ^ 3 + a2 * exp(-(k * x) + 0.5e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.5e1 * b) * a12 ^ 4 + 0.30e2 * a1 ^ 2 * exp(0.2e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.2e1 * b) * a2 ^ 2 * (k ^ 2) * a12 - a1 ^ 2 * exp((2 * k * x) + 0.4e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.4e1 * b) * a12 ^ 3 + a1 * exp((k * x) + 0.5e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.5e1 * b) * a12 ^ 4 - a1 ^ 3 * exp((k * x) + sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + b) * a2 ^ 2 - a2 ^ 3 * exp(-(k * x) + sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + b) * a1 ^ 2 - a1 ^ 2 * exp((2 * k * x) + 0.4e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.4e1 * b) * a12 ^ 3 * (k ^ 4) + a2 * exp(-(k * x) + 0.5e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.5e1 * b) * a12 ^ 4 * (k ^ 4) - a2 ^ 2 * exp(-(2 * k * x) + 0.4e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.4e1 * b) * a12 ^ 3 * (k ^ 4) + a1 * exp((k * x) + 0.5e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.5e1 * b) * a12 ^ 4 * (k ^ 4)) ./ (0.4e1 * a1 * a2 * (k ^ 2) - a1 * a2 - a1 * exp((k * x) + sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + b) * a12 + a1 * exp((k * x) + sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + b) * a12 * (k ^ 2) - a2 * exp(-(k * x) + sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + b) * a12 + a2 * exp(-(k * x) + sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + b) * a12 * (k ^ 2) - a12 ^ 2 * exp(0.2e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.2e1 * b) + a12 ^ 2 * exp(0.2e1 * sqrt(-(k ^ 2 * (-1 + k ^ 2))) * t + 0.2e1 * b) * (k ^ 2)) .^ 3;

  end
  
  function [u, dudt ] = GetInitialCondition2w(this, x, t, k1, k2, b1, b2, a1, a2, a12)
    if(nargin == 3)
        k1 = 1/3;
        k2 = -1/2;
        b1 = 0;
        b2 = 0;
        a1 = 1; 
        a2 = 1; 
        a12 = 1;
    end
     u = -(-0.2e1 .* (a1 * (k1 ^ 2) .* exp((k1 .* x) + sqrt((k1 .^ 2 .* (1 - k1 ^ 2))) .* t + b1) + a2 .* (k2 .^ 2) .* exp((k2 .* x) +...
           sqrt((k2 .^ 2 .* (1 - k2 .^ 2))) .* t + b2) + a12 .* ((k1 + k2) .^ 2) .* exp(((k1 + k2) .* x) + (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) +...
           sqrt((k2 .^ 2 .* (1 - k2 .^ 2)))) .* t + b1 + b2)) ./ (a1 .* a2 .* ((4 * k1 ^ 4) + (4 .* k2 .^ 4) - (3 .* k1 ^ 2) - (3 .* k2 ^ 2) -...
           (2 .* k1 ^ 2 .* k2 ^ 2) + 0.6e1 .* sqrt((-k1 ^ 4 + k1 ^ 2)) .* sqrt((-k2 .^ 4 + k2 ^ 2))) ./ a12 / ((k1 + k2) ^ 2) ./ ...
           (4 * k1 ^ 2 + 4 * k1 * k2 + 4 * k2 ^ 2 - 3) + a1 .* exp((k1 .* x) + sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* t + b1) + a2 .* exp((k2 .* x) +...
           sqrt((k2 ^ 2 .* (1 - k2 ^ 2))) .* t + b2) + a12 .* exp(((k1 + k2) .* x) + (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 .* (1 - k2 ^ 2)))) .* t +...
           b1 + b2)) + 0.2e1 .* (a1 .* k1 .* exp((k1 .* x) + sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* t + b1) + a2 .* k2 .* exp((k2 .* x) + ...
           sqrt((k2 ^ 2 .* (1 - k2 ^ 2))) .* t + b2) + a12 .* (k1 + k2) .* exp(((k1 + k2) .* x) + (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + ...
           sqrt((k2 ^ 2 .* (1 - k2 ^ 2)))) .* t + b1 + b2)) .^ 2 ./ (a1 .* a2 .* ((4 * k1 ^ 4) + (4 .* k2 .^ 4) - (3 .* k1 ^ 2) - (3 .* k2 ^ 2) -...
           (2 .* k1 ^ 2 .* k2 ^ 2) + 0.6e1 .* sqrt((-k1 ^ 4 + k1 ^ 2)) .* sqrt((-k2 .^ 4 + k2 ^ 2))) ./ a12 ./ ((k1 + k2) ^ 2) ./ (4 * k1 ^ 2 + 4 * k1 * k2 + 4 * k2 ^ 2 - 3) +...
           a1 .* exp((k1 .* x) + sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* t + b1) + a2 .* exp((k2 .* x) + sqrt((k2 ^ 2 * (1 - k2 ^ 2))) .* t + b2) + a12 .* exp(((k1 + k2) .* x) + ...
           (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 * (1 - k2 ^ 2)))) .* t + b1 + b2)) .^ 2);
        
     dudt = -(-0.2e1 .* (a1 * (k1 ^ 2) .* sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* exp((k1 .* x) + sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* t + b1) + a2 .* (k2 ^ 2) .* sqrt((k2 ^ 2 * (1 - k2 ^ 2))) .* exp((k2 .* x) + sqrt((k2 ^ 2 * (1 - k2 ^ 2))) .* t + b2) + a12 .* ((k1 + k2) ^ 2) .* (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 * (1 - k2 ^ 2)))) .* exp(((k1 + k2) .* x) + (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 * (1 - k2 ^ 2)))) .* t + b1 + b2)) ./ (a1 .* a2 .* ((4 * k1 ^ 4) + (4 .* k2 .^ 4) - (3 .* k1 ^ 2) - (3 .* k2 ^ 2) - (2 .* k1 ^ 2 .* k2 ^ 2) + 0.6e1 .* sqrt((-k1 ^ 4 + k1 ^ 2)) .* sqrt((-k2 .^ 4 + k2 ^ 2))) ./ a12 ./ ((k1 + k2) ^ 2) ./ (4 * k1 ^ 2 + 4 * k1 * k2 + 4 * k2 ^ 2 - 3) + a1 .* exp((k1 .* x) + sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* t + b1) + a2 .* exp((k2 .* x) + sqrt((k2 ^ 2 * (1 - k2 ^ 2))) .* t + b2) + a12 .* exp(((k1 + k2) .* x) + (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 * (1 - k2 ^ 2)))) .* t + b1 + b2)) + 0.2e1 .* (a1 * (k1 ^ 2) .* exp((k1 .* x) + sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* t + b1) + a2 .* (k2 ^ 2) .* exp((k2 .* x) + sqrt((k2 ^ 2 * (1 - k2 ^ 2))) .* t + b2) + a12 .* ((k1 + k2) ^ 2) .* exp(((k1 + k2) .* x) + (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 * (1 - k2 ^ 2)))) .* t + b1 + b2)) ./ (a1 .* a2 .* ((4 * k1 ^ 4) + (4 .* k2 .^ 4) - (3 .* k1 ^ 2) - (3 .* k2 ^ 2) - (2 .* k1 ^ 2 .* k2 ^ 2) + 0.6e1 .* sqrt((-k1 ^ 4 + k1 ^ 2)) .* sqrt((-k2 .^ 4 + k2 ^ 2))) ./ a12 ./ ((k1 + k2) ^ 2) ./ (4 * k1 ^ 2 + 4 * k1 * k2 + 4 * k2 ^ 2 - 3) + a1 .* exp((k1 .* x) + sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* t + b1) + a2 .* exp((k2 .* x) + sqrt((k2 ^ 2 * (1 - k2 ^ 2))) .* t + b2) + a12 .* exp(((k1 + k2) .* x) + (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 * (1 - k2 ^ 2)))) .* t + b1 + b2)) .^ 2 .* (a1 .* sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* exp((k1 .* x) + sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* t + b1) + a2 .* sqrt((k2 ^ 2 * (1 - k2 ^ 2))) .* exp((k2 .* x) + sqrt((k2 ^ 2 * (1 - k2 ^ 2))) .* t + b2) + a12 .* (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 * (1 - k2 ^ 2)))) .* exp(((k1 + k2) .* x) + (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 * (1 - k2 ^ 2)))) .* t + b1 + b2)) + 0.4e1 .* (a1 .* k1 .* exp((k1 .* x) + sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* t + b1) + a2 .* k2 .* exp((k2 .* x) + sqrt((k2 ^ 2 * (1 - k2 ^ 2))) .* t + b2) + a12 .* (k1 + k2) .* exp(((k1 + k2) .* x) + (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 * (1 - k2 ^ 2)))) .* t + b1 + b2)) ./ (a1 .* a2 .* ((4 * k1 ^ 4) + (4 .* k2 .^ 4) - (3 .* k1 ^ 2) - (3 .* k2 ^ 2) - (2 .* k1 ^ 2 .* k2 ^ 2) + 0.6e1 .* sqrt((-k1 ^ 4 + k1 ^ 2)) .* sqrt((-k2 .^ 4 + k2 ^ 2))) ./ a12 ./ ((k1 + k2) ^ 2) ./ (4 * k1 ^ 2 + 4 * k1 * k2 + 4 * k2 ^ 2 - 3) + a1 .* exp((k1 .* x) + sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* t + b1) + a2 .* exp((k2 .* x) + sqrt((k2 ^ 2 * (1 - k2 ^ 2))) .* t + b2) + a12 .* exp(((k1 + k2) .* x) + (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 * (1 - k2 ^ 2)))) .* t + b1 + b2)) .^ 2 .* (a1 .* k1 .* sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* exp((k1 .* x) + sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* t + b1) + a2 .* k2 .* sqrt((k2 ^ 2 * (1 - k2 ^ 2))) .* exp((k2 .* x) + sqrt((k2 ^ 2 * (1 - k2 ^ 2))) .* t + b2) + a12 .* (k1 + k2) .* (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 * (1 - k2 ^ 2)))) .* exp(((k1 + k2) .* x) + (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 * (1 - k2 ^ 2)))) .* t + b1 + b2)) - 0.4e1 .* (a1 .* k1 .* exp((k1 .* x) + sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* t + b1) + a2 .* k2 .* exp((k2 .* x) + sqrt((k2 ^ 2 * (1 - k2 ^ 2))) .* t + b2) + a12 .* (k1 + k2) .* exp(((k1 + k2) .* x) + (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 * (1 - k2 ^ 2)))) .* t + b1 + b2)) .^ 2 ./ (a1 .* a2 .* ((4 * k1 ^ 4) + (4 * k2 ^ 4) - (3 .* k1 ^ 2) - (3 .* k2 ^ 2) - (2 .* k1 ^ 2 .* k2 ^ 2) + 0.6e1 .* sqrt((-k1 ^ 4 + k1 ^ 2)) .* sqrt((-k2 .^ 4 + k2 ^ 2))) ./ a12 ./ ((k1 + k2) ^ 2) ./ (4 * k1 ^ 2 + 4 * k1 * k2 + 4 * k2 ^ 2 - 3) + a1 .* exp((k1 .* x) + sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* t + b1) + a2 .* exp((k2 .* x) + sqrt((k2 ^ 2 * (1 - k2 ^ 2))) .* t + b2) + a12 .* exp(((k1 + k2) .* x) + (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 * (1 - k2 ^ 2)))) .* t + b1 + b2)) .^ 3 .* (a1 .* sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* exp((k1 .* x) + sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) .* t + b1) + a2 .* sqrt((k2 ^ 2 * (1 - k2 ^ 2))) .* exp((k2 .* x) + sqrt((k2 ^ 2 * (1 - k2 ^ 2))) .* t + b2) + a12 .* (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 * (1 - k2 ^ 2)))) .* exp(((k1 + k2) .* x) + (sqrt((k1 ^ 2 .* (1 - k1 ^ 2))) + sqrt((k2 ^ 2 * (1 - k2 ^ 2)))) .* t + b1 + b2)));

  end
end
    
end