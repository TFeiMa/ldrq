# ldrq
## Abstract
We present a linear-time approximate algorithm for the problem(RQ)of minimizing the ratio of quadratic functions over an ellipsoid. The algorithmfrstly utilizes the bisection search to fnd a solution of Lagrange dual problemof (RQ), and then the global mimimum of (RQ) can be attained. Moreover,we prove the strong duality between primary problem and dual problem as well as give the interval of optimal Lagrange multiplier.
## Introduction
we study the problem of minimizing the ratio of two quadratic functions subject to a single quadratic constraint:
<img src="http://www.forkosh.com/mathtex.cgi? \begin{equation}
(\text{RQ}): \quad
\begin{array}{c l}
\min &f(x) = \displaystyle \frac{ f_1(x)}{1+ \| x\|^2}, \\
\text{s.t.} & c(x) \leq 0,\\
\end{array}
\label{rq1}
\end{equation}">
