fastQR: An R Package for fast QR decomposition and update
================
Mauro Bernardi
November 13, 2024

[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/0.1.0/active.svg)](http://www.repostatus.org/#active)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/fdaSP)](https://CRAN.R-project.org/package=fdaSP)
![downloads](http://cranlogs.r-pkg.org/badges/grand-total/fdaSP)

# QR decomposition

The QR decomposition algorithm decomposes a matrix
$\mathbf{X}\in\mathbb{R}^{n\times p}$ into the product of two matrices,
namely, an orthonormal matrix $\mathbf{Q}\in\mathbb{R}^{n\times n}$ and
a trapezoidal matrix $\mathbf{R}\in\mathbb{R}^{n\times p}$,
i.e. $\mathbf{X}=\mathbf{Q}\mathbf{R}$.

**The Cauchy-Schwarz Inequality**
$$\left( \sum_{k=1}^n a_k b_k \right)^2 \leq \left( \sum_{k=1}^n a_k^2 \right) \left( \sum_{k=1}^n b_k^2 \right)$$

## Motivation

The QR decomposition is the foundation algorithm of least squares
regression, and is also used in the QR eigenvalue algorithm.

## Usage
