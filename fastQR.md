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
$\mathbf{A}\in\mathbb{R}^{n\times p}$ into the product of two matrices,
namely, an orthonormal matrix $\mathbf{Q}\in\mathbb{R}^{n\times n}$ and
a trapezoidal matrix $\mathbf{R}\in\mathbb{R}^{n\times p}$,
i.e. $\mathbf{A}=\mathbf{Q}\mathbf{R}$.

## Motivation

The QR decompositoin is the foundation algorithm of least squares
regression, and is also used in the QR eigenvalue algorithm.

## Usage
