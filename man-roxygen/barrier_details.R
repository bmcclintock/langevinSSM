#' @details
#' \strong{Barrier Constraints:}
#' When a \code{barrier} is specified, the model restricts animal movement using a physical penalty method for stochastic differential equations (SDEs).
#' Internally, the binary mask provided in \code{spatialCovs[[barrier]]} (where 1 indicates unrestricted areas and 0 indicates restricted areas) is converted into a continuous signed distance function (SDF), denoted as \eqn{d(\boldsymbol{\mu})}. Distance is calculated as positive in the unrestricted areas (allowed space) and negative in the restricted areas.
#'
#' To prevent the animal from entering the restricted zone (where \eqn{d(\boldsymbol{\mu}) \le 0}), the model applies a harmonic oscillator penalty, treating the boundary like a spring. The potential energy (\eqn{U}) of this spring is a quadratic function of the distance traveled into the restricted zone:
#' \deqn{U(\boldsymbol{\mu}) = \frac{1}{2} \lambda d(\boldsymbol{\mu})^2}
#'
#' In Langevin dynamics, the force acting on the animal is the negative spatial gradient of the potential energy (\eqn{ -\nabla U(\boldsymbol{\mu})}). Applying the chain rule to the potential energy function yields the exact force applied during simulation and model fitting:
#' \deqn{\text{Force} = -\lambda d(\boldsymbol{\mu}) \nabla d(\boldsymbol{\mu})}
#'
#' Because \eqn{d(\boldsymbol{\mu})} is negative in the restricted zone, the equation produces a positive force that pushes the animal in the direction of the shortest path back to the unrestricted area (\eqn{\nabla d(\boldsymbol{\mu})}).
#'
#' The \code{lambda} parameter acts as the spring constant (\eqn{\lambda}); higher values create a stiffer, harder boundary. This boundary force is added directly to the spatial gradient of the habitat covariates during optimiziation, with a corresponding selection coefficient (e.g., \code{beta_barrier} if \code{barrier = "barrier"}). Because \code{lambda} is a mathematical property of the SDE rather than a biological parameter, it cannot be safely estimated via maximum likelihood estimation. It must be provided as a fixed hyperparameter (see \code{\link{checkBarrier}} for a strategy to determine the optimal penalty).
