//TODO: I should create a static function that computes the coefficients.
//Then, having a basic contributions called compute_linearized_operators_basic,
//I create static function that calls that one from within. Adding the necessary
//contributions.
/**
Header file used for the definition of the [TFE](#thinFilm.h) Roughness Terms.

First we define the roughness $s$, its gradient $\mathbf \nabla s$ and
Laplacian $\Delta s$.
*/

face vector gs[];
scalar s[];
scalar lapS[];


/**
Macros for the computation of the free-surface positions are provided
*/

#define z ( h[] + s[] )
#define lapZ ( lapH[] + lapS[] )

/**
Macros that compute the roughness Jacobians
*/
#define roughLapJTerm ( pow_der(hf, 3.0) * face_gradient_x(lapS, 0) )

#if INCLINED
#define roughInclJTerm ( pow_der(hf, 3.0) * grav1 * gs.x[] )
#endif
