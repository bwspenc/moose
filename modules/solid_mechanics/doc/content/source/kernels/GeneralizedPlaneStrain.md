# Generalized Plane Strain Kernel

!syntax description /Kernels/GeneralizedPlaneStrain

The `GeneralizedPlaneStrain` kernel assembles the out-of-plane equilibrium
residual for the scalar out-of-plane strain variable, together with the
off-diagonal Jacobian coupling that variable to its own attached in-plane
displacement (or temperature) variable. This object is usually set up by the
[GeneralizedPlaneStrainAction](SolidMechanics/GeneralizedPlaneStrain/index.md),
which creates one instance per in-plane displacement (plus one for
temperature, if coupled); exactly one of those instances (the "primary" one)
additionally owns the scalar variable's own residual and diagonal Jacobian
contributions, since that data must only be assembled once. See
[ADGeneralizedPlaneStrain.md] for the automatic-differentiation counterpart of
this kernel, which needs only a single instance because AD captures all of
the coupling automatically.

The equilibrium condition when the out-of-plane direction is the $x$-direction is given as
\begin{equation}
	\int_{A}{\sigma_{xx}dA} = \bar{N}_{xx}
\end{equation}

The equilibrium condition when the out-of-plane direction is the $y$-direction is given as
\begin{equation}
	\int_{A}{\sigma_{yy}dA} = \bar{N}_{yy}
\end{equation}

The equilibrium condition when the out-of-plane direction is the $z$-direction is given as
\begin{equation}
	\int_{A}{\sigma_{zz}dA} = \bar{N}_{zz}
\end{equation}


A detailed description of generalized plane strain formulation can be found in [here](solid_mechanics/generalized_plane_strain.md).

!syntax parameters /Kernels/GeneralizedPlaneStrain

!syntax inputs /Kernels/GeneralizedPlaneStrain

!syntax children /Kernels/GeneralizedPlaneStrain

!bibtex bibliography
