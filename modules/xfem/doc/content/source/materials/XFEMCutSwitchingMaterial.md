# XFEMCutSwitchingMaterial

!syntax description /Materials/XFEMCutSwitchingMaterialReal

To allow using XFEM to model interfaces between materials, the XFEM system uses cut subdomain IDs to denote the subset of a standard MOOSE subdomain (element block) that a material point belongs to. Based on these cut subdomains, `XFEMCutSwitchingMaterial` switches between material properties depending on which cut subdomain a point is located within. This class reads in multiple versions of a material property, all with different values of `base_name`, and sets the value of the material property without that base name to the name of the property that applies to the current cut subdomain.

This class has four instantiations described below that are used to switch between properties with specific types.
Automatic-differentiation (AD) counterparts are provided for each instantiation; they mirror the
behavior of the non-AD materials but expose AD material properties so that derivatives are computed
automatically when the material is used with AD-capable kernels or other objects.

## XFEMCutSwitchingMaterialReal

`XFEMCutSwitchingMaterialReal` switches between two material properties of type `Real` that coexist on the same element with different base_names.

## Example Input File Syntax

!syntax parameters /Materials/XFEMCutSwitchingMaterialReal

!syntax inputs /Materials/XFEMCutSwitchingMaterialReal

## XFEMCutSwitchingMaterialRankTwoTensor

`XFEMCutSwitchingMaterialRankTwoTensor` switches between two material properties of type `RankTwoTensor` that coexist on the same element with different base_names.

## Example Input File Syntax

!syntax parameters /Materials/XFEMCutSwitchingMaterialRankTwoTensor

## XFEMCutSwitchingMaterialRankThreeTensor

`XFEMCutSwitchingMaterialRankThreeTensor` switches between two material properties of type `RankThreeTensor` that coexist on the same element with different base_names.

## Example Input File Syntax

!syntax parameters /Materials/XFEMCutSwitchingMaterialRankThreeTensor

## XFEMCutSwitchingMaterialRankFourTensor

`XFEMCutSwitchingMaterialRankFourTensor` switches between two material properties of type `RankFourTensor` that coexist on the same element with different base_names.

## Example Input File Syntax

!syntax parameters /Materials/XFEMCutSwitchingMaterialRankFourTensor

## ADXFEMCutSwitchingMaterialReal

!syntax description /Materials/ADXFEMCutSwitchingMaterialReal

`ADXFEMCutSwitchingMaterialReal` is the AD equivalent of
`XFEMCutSwitchingMaterialReal`, exposing the switched `Real` material property as
an AD value.

!syntax parameters /Materials/ADXFEMCutSwitchingMaterialReal

## ADXFEMCutSwitchingMaterialRankTwoTensor

!syntax description /Materials/ADXFEMCutSwitchingMaterialRankTwoTensor

`ADXFEMCutSwitchingMaterialRankTwoTensor` is the AD counterpart for
`XFEMCutSwitchingMaterialRankTwoTensor`, switching between `RankTwoTensor`
material properties while preserving derivative information.

!syntax parameters /Materials/ADXFEMCutSwitchingMaterialRankTwoTensor

## ADXFEMCutSwitchingMaterialRankThreeTensor

!syntax description /Materials/ADXFEMCutSwitchingMaterialRankThreeTensor

`ADXFEMCutSwitchingMaterialRankThreeTensor` provides AD material property
switching for `RankThreeTensor` quantities.

!syntax parameters /Materials/ADXFEMCutSwitchingMaterialRankThreeTensor

## ADXFEMCutSwitchingMaterialRankFourTensor

!syntax description /Materials/ADXFEMCutSwitchingMaterialRankFourTensor

`ADXFEMCutSwitchingMaterialRankFourTensor` offers the AD equivalent for
`RankFourTensor` properties, enabling automatic differentiation support when
switching between the underlying base-name properties.

!syntax parameters /Materials/ADXFEMCutSwitchingMaterialRankFourTensor
