# elmer-physics-meshtools

The mesh utility solvers of Elmer's `fem/src/modules` for the Rust-native
Elmer rewrite (lane `w6mappers`), as `Solver` implementations against
`elmer-solver`. `register(&mut registry)` binds all four.

| Procedure | Fortran | Does |
|---|---|---|
| `MeshChecksum` | MeshChecksum.F90 | eight mesh sums and a pseudo-norm against `Reference Values` |
| `RigidMeshMapper` | RigidMeshMapper.F90 | rotate / scale / translate bodies (`Mesh Rotate`, `Mesh Scale`, `Mesh Translate`, `Mesh Origin`, `Mesh Matrix`), weighted by `Mesh Relax` or a Laplace relaxation field |
| `StructuredMeshMapper` | StructuredMeshMapper.F90 | map an extruded mesh between `Top Surface` / `Bottom Surface` (and `Mid Surface`), or between fixed layers |
| `StructuredProjectToPlane` | StructuredProjectToPlane.F90 | reduce fields of an extruded mesh to its top plane, or compute `depth`, `height`, `distance`, `index` |

`detect_extruded_structure` is `DetectExtrudedStructure` (MeshTransform.F90).

## Moving the mesh

The Fortran mappers assign `Mesh % Nodes % x/y/z` in place. A Rust solver sees
the mesh as `&dyn MeshView` and cannot. `motion::MeshMotion` is the one thing
the mappers need from the runner: install it for the simulation's thread with
`motion::install(Some(Arc<dyn MeshMotion>))`, and make `MeshView::coord` read
the same coordinates (`motion::SharedCoords` is a ready-made store). Without
it, a mapper that would move a node fails with `NotWired`. When `elmer-core`
grows `MeshView::set_coords` this seam can go.

## Verified against the pristine ElmerSolver

Mapped coordinates are compared with the pristine run's VTU points, solved
fields with its VTU point data, and norms with the `.sif`'s `Reference Norm`.

| Quick test | Solvers | Checked |
|---|---|---|
| `ElmerGridConvertUnv` | MeshChecksum | all eight sums bit for bit, pseudo-norm 1 |
| `RigidMeshMapper1` | RigidMeshMapper + HeatSolve | nodes 1e-12, temperature, `Reference Norm` |
| `RigidMeshMapper2` | RigidMeshMapper (Laplace relaxation) + HeatSolve | `MeshRelax`, nodes, `Reference Norm` |
| `structmap` | StructuredMeshMapper + HeatSolve | nodes 1e-12, temperature, `Reference Norm` |
| `structmap5` | StructuredMeshMapper (mid surface) + StructuredProjectToPlane + HeatSolve | nodes, eleven projected fields 1e-12, both `Reference Norm`s |

The test meshes are the ElmerGrid conversions of the pristine tests' inputs.

## Not ported (refused with `Unsupported`)

- `RigidMeshMapper`: `Rotor Mode`, `Calculate Mesh Velocity`, `Target Mesh
  Solver Index`, `mmg remesh`, `Mesh Youngs Modulus`.
- `StructuredMeshMapper`, `StructuredProjectToPlane`: `Mapping Mask
  Variable` / `Mapping Mask Name`, `Mapping Coordinate Transformation`,
  `Mapping Original Coordinates`, `Base Displacement Variable`, `Mapped Mesh
  Name`, `Mesh Height Map`, `Correct Surface Mask`, the `error ...` operators
  and discontinuous Galerkin input variables.

Quick tests that need these (`structmap6`, `structmap7`, `structmap8`) and the
other conversions behind `ElmerGridConvert{Comsol,Fvcom,GmshBinary}` (which
need the `elmer-mesh` readers) are not covered yet.
