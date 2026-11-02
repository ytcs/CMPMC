# Materials and inputs

## Built-in material presets

`#include <cmpmc/materials/presets.hpp>`, namespace `cmpmc::materials`:

| Function | Material | Crystal system |
|---|---|---|
| `Germanium()` | Ge | cubic |
| `Silicon()` | Si | cubic |
| `GalliumArsenide()` | GaAs | cubic |
| `CalciumFluoride()` | CaF₂ | cubic |
| `LithiumFluoride()` | LiF | cubic |
| `Sapphire()` | Al₂O₃ | rhombohedral |
| `AlphaQuartz()` | α-SiO₂ | rhombohedral |
| `CalciumTungstate()` | CaWO₄ | tetragonal |
| `AmorphousSilica()` | a-SiO₂ | amorphous |
| `AluminumCrystal()` | Al | cubic |
| `NiobiumCrystal()` | Nb | cubic |
| `Copper()` | Cu | cubic |

Thin superconducting films (`SuperconductorFilmConfig`): `AluminumFilm()` and
`TungstenTES()`.

The Si, Ge and GaAs phonon constants follow arXiv:2408.04732. The Si pair energy and
Fano factor are from arXiv:2004.11499. The primary-phonon frequency `debye_freq` is the
maximum optical-phonon frequency. [validation.md](validation.md) lists which constants
are independently checked.

## `CrystalConfig`

`CrystalConfig` (`materials/crystal_config.hpp`) is a plain struct. Every field has a
unit comment, and the defaults describe germanium. The field groups are:

| Group | Fields |
|---|---|
| Identity | `name`, `group` (`BravaisGroup`), `atomic_Z`, `atomic_A`, `density` |
| Elasticity | `C11 … C66`, or `set_stiffness(i, j, value)` for any 1-based Voigt component; `beta, gamma, lambda, mu` (third-order constants) |
| Phonon rates | `B_iso` (s³), `A_anh` (s⁴), `fTT_frac`, `polycrystal_mfp_m` |
| Density of states | `fLDOS`, `fSTDOS`, `fFTDOS` |
| Sound speeds | `v_sound`, `v_trans`, `debye_freq` |
| Band structure | `bandgap`, `pair_energy`, `fano_factor`, `epsilon_rel`, `h_mass_ratio`, `e_mass_tensor`, `valley_dirs` |
| Intervalley scattering | `iv_model` (`IVRateModel`), `iv_quad_*`, `iv_lin_*`, `iv_deform_eV_m`, `iv_energy_eV` |
| Trapping | `e_trapping_mfp_m`, `h_trapping_mfp_m`, `*_trap_ion_*_mfp_m`, trap energies |
| Bulk superconductor | `sc_delta0_eV`, `sc_tau0_qp_s`, `sc_tau0_ph_s`, `sc_dn_m2_s`, `sc_teff_K`, `sc_trap_tau_s` |
| Orientation | `orient` (a rotation matrix), or `set_miller_orientation(h, k, l)` |

Use unit constants from `cmpmc::units` when filling fields, for example
`c.C11 = 165.6 * units::GPa`.

## Importing G4CMP `config.txt` files

```cpp
#include <cmpmc/io/g4cmp_config_import.hpp>

cmpmc::CrystalConfig ge = cmpmc::io::import_g4cmp_config(
    "CrystalMaps/Ge/config.txt", /*density_kg_m3=*/5323.0, /*name=*/"Ge");
```

The importer reads G4CMP's whitespace-separated `key values [unit]` format (`#` starts a
comment). It recognises these keys: the crystal-system keywords (`cubic`, `tetragonal`,
…), `stiffness`, `dyn`, `scat`, `decay`, `decayTT`, `LDOS`, `STDOS`, `FTDOS`, `Debye`,
`vsound`, `vtrans`, `bandgap`, `pairEnergy`, `fanoFactor`, `hmass`, `emass`, `valleyDir`,
`acDeform_e`, `acDeform_h`, `sc_tau0_qp` and `sc_tau0_ph`. It ignores keys it does not
recognise. The file format does not carry the density, so pass it as an argument.
Every field not set by the file keeps its `CrystalConfig` default.

The `test_g4cmp_import` test checks the built-in presets against a G4CMP
`CrystalMaps` tree when you configure with
`-DCMPMC_G4CMP_CRYSTALMAPS_DIR=/path/to/CrystalMaps`.

## Electric-field meshes

`#include <cmpmc/geometry/mesh_field.hpp>`, namespace `cmpmc::geometry`.

- **`TriLinearMesh3D::load_from_files(coords, tets, coord_unit_m, potential_unit_V)`**
  takes two files:
  - the vertex file: one `x y z V` line per node
  - the tetrahedron file: one `i0 i1 i2 i3` line per cell (0-based node indices)

  Lines starting with `#` or `%` are ignored. Use the unit arguments to convert the
  coordinates to metres and the potential to volts.
- **`TriLinearMesh3D::build_cartesian_grid(...)`** and
  **`BiLinearMesh2D::build_grid(..., axisymmetric_rz)`** build meshes from a potential
  function, which is useful for analytic fields and tests.

Assign a mesh to `DetectorModel::mesh_field_3d` or `mesh_field_2d`. It takes precedence
over the uniform field.

## Electrode and sensor layouts

`geometry::ElectrodePatternMask` describes the active metal on a face as a list of
pads. A pad is a disk, an annular sector or a rectangle, with an optional fill fraction
and a channel id. Two presets are provided:
- `supercdms_qet_mask(R, fill)`: the SuperCDMS five-channel QET layout, with an inner disk, three 120° outer sectors, and a
  guard ring
- `rectangular_grid_mask(nx, ny, ...)`: an nx × ny grid of pads

Attach a mask to `SurfaceProperty::electrode_mask`. A hit on a pad records its channel
in `SensorHit::channel_id`. A hit between pads reflects as from a bare surface.
