# Experimental modules

The modules below run, but no validated reference comparison exists for them yet. Each
entry point calls `cmpmc::require_validated`, which throws `cmpmc::UnvalidatedModuleError`
unless you opt in:

```cpp
cmpmc::DetectorModel detector;
detector.allow_unvalidated = true;   // for benchmarking and development only
```

Results from these modules should not be used as physics results. Their known defects
and the tests each one needs before it can be promoted are listed in
[limitations.md](limitations.md).

| Module | Location | How to enable |
|---|---|---|
| SoA wavefront engine | `cmpmc/experimental/wavefront_soa.hpp`, `experimental::WavefrontSoAEngine` | construct it with `allow_unvalidated` |
| Queue-based wavefront pipeline | `cmpmc/experimental/wavefront_pipeline.hpp`, `experimental::WavefrontPipelineEngine` | construct it with `allow_unvalidated` |
| Jump-to-ballistic down-conversion | `cmpmc/experimental/jump_to_ballistic.hpp`, `DetectorModel::jump_to_ballistic.enabled` | `allow_unvalidated` |
| Next-event estimation | `cmpmc/experimental/next_event_estimation.hpp`, `experimental::evaluate_sparse_sensor_nee` | pass `allow_unvalidated = true` as the last argument |
| Bulk 3D quasiparticle transport | `SimulationEngine::transport_superconductor_coupled` | `allow_unvalidated` |
| Intervalley Matrix rate | `CrystalConfig::iv_model = IVRateModel::Matrix` | `allow_unvalidated` |
| Polynomial surface coefficients | `SurfaceProperty::demo_polynomial_sensor_surface()`, `demo_polynomial_bare_surface()` | `allow_unvalidated` |

The umbrella header `cmpmc/cmpmc.hpp` does not include `cmpmc/experimental/*`. Include
those headers explicitly.
