# ODEInterface

## (Planned) Version Numbers

Here is a picture, describing the (planned) version numbers:

```
 v0.1.4  v0.4.0  ...    [branch: master] (for current/nightly)
───┴─────┬─┴──────┴───────────────────────────────────────────────────────
         │
         │ [branch: julia_v0.4-v0.6]           [branch: julia_v0.4]
         └───┬──────┬────────────────┬────────────┬──────────┬────────────
           ...   v0.1.5              │        v0.1.6        ...
                                     │
                                     │
                                     │         [branch: julia_v0.5] 
                                     └┬───┬──────┬────────────────────────
                                      │ v0.2.0  ...  
                                      │
                                      │
                                      │        [branch: julia_v0.6]
                                      └───┬──────┬────────────────────────
                                       v0.3.0   ...

────────────────────────────────────────────────────────────────────> time
```

* v0.4.0, v0.4.1 are the new versions (with new features) in 
  the master branch for julia nightly.
* Versions until v0.1.5 are the versions 
  where new features (of the master branch) are backported
  for 0.4 ≤ julia ≤ 0.6.
  v0.1.5 is the last "feature" version for julia 0.4.
* v0.1.6, etc. are bug fixes (for the julia 0.4 version).
  Unfortunately this versions cannot be found in the Julia Metadata.
* The v0.2.* versions were planned for new features (of the master branch) 
  that are backported for 0.5 ≤ julia ≤ 0.6. 
  But the end of the supoort for julia 0.5 came first. Hence:
  v0.1.5 is the last "feature" version for julia 0.5.
* v0.2.0, etc are only bug fixes (for the julia 0.5 version).
  Unfortunately this versions cannot be found in the Julia Metadata.
* v0.3.0, v0.3.1 are the versions
  where new features (of the master branch) are backported
  for julia 0.6.

