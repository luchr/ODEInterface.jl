# ODEInterface

## (Planned) Version Numbers

Here is a picture, describing the (planned) version numbers:

```
 v0.1.4  v0.4.0         [branch: master] (for current/nightly)
───┴─────┬─┴──────────────────────────────────────────────────────────────
         │
         │              [branch: julia_v0.4-v0.6] 
         └───┬──────┬─────┬──────┬────────────┬──────────┬────────────────
           v0.1.5  ...   v0.1.a  │        v0.1.(a+1)    ...
                                 │
                                 │   [branch: julia_v0.5-v0.6] 
                                 └───┬──────┬────┬────┬────┬──────────┬───
                                   v0.2.0  ... v0.2.b │ v0.2.(b+1)   ...
                                                      │
                                                      │  [branch: julia_v0.6]
                                                      └───┬──────┬────────
                                                       v0.3.0   ...

──────────────────────────────────────────────────────────────────────> time
```

* v0.4.0, v0.4.1 are the new versions (with new features) in 
  the master branch for julia nightly.
* v0.1.5, v0.1.6 until v0.1.a are the versions 
  where new features (of the master branch) are backported
  for 0.4 ≤ julia ≤ 0.6.
  v0.1.a is the last "feature" version for julia 0.4.
* v0.1.(a+1), etc. are only bug fixes (for the julia 0.4 version).
* v0.2.0, v0.2.1 until v0.2.b are the versions
  where new features (of the master branch) are backported
  for 0.5 ≤ julia ≤ 0.6.
  v0.2.b is the last "feature" version for julia 0.5.
* v0.2.(b+1), etc are only bug fixes (for the julia 0.5 version).
* v0.3.0, v0.3.1 are the versions
  where new features (of the master branch) are backported
  for julia 0.6.

