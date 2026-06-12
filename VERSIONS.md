# ODEInterface

## (Planned) Version Numbers


| Julia Version |  ODEInterface Versions |
| ------------: | ---------------------- |
| v0.4          | `v0.1.*`               |
| v0.5          | `v0.2.*`               |
| v0.6          | `v0.3.*`               |
| v0.7          | `v0.4.*`               |
| ≥ v1.0        | `v0.5.*`               |

Here is a picture, describing the (planned) version numbers:

```
 v0.1.4  v0.4.0  ...  v0.4.8    v0.5.0               [branch: master]
───┴─────┬─┴──────┴──────┴───┬─────┴─────────────────────────────────
         │                   │
         │                   │                   [branch: julia_v0.7]
         │                   └─────┬────────┬────────────────────────
         │                      v0.4.9     ...
         │
         │
         │                                       [branch: julia_v0.6]
         │                            ┌────┬──────┬──────────────────
         │                            │ v0.3.0   ...
         │                            │
         │                            │
         │                            │          [branch: julia_v0.5]
         │                           ┌┴───┬──────┬───────────────────
         │                           │  v0.2.0  ...
         │                           │
         │                           │
         │ [branch: julia_v0.4-v0.6] │           [branch: julia_v0.4]
         └───┬──────┬────────────────┴────┬──────┬───────────────────
           ...   v0.1.5                 v0.1.6  ...

```

* julia ≥ 1.0: v0.5.0, v0.5.1 are the versions (with new features) 
  in the master branch
* julia 0.7: v0.4.8 is the last "feature" version; v0.4.9, ... only bugfixes
* julia 0.6: v0.3.1 is the last "feature" version; v0.3.2, ... only bugfixes
  Unfortunately this versions cannot be found in the Julia Metadata.
* julia 0.5: v0.1.5 is the last "feature" version; v0.2.0, ... only bugfixes
  Unfortunately this versions cannot be found in the Julia Metadata.
* julia 0.4: v0.1.5 is the last "feature" version; v0.1.6, ... only bugfixes
  Unfortunately this versions cannot be found in the Julia Metadata.

