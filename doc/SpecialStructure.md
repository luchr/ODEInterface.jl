[ This file was auto-generated from the module's documentation included in the doc-strings. Use julia's help system to get these informations in a nicer output format. ]

## Special Structure

Some solvers (e.g. radau5 and radau) supports ODEs with a "special structure". In this context an ODE has special structure  if there exists M1>0 and M2>0 and M1 = M⋅M2 for some integer M and

```
             x'(k) = x(k+M2)   for all k = 1,…,M1            (*)
```

There are the options `OPT_M1` and `OPT_M2` to tell the solvers, if the ODE has this special structure. In this case only the non-trivial parts of the right-hand side, of the mass matrix and of the Jacobian matrix  had to be supplied.

If an ODE has the special structure, then the right-hand side  for `OPT_RHS_CALLMODE == RHS_CALL_RETURNS_ARRAY` has to return a vector of length d-M1 because the right-hand side for the first M1 entries are known from (*). For `OPT_RHS_CALLMODE == RHS_CALL_INSITU` the right-hand side gets (a reference) to the array of length d, but only the last d-M1 components need to be filled in.

The mass matrix has the form

```
     ⎛ 1    │    ⎞ ⎫
     ⎜  ⋱   │    ⎟ ⎪
     ⎜   ⋱  │    ⎟ ⎬ M1
     ⎜    ⋱ │    ⎟ ⎪
M =  ⎜     1│    ⎟ ⎭
     ⎜──────┼────⎟
     ⎜      │    ⎟ ⎫
     ⎜      │ M̃  ⎟ ⎬ d-M1
     ⎝      │    ⎠ ⎭
```

Then there has to be only the (d-M1)×(d-M1) matrix M̃ in `OPT_MASSMATRIX`. Of course, M̃ can be banded. Then `OPT_MASSMATRIX` is a `BandedMatrix` with lower bandwidth < d-M1.

If an ODE has the special structure, then the Jacobian matrix has the form

```
     ⎛0   1      ⎞ ⎫
     ⎜ ⋱   ⋱     ⎟ ⎪
     ⎜  ⋱   ⋱    ⎟ ⎬ M1
     ⎜   ⋱   ⋱   ⎟ ⎪
J =  ⎜    0   1  ⎟ ⎭
     ⎜───────────⎟
     ⎜           ⎟ ⎫
     ⎜      J̃    ⎟ ⎬ d-M1
     ⎝           ⎠ ⎭
```

Then there has to be only the (d-M1)×d matrix J̃ in `OPT_JACOBIMATRIX`. In this case banded Jacobian matrices are only supported for the case M1 + M2 == d. Then in this case J̃ is divided into d/M2 = 1+(M1/M2) blocks of the size M2×M2. All this blocks can be banded with a (common) lower bandwidth < M2. 

The option `OPT_JACOBIBANDSTRUCT` is used to describe the banded structure of the Jacobian. It is eigher `nothing` if the Jacobian is full or a tuple `(l,u)` with the lower and upper bandwidth.

The function for providing the Jacobian for ∂f/∂x can have the following forms:

```
 function (t,x,J)       -> nothing       (A)
 function (t,x,J1,…,JK) -> nothing       (B)
```

The following table shows when `OPT_JACOBIMATRIX` has the form (A)  and when it has the form (B):

<table>
<tr><th><pre></pre></th>
<th><pre> JACOBIBANDSTRUCT &#61;&#61; nothing
</pre></th>
<th><pre> JACOBIBANDSTRUCT &#61;&#61; &#40;l&#44;u&#41;
</pre></th>
</tr>
<tr><td><pre> M1&#61;&#61;0
</pre></td>
<td><pre> &#40;A&#41;&#44; J full&#44; size&#58; d&#215;d
</pre></td>
<td><pre> &#40;A&#41; J &#40;l&#44;u&#41;&#45;banded&#44; size d&#215;d
</pre></td>
</tr>
<tr><td><pre> M1&#62;0
</pre></td>
<td><pre> &#40;A&#41;&#44; J full&#44; size&#58; &#40;d&#45;M1&#41;&#215;d
</pre></td>
<td><pre> &#40;B&#41; J1&#44;&#8230;&#44;JK &#40;l&#44;u&#41;&#45;banded
     each with size M2&#215;M2 and
     K &#61; 1 &#43; M1&#47;M2
     M1 &#43; M2 &#61;&#61; d
</pre></td>
</tr>
</table>



