[ This file was auto-generated from the module's documentation included in the doc-strings. Use julia's help system to get these informations in a nicer output format. ]

## The `opt` Argument: the Options

All options are handled by `OptionsODE`. See `OptionsODE` how to query, set and change options.

There are the following classes of options.

1. Options for this ODEInterface (common for all solvers)
2. Options for the ODE solvers

### Options of this ODEInterface

  * `OPT_RHS_CALLMODE`: There are two possible ways to call the Julia right-hand side: `RHS_CALL_RETURNS_ARRAY` and `RHS_CALL_INSITU`, see `help_callsolvers` for an explanation. difference.
  * `OPT_LOGIO`: This option sets the `IO` that is used for logging messages
  * `OPT_LOGLEVEL`: This is a bitmask for activating different logging messages. The following bitmasks are available.

    ```
       LOG_NOTHING     log nothing
       LOG_GENERAL     log some general information,
                       especially the main julia call of the solver
       LOG_RHS         log all calls of the right-hand side
       LOG_SOLVERARGS  log the arguments for the C-/Fortran-calls
                       before and after the call
       LOG_OUTPUTFCN   log calls of the julia output function
       LOG_SOLOUT      log calls of the solution output routine
       LOG_EVALSOL     log calls of the eval_sol_fcn
       LOG_MASS        log call(s) of the mass function
       LOG_JAC         log calls of the jacobian function of RHS
       LOG_BC          log calls of the boundary condition function
       LOG_BVPIVPSOL   log (during boundary value problems) calls to
                       initial value solvers
       LOG_RHSDT       log calls of the right-hand side time-derivative
       LOG_JACBC       log calls of the jacobian of the boundary condition
       LOG_GUESS       log calls to the guess function
       LOG_ALL         all of the above
    ```

### Options for the solvers

Different solvers support different options. All the options a solver supports are listed in the help-command of the specific solver, e.g. `help_dopri5`.

To get an overview, what options are supported by what solvers, call `ODEInterface.help_options()`.




## Solvers

1. bvpm2
2. bvpsol
3. colnew
4. ddeabm
5. ddebdf
6. dop853
7. dopri5
8. odex
9. radau
10. radau5
11. rodas
12. seulex

## Options used by Solvers

<table>
<tr><th><pre>        Option
</pre></th>
<th><pre> 1
</pre></th>
<th><pre> 2
</pre></th>
<th><pre> 3
</pre></th>
<th><pre> 4
</pre></th>
<th><pre> 5
</pre></th>
<th><pre> 6
</pre></th>
<th><pre> 7
</pre></th>
<th><pre> 8
</pre></th>
<th><pre> 9
</pre></th>
<th><pre>10
</pre></th>
<th><pre>11
</pre></th>
<th><pre>12
</pre></th>
</tr>
<tr><td><pre> OPT&#95;ADDGRIDPOINTS
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;ATOL
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;BVPCLASS
</pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;COARSEGUESSGRID
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;COLLOCATIONPTS
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;DENSEOUTPUTWOEE
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;DIAGNOSTICOUTPUT
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;DIMOFIND1VAR
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;DIMOFIND2VAR
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;DIMOFIND3VAR
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;EPS
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;ERRORCONTROL
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;FREEZEINTERVALS
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;FREEZESSLEFT
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;FREEZESSRIGHT
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;INITIALSS
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;INITSTAGES
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;INTERPOLDEGREE
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;IVPOPT
</pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;JACOBIBANDSTRUCT
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;JACOBIMATRIX
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;JACRECOMPFACTOR
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;LAMBDADENSE
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;M1
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;M2
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;MASSMATRIX
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;MAXEXCOLUMN
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;MAXNEWTONITER
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;MAXS
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;MAXSS
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;MAXSTABCHECKLINE
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;MAXSTABCHECKS
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;MAXSTAGES
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;MAXSTEPS
</pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;MAXSUBINTERVALS
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;METHODCHOICE
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;MINSTAGES
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;NEWTONSTARTZERO
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;NEWTONSTOPCRIT
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;ORDERDECFACTOR
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;ORDERDECFRAC
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;ORDERDECSTEPFAC1
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;ORDERDECSTEPFAC2
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;ORDERINCFACTOR
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;ORDERINCFRAC
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;OUTPUTATTIMES
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;OUTPUTFCN
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;OUTPUTMODE
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;RHO
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;RHO2
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;RHSAUTONOMOUS
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;RHSTIMEDERIV
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;RHS&#95;CALLMODE
</pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;RTOL
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;SINGULARTERM
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;SOLMETHOD
</pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;SSBETA
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;SSMAXSEL
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;SSMINSEL
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;SSREDUCTION
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;SSSELECTPAR1
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;SSSELECTPAR2
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;STEPSIZESEQUENCE
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;STEPSIZESTRATEGY
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;STEST
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;SUBINTERVALS
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;TRANSJTOH
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;TSTOP
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre> &#10004;
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
<tr><td><pre> OPT&#95;WORKFORDEC
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;WORKFORJAC
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;WORKFORRHS
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre> OPT&#95;WORKFORSOL
</pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre> &#10004;
</pre></td>
</tr>
<tr><td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
<td><pre></pre></td>
</tr>
</table>


