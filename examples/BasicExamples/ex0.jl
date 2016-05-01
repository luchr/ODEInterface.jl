using ODEInterface

# try to load shared libraries with solvers
ODEInterface.loadODESolvers()

# give summary
display(ODEInterface.help_solversupport())

# vim:syn=julia:cc=79:fdm=indent:
