# boundary conditions for the auxiliary problem of the reentry problem
# following Stoer/Bulirsch 1973

# authors: Folkmar Bornemann, Vishal Sontakke, 2016/04/23

function aux_bc(xa,xb,r)
    if dir == "forward"
        aux_bc_(xa,xb,r);
    elseif dir == "backward"
        aux_bc_(xb,xa,r);
    end
    return nothing
end

function aux_bc_(xa,xb,r)
    r[1] = (xa[1]-v0);
    r[2] = xa[2]-gamma0;
    r[3] = xa[3]-h0;

    r[4] = (xb[1]-v1);
    r[5] = xb[2]-gamma1;
    r[6] = xb[3]-h1;
    return nothing
end
