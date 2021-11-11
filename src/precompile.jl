function capture_stdout(f)
    #https://discourse.julialang.org/t/consistent-way-to-suppress-solver-output/20437/6
    stdout_orig = stdout
    (rd, wr) = redirect_stdout()
    f()
    close(wr)
    redirect_stdout(stdout_orig)
    read(rd, String)
end

function __precompile__()
    System = SquareLattice.getSquareLattice(4,[1,0.5,0.2],test=true)
    System = Pyrochlore.getPyrochlore(4,[1,0.5,0.2],test=true)
    return
end

function __precompile__quiet__()
    capture_stdout(__precompile__)
    return
end