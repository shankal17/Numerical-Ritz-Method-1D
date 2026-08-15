using Weave

for dir in filter(isdir, readdir("scripts"; join=true))
    jmd = joinpath(dir, "report.jmd")
    isfile(jmd) || continue
    weave(jmd; doctype="md2html", out_path=joinpath(dir, "report"))
    @info "Wove $jmd"
end
