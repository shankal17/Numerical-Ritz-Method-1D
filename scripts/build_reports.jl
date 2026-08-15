using Weave

const PDF_TEMPLATE = joinpath(@__DIR__, "pdf_template.tpl")

for dir in filter(isdir, readdir("scripts"; join=true))
    jmd = joinpath(dir, "report.jmd")
    isfile(jmd) || continue
    weave(jmd; doctype="md2html", out_path=joinpath(dir, "report"))
    weave(jmd; doctype="md2pdf", out_path=joinpath(dir, "report"), template=PDF_TEMPLATE)
    @info "Wove $jmd"
end
