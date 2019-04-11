using Weave
cssfile = joinpath(@__DIR__,"docs/weave", "skeleton_css.css")
weave(joinpath(@__DIR__,"docs/weave/1dharmonic.jmd"),out_path="./docs/html",doctype = "md2html";css=cssfile)
