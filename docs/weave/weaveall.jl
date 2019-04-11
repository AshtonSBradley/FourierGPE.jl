using Weave
cssfile = joinpath(@__DIR__,"templates","skeleton_css.css")
latexfile = joinpath(@__DIR__, "..", "templates", "julia_tex.tpl")

source_path = @__DIR__
html_path = joinpath(source_path,"..","html")
pdf_path = joinpath(source_path,"..","pdf")
notebook_path = joinpath(source_path,"..","notebooks")

source_list = ["1dharmonic"]

function weavedocs(file)
    tmp = joinpath(source_path,file)*".jmd"
    weave(tmp,out_path=html_path,doctype = "md2html")
    weave(tmp,out_path=pdf_path,doctype="md2pdf"; template=latexfile)
    Weave.convert_doc(tmp,joinpath(notebook_path,file)*".ipynb")
end

for file ∈ source_list
    weavedocs(file)
end
