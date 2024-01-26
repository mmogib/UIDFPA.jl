using UIDFPA
using Documenter

DocMeta.setdocmeta!(UIDFPA, :DocTestSetup, :(using UIDFPA); recursive=true)

makedocs(;
    modules=[UIDFPA],
    authors="Mohammed Alshahrani <mmogib@gmail.com> and contributors",
    sitename="UIDFPA.jl",
    format=Documenter.HTML(;
        canonical="https://mmogib.github.io/UIDFPA.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mmogib/UIDFPA.jl",
    devbranch="master",
)
