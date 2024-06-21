cd("notes") do

#load macros from latex_macros file
macros = open("latex_macros.jl", "r") do f
    macros = read(f,String)
    macros = split(macros,"```")[2][6:end-1] #remove math and linebreaks
    macros = "md\"\"\"\n \$" * macros * "\$\n\"\"\""
    macros = replace(macros, "\n\n" => "\n") # remove empty lines
end

#load toc (to put in sidebar) from sidebar.md
toc = open("sidebar.md", "r") do f
    "md\"\"\"\n" * read(f,String) * "\n\"\"\""
end

for filepath in readdir()

    #exceptions and file loading
    if filepath[end-2:end] != ".jl"
        continue
    end
    
    file = read(filepath,String)

    if file[1:27] != "### A Pluto.jl notebook ###"
        println("$filepath is not a Pluto notebook\n")
        continue
    end

    do_write = false

    # latex macros
    macros_delimeter = "# latex macros --- DO NOT TOUCH THIS LINE"
    file_split = split(file,macros_delimeter)
    if length(file_split) == 3
        file_split[2] = "\n" * macros * " \n"
        file = join(file_split,macros_delimeter)
        println("$filepath : Latex macros added")
        do_write = true
    else 
        println("$filepath : No latex macros delimeted")
    end
    
    # Sidebar
    sidebar_delimeter = "# sidebar --- DO NOT TOUCH THIS LINE"
    file_split = split(file,sidebar_delimeter)
    if length(file_split) == 3
        file_split[2] = "\n" * toc * " \n"
        file = join(file_split,sidebar_delimeter)
        println("$filepath : Sidebar added")
        do_write = true
    else 
        println("$filepath : No sidebar (course TOC) delimited")
    end
    
    #write new file
    if do_write 
        open(filepath, "w") do f
            write(f,file)
            println("$filepath written")
        end
    end 

    println()

end

println("---")

end
