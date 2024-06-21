cd("src") do

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
        file_split[2] = "\ninclude(\"latex_macros.jl\") \n"
        file = join(file_split,macros_delimeter)
        println("$filepath : Latex macros changed to include(latex_macros.jl)")
        do_write = true
    else 
        println("$filepath : No latex macros delimeted")
    end
    
    # Sidebar
    sidebar_delimeter = "# sidebar --- DO NOT TOUCH THIS LINE"
    file_split = split(file,sidebar_delimeter)
    if length(file_split) == 3
        file_split[2] = "\nMarkdown.parse( \"**Error control in scientific modeling** \n\" * read(\"sidebar.md\",String)) \n"
        file = join(file_split,sidebar_delimeter)
        println("$filepath : toc changed to read(\"sidebar.md\",String)")
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
