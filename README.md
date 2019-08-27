
<center> <h1> multicrispr </h1> </center>


![](https://gitlab.gwdg.de/loosolab/software/multicrispr/wikis/uploads/43b432cd32eb156af2ac217efd98aceb/workflow.png)


1. Read a set of genomics ranges 

2. Find cas9 sequences which 

    - target

        + **thousands** of such **ranges**, their **flanks**, or **slopped extensions**
        + **efficiently** with **data.table** and **Biostrings** based **c-level** looping and matching
        + enabling **crispr<sub>ko</sub>**, **crispr<sub>i</sub>**, and **crispr<sub>a</sub>** design.

   - are free of offtarget (mis)matches

   - bind well
   
        + as determined by Doench 2016 ontargetscore (interfacing to original python module **azimuth** when available)
        + or Doench 2014 ontargetscore otherwise
