---------------------------------------------------------------------------
Instructions for running the bash-ified version of Tim Woods' python script
---------------------------------------------------------------------------

The working directory should contain the following files at a bare minimum:
Start.sh
move.sh
pipeToPython.sh
Pipeable.py

The Start.sh script calls the other two bash scripts in its execution. move.sh is
reponsible for moving the pertainent files to the "input" directory in the working
directory (which is created if it doesn't exist). The pipeToPython.sh script is
responsible for sorting the input files and passing them to Pipeable.py for
evaluation.

1) Place a .gz file containing the protein data  into the same directory as
the above files.

2) Execute the Start.sh file in the command line.

3) Relax and enjoy as the program executes...

4) Results will initally be placed in to the working directory, then zipped
and moved in to the "output" directory in batches of 5 output files. As a
result there may be some non-zipped directories left over in the working
directory after the program has run. This could be changed later...

5) The bash scripts remove the un-tar'd input data after running, but leave
the tar.gz file alone. 

