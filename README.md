# ParallelFlow

15418 Spring 2025 Project
https://github.com/garygao33/ParallelFlow

## Deliverables:

[Project Proposal](15418_Project_Proposal.pdf)

[Project Milestone Report](15418_Project_Milestone.pdf)

[Project Report](15418_Project_Report.pdf)

[Project Presentation](15418_Project_Slides.pdf)

## How to run the code:

We kept some legacy versions of our implementations. Before compiling, remove every file except `parallel.cpp` (or any version you would like to keep) from the `/src` directory

Run `make` in the root directory. The run the program by `./flow -n num_threads -m mode < path/to/input`.

## More on inputs:

Due to restrictions on file sizes, we are unable to upload the larger inputs (such as ones that should be in `inputs/sparse`) to this repository. If you are reading this from Autolab submission, we had to remove all inputs due to submission size limit. We only uploaded some basic tests. 

If you would like to generate inputs similar to the ones used in the report, you can use the generation scripts in `/inputs`. If you would like to obtain the exact same inputs, please contact the owner of this repository via email.