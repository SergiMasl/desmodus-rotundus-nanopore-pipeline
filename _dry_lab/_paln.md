# HPC:

## Step 1: Copy the project files

The project files are in the `bat_spleens` directory. Copy them to your personal workspace on the HPC:

```
cp -r /projects/class/binf6350_001/bat_spleens /users/<user_id>/
```

Replace `<user_id>` with your HPC user ID.

**Output:** You have a personal copy of the project to work in, leaving the shared source files unchanged.

## Step 2: Read the first-stage instructions

Read **Stage A1: First Contact (Gate).md**. It contains the instructions for the first stage of the workflow; complete that stage before moving on.

## Poly-A Tail Trimming:

Do not use Porechop, and cutadapt --poly-a is not the right call for our particular project. We are aiming to find the fraction of reads that span a complete transcript, so using 
minimap2 later in the computational pipeline is a far better option. 
