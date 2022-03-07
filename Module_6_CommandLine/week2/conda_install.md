# conda

Miniconda is the command line version. Anaconda is graphical but is too big for its own good.

These instructions are for **Latest - Conda 4.11.0 Python 3.9.7 released February 15, 2022.**

Full doc: https://docs.conda.io/en/latest/miniconda.html

## Mac

Do:

`curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh > ~/miniconda.sh`

`shasum -a 256  ~/miniconda.sh`

Does it say 7717253055e7c09339cd3d0815a0b1986b9138dcfcb8ec33b9733df32dd40eaa?

`bash ~/miniconda.sh -b -p $HOME/miniconda`

Follow the prompts. If you are unsure about any setting, accept the defaults. You can change them later.


The installer prompts “Do you wish the installer to initialize Miniconda3 by running conda init?”  **YES**

Full document: https://conda.io/projects/conda/en/latest/user-guide/install/macos.html

## Ubuntu

`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

`sha256sum  Miniconda3-latest-Linux-x86_64.sh`

Does it say 4ee9c3aa53329cd7a63b49877c0babb49b19b7e5af29807b793a76bdb1d362b4?

`bash Miniconda3-latest-Linux-x86_64.sh`

Full doc: https://conda.io/projects/conda/en/latest/user-guide/install/linux.html
