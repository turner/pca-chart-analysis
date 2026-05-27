# Conda Cheat Sheet

## Environments

```bash
conda env list                              # list all envs
conda create -n NAME python=3.11            # new env with Python
conda create -n NAME python=3.11 numpy pandas matplotlib   # with packages
conda activate NAME                         # switch in
conda deactivate                            # exit current
conda env remove -n NAME                    # delete env
conda rename -n OLD NEW                     # rename env
```

## Packages

```bash
conda install PKG                           # install (current env)
conda install PKG=1.2.3                     # pinned version
conda install -c conda-forge PKG            # from conda-forge channel
conda update PKG                            # upgrade one
conda update --all                          # upgrade everything
conda remove PKG                            # uninstall
conda list                                  # what's installed
conda list PKG                              # check one package
conda search PKG                            # available versions
```

## pip inside conda
When a package isn't on conda, pip is fine — but install conda packages first.
```bash
conda install pip                           # ensure pip is in this env
pip install PKG                             # then pip-install
```

## Export / reproduce

```bash
conda env export > environment.yml          # full lockfile (platform-specific)
conda env export --from-history > env.yml   # only what YOU asked for (portable)
conda env create -f environment.yml         # recreate from file
conda env update -f environment.yml         # apply changes
```

## Channels

```bash
conda config --show channels                # current priority order
conda config --add channels conda-forge     # prepend channel
conda config --set channel_priority strict  # recommended for conda-forge
```

## Info / cleanup

```bash
conda info                                  # versions, paths, config
conda info --envs                           # same as `conda env list`
conda clean --all                           # purge caches (saves GB)
which python ; python -V                    # confirm active env
```

## Common gotchas

- **Always `conda activate NAME` before `pip install`** — otherwise it installs into base.
- **`conda-forge` + `defaults` mixed** can cause solver pain. Pick one (usually `conda-forge`) and set `channel_priority: strict`.
- **`environment.yml` from `--from-history`** is portable across OSes; the default export pins platform-specific builds and only works on the same OS.
- **`base` env is for conda itself** — don't install project packages there.
- **`mamba` is a drop-in faster replacement** for `conda install/create/update`: `mamba install PKG`.

## One-liner: new project env

```bash
conda create -n myproj -c conda-forge python=3.11 numpy pandas matplotlib scikit-image pillow tqdm && conda activate myproj
```
