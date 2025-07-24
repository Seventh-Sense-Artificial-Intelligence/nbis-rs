### Pytest for the python wheel of NBIS-RS Library

### Setup
First generate the python wheels by running the following command in the root directory:

```bash
./build_python.sh
```

This will create a virtual environment `.venv` and the python wheel in the `dist` folder.

Now, in the root directory, activate the virtual environment and install the wheel:

```bash
source ./.venv/bin/activate
pip install ./dist/nbis_py-0.1.2-py3-none-manylinux_2_31_x86_64.whl
```