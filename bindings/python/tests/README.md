### Pytest for the python wheel of NBIS-RS Library

### Setup
First generate the python wheels by running the following command in the root directory:

```bash
./build_python.sh
```

This will create a virtual environment `.venv` and the python wheel in the `dist` folder.

Next, in the root directory, activate the virtual environment and install the wheel. Modify the wheel name if necessary:

```bash
source ./.venv/bin/activate
pip install ./dist/nbis_py-0.1.2-py3-none-manylinux_2_34_x86_64.whl
```

Next, change the directory to the test folder:
```bash
cd bindings/python/tests
```

### Run the pytests

Install the pytest and other dependencies as:
```bash
pip install -r requirements.txt
```

And run the test:
```bash
pytest test_nbis.py 
```