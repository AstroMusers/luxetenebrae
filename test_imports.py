# test_imports.py
# Run this script from the parent directory (where luxetenebrae/ is located)

from luxetenebrae.utils import utils, calculations

print("utils module loaded:", utils)
print("calculations module loaded:", calculations)

# Example: call a function if it exists
if hasattr(calculations, 'some_function'):
    print("Result from some_function:", calculations.some_function())
else:
    print("No 'some_function' in calculations module.")
