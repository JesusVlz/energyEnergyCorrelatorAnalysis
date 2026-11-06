import pathlib
import subprocess
import tempfile
import unittest


class ClosureMatchingTest(unittest.TestCase):
    def test_cpp_contract(self):
        source = pathlib.Path(__file__).with_suffix('.cxx')
        with tempfile.TemporaryDirectory(prefix='oo_closure_matching_') as temporary:
            binary = pathlib.Path(temporary) / 'matching'
            subprocess.run(['g++', '-std=c++17', '-O2', '-Wall', '-Wextra', '-pedantic',
                            str(source), '-o', str(binary)], check=True)
            subprocess.run([str(binary)], check=True)


if __name__ == '__main__':
    unittest.main()
