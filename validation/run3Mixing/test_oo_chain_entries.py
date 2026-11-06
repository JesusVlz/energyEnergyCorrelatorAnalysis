import pathlib
import shlex
import subprocess
import tempfile
import unittest


class MultiFileChainCount(unittest.TestCase):
    def test_native_generator_and_reco_total_count(self):
        repo = pathlib.Path(__file__).resolve().parents[2]
        flags = shlex.split(subprocess.check_output(['root-config', '--cflags', '--libs'], text=True))
        with tempfile.TemporaryDirectory(prefix='oo_chain_entries_') as temporary:
            binary = pathlib.Path(temporary)/'reader'
            subprocess.run(['g++', '-Isrc', str(pathlib.Path(__file__).with_suffix('.cxx')),
                            'src/ForestReader.cxx', 'src/HighForestReader.cxx',
                            'src/GeneratorLevelForestReader.cxx', *flags, '-lEG', '-o', str(binary)],
                           cwd=repo, check=True)
            subprocess.run([str(binary), str(pathlib.Path(temporary)/'fixtures')], cwd=repo, check=True)


if __name__ == '__main__':
    unittest.main()
