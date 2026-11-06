"""Exercise plotting build invalidation without requiring a ROOT installation."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import textwrap
import time
import unittest


REPOSITORY = Path(__file__).resolve().parents[2]


class PlottingBuildTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="eec-plotting-build-")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        shutil.copyfile(REPOSITORY / "plotting/Makefile", self.directory / "Makefile")
        self.tools = self.directory / "root/bin"
        self.tools.mkdir(parents=True)
        self.path_tools = self.directory / "path-bin"
        self.path_tools.mkdir()
        self.log = self.directory / "events.jsonl"
        self.identity = self.directory / "root-version"
        self.identity.write_text("6.32.fixture")
        self.environment = dict(os.environ, PATH=str(self.path_tools) + os.pathsep + os.environ["PATH"])
        self.environment["BUILD_EVENT_LOG"] = str(self.log)
        self.environment["FAKE_ROOT_VERSION"] = str(self.identity)
        for name in ("SplitCanvas.h", "SplitCanvas.cxx", "Helper.h", "Helper.cxx", "Common.h", "LinkDef.h"):
            path = self.directory / name
            path.write_text("// fixture source\n")
            os.utime(path, (time.time() - 100, time.time() - 100))
        self.write_tool(self.tools / "root-config", """
            version = Path(os.environ["FAKE_ROOT_VERSION"]).read_text()
            prefix = Path(__file__).resolve().parents[1]
            answers = {
                "--version": version, "--prefix": str(prefix),
                "--bindir": str(prefix / "bin"), "--libdir": str(prefix / "lib"),
                "--incdir": str(prefix / "include"), "--cxx": "fixture-c++",
                "--cflags": "-I" + str(prefix / "include") + " -DROOT_VERSION=" + version,
                "--libs": "-L" + str(prefix / "lib") + " -lCore",
            }
            print(" ".join(answers[arg] for arg in sys.argv[1:]))
        """)
        self.write_tool(self.tools / "rootcling", """
            if "--version" in sys.argv:
                print(Path(os.environ["FAKE_ROOT_VERSION"]).read_text())
                sys.exit(0)
            record("dictionary")
            output = Path(sys.argv[sys.argv.index("-f") + 1])
            output.write_text(Path(os.environ["FAKE_ROOT_VERSION"]).read_text())
            output.with_name(output.stem + "_rdict.pcm").write_text("fixture pcm")
        """)
        self.write_tool(self.path_tools / "rootcling", """
            record("wrong-rootcling")
            sys.exit("rootcling must come from root-config --bindir")
        """)
        self.compiler = self.path_tools / "fixture-c++"
        self.write_tool(self.compiler, """
            if any(arg in ("--version", "-dumpversion", "-dumpfullversion", "-dumpmachine") for arg in sys.argv):
                print("fixture-compiler-1")
                sys.exit(0)
            output = Path(sys.argv[sys.argv.index("-o") + 1])
            compiling = "-c" in sys.argv
            if compiling:
                source = next(Path(arg) for arg in sys.argv[1:] if arg.endswith((".cxx", ".cc")))
                if source.name == "nanoDict.cc":
                    expected = Path(os.environ["FAKE_ROOT_VERSION"]).read_text()
                    if source.read_text() != expected:
                        sys.exit("stale dictionary reached compiler")
                record("compile", output=str(output))
                if "-MMD" in sys.argv:
                    dependency = Path(sys.argv[sys.argv.index("-MF") + 1]) if "-MF" in sys.argv else output.with_suffix(".d")
                    headers = [source.with_suffix(".h"), Path("Common.h")]
                    dependency.write_text(str(output) + ": " + str(source) + " " + " ".join(str(header) for header in headers if header.exists()) + "\\n")
            else:
                objects = [Path(arg) for arg in sys.argv[1:] if arg.endswith(".o")]
                if not all(path.exists() for path in objects):
                    sys.exit("link started without all objects")
                record("link", output=str(output))
            output.write_text("fixture output")
        """)
        # Delaying removal makes an incorrectly parallel clean/build observable.
        self.write_tool(self.path_tools / "rm", """
            cleaning = any(argument.endswith((".o", ".so")) for argument in sys.argv[1:])
            if cleaning:
                record("clean-start")
                time.sleep(0.1)
            for argument in sys.argv[1:]:
                if not argument.startswith("-"):
                    Path(argument).unlink(missing_ok=True)
            if cleaning:
                record("clean-end")
        """)

    def write_tool(self, path, body):
        prelude = """
import json
import os
from pathlib import Path
import sys
import time

def record(kind, **details):
    descriptor = os.open(os.environ["BUILD_EVENT_LOG"], os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
    try:
        os.write(descriptor, (json.dumps(dict(kind=kind, **details)) + "\\n").encode())
    finally:
        os.close(descriptor)
"""
        path.write_text("#!" + sys.executable + "\n" + prelude + textwrap.dedent(body))
        path.chmod(0o755)

    def make(self, *targets):
        result = subprocess.run(
            ["make", "-j3", "--no-print-directory", "ROOT_CONFIG=" + str(self.tools / "root-config"),
             "CXX=" + str(self.compiler), "HDRS=SplitCanvas.h Helper.h", *targets],
            cwd=self.directory, env=self.environment, text=True, capture_output=True, timeout=30,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def events(self):
        return [json.loads(line) for line in self.log.read_text().splitlines()] if self.log.exists() else []

    def reset_events(self):
        self.log.write_text("")

    def assert_complete_build(self):
        events = self.events()
        self.assertEqual([event["kind"] for event in events].count("dictionary"), 1, events)
        self.assertEqual({event["output"] for event in events if event["kind"] == "compile"},
                         {"SplitCanvas.o", "Helper.o", "nanoDict.o"}, events)
        self.assertEqual([event["kind"] for event in events].count("link"), 1, events)
        self.assertTrue((self.directory / "DrawingClasses.so").exists())
        self.assertTrue((self.directory / "nanoDict_rdict.pcm").exists())

    def test_stale_dictionary_rebuilt_with_selected_root_and_second_make_is_idle(self):
        (self.directory / "nanoDict.cc").write_text("obsolete ROOT dictionary")
        (self.directory / "nanoDict_rdict.pcm").write_text("obsolete pcm")
        self.make()
        self.assert_complete_build()
        self.assertNotIn("wrong-rootcling", [event["kind"] for event in self.events()])
        self.reset_events()
        self.make()
        self.assertEqual(self.events(), [])

    def test_root_change_at_same_config_path_rebuilds_dictionary_and_all_objects(self):
        self.make()
        self.reset_events()
        self.identity.write_text("6.40.fixture")
        self.make()
        self.assert_complete_build()

    def test_missing_dictionary_outputs_are_recreated(self):
        self.make()
        for output in ("nanoDict_rdict.pcm", "nanoDict.cc"):
            with self.subTest(output=output):
                (self.directory / output).unlink()
                self.reset_events()
                self.make()
                kinds = [event["kind"] for event in self.events()]
                self.assertEqual(kinds.count("dictionary"), 1, self.events())
                self.assertTrue((self.directory / "nanoDict.cc").exists())
                self.assertTrue((self.directory / "nanoDict_rdict.pcm").exists())

    def test_transitive_header_change_recompiles_dependent_objects(self):
        self.make()
        self.reset_events()
        (self.directory / "Common.h").write_text("// changed shared header\n")
        self.make()
        self.assertEqual({event["output"] for event in self.events() if event["kind"] == "compile"},
                         {"SplitCanvas.o", "Helper.o", "nanoDict.o"})
        self.assertEqual([event["kind"] for event in self.events()].count("link"), 1)

    def test_parallel_cl_finishes_cleanup_before_rebuilding(self):
        self.make()
        self.reset_events()
        self.make("cl")
        kinds = [event["kind"] for event in self.events()]
        self.assertIn("clean-end", kinds)
        self.assertLess(kinds.index("clean-end"), kinds.index("dictionary"))
        self.assertLess(kinds.index("clean-end"), kinds.index("compile"))
        self.assert_complete_build()


if __name__ == "__main__":
    unittest.main()
