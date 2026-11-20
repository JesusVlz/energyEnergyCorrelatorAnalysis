"""End-to-end verification of MixingLookupTable, donor file-open caps, sharding, and CLI tools."""

from collections import Counter
import os
import pathlib
import shlex
import subprocess
import tempfile
import unittest


OPEN_TRACER_C = r"""
#define _GNU_SOURCE
#include <dlfcn.h>
#include <fcntl.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static void log_root_open(const char* path) {
  if (!path) return;
  const size_t len = strlen(path);
  if (len < 5 || strcmp(path + len - 5, ".root") != 0) return;
  const char* log_path = getenv("EEC_ROOT_OPEN_LOG");
  if (!log_path || !log_path[0]) return;
  int (*real_open)(const char*, int, ...) = dlsym(RTLD_NEXT, "open");
  if (!real_open) return;
  int fd = real_open(log_path, O_WRONLY | O_CREAT | O_APPEND, 0644);
  if (fd >= 0) {
    ssize_t unused1 = write(fd, path, len);
    ssize_t unused2 = write(fd, "\n", 1);
    (void)unused1;
    (void)unused2;
    close(fd);
  }
}

int open(const char* pathname, int flags, ...) {
  mode_t mode = 0;
  if (flags & (O_CREAT | O_TMPFILE)) {
    va_list ap;
    va_start(ap, flags);
    mode = va_arg(ap, int);
    va_end(ap);
  }
  log_root_open(pathname);
  int (*real_open)(const char*, int, ...) = dlsym(RTLD_NEXT, "open");
  return real_open(pathname, flags, mode);
}

int open64(const char* pathname, int flags, ...) {
  mode_t mode = 0;
  if (flags & (O_CREAT | O_TMPFILE)) {
    va_list ap;
    va_start(ap, flags);
    mode = va_arg(ap, int);
    va_end(ap);
  }
  log_root_open(pathname);
  int (*real_open64)(const char*, int, ...) = dlsym(RTLD_NEXT, "open64");
  return real_open64(pathname, flags, mode);
}

int openat(int dirfd, const char* pathname, int flags, ...) {
  mode_t mode = 0;
  if (flags & (O_CREAT | O_TMPFILE)) {
    va_list ap;
    va_start(ap, flags);
    mode = va_arg(ap, int);
    va_end(ap);
  }
  log_root_open(pathname);
  int (*real_openat)(int, const char*, int, ...) = dlsym(RTLD_NEXT, "openat");
  return real_openat(dirfd, pathname, flags, mode);
}

int openat64(int dirfd, const char* pathname, int flags, ...) {
  mode_t mode = 0;
  if (flags & (O_CREAT | O_TMPFILE)) {
    va_list ap;
    va_start(ap, flags);
    mode = va_arg(ap, int);
    va_end(ap);
  }
  log_root_open(pathname);
  int (*real_openat64)(int, const char*, int, ...) = dlsym(RTLD_NEXT, "openat64");
  return real_openat64(dirfd, pathname, flags, mode);
}
"""


def write_oo_mixing_card(base_card: pathlib.Path, out_card: pathlib.Path, overrides: dict[str, str]) -> None:
    lines = base_card.read_text().splitlines()
    filtered = []
    override_keys = set(overrides.keys())
    for line in lines:
        stripped = line.strip()
        if stripped and not stripped.startswith('#'):
            key = stripped.split()[0]
            if key in override_keys:
                continue
        filtered.append(line)
    for key, val in overrides.items():
        filtered.append(f"{key} {val}")
    out_card.write_text("\n".join(filtered) + "\n")


class MixingLookupTableTest(unittest.TestCase):
    def test_lookup_table_open_caps_sharding_and_equivalence(self):
        repo = pathlib.Path(__file__).resolve().parents[2]
        subprocess.run(['make', '-j4'], cwd=repo, check=True)

        flags = shlex.split(subprocess.check_output(['root-config', '--cflags', '--libs'], text=True))
        objects = [str(p.relative_to(repo)) for p in sorted((repo / 'src').glob('*.o'))]

        with tempfile.TemporaryDirectory(prefix='eec_mixing_lookup_test_') as tmp:
            tmpdir = pathlib.Path(tmp)
            fixtures = tmpdir / 'fixtures'

            # 1. Compile and run C++ equivalence + leak-check suite
            cpp_bin = tmpdir / 'test_mixing_lookup_cpp'
            subprocess.run(
                ['g++', '-g', '-Wall', '-fno-access-control', '-Isrc', '-I.',
                 str(pathlib.Path(__file__).with_suffix('.cxx')),
                 *objects, *flags, '-lEG', '-o', str(cpp_bin)],
                cwd=repo, check=True
            )
            subprocess.run([str(cpp_bin), '--all', str(fixtures)], cwd=repo, check=True)

            # 2. Compile LD_PRELOAD .root open tracer to count exact kernel-level ROOT file opens
            tracer_c = tmpdir / 'open_tracer.c'
            tracer_so = tmpdir / 'open_tracer.so'
            tracer_c.write_text(OPEN_TRACER_C)
            subprocess.run(['gcc', '-O2', '-shared', '-fPIC', str(tracer_c), '-ldl', '-o', str(tracer_so)], check=True)

            def run_traced(cmd: list[str], log_name: str, check: bool = True):
                log_file = tmpdir / log_name
                if log_file.exists():
                    log_file.unlink()
                env = dict(os.environ)
                env['LD_PRELOAD'] = str(tracer_so)
                env['EEC_ROOT_OPEN_LOG'] = str(log_file)
                res = subprocess.run(cmd, cwd=repo, env=env, capture_output=True, text=True)
                if check and res.returncode != 0:
                    self.fail(f"Command {cmd} failed ({res.returncode}):\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
                opens = log_file.read_text().splitlines() if log_file.exists() else []
                donor_opens = [os.path.basename(p) for p in opens if os.path.basename(p).startswith('donor_')]
                return res, Counter(donor_opens)

            donor_list = fixtures / 'oo_donors.txt'
            sidecar_lookup = pathlib.Path(str(donor_list) + '.lookup.tsv')
            if sidecar_lookup.exists():
                sidecar_lookup.unlink()

            base_card = repo / 'cardEECOOInclusiveR04.input'
            card_default = tmpdir / 'card_oo_mix.input'
            write_oo_mixing_card(base_card, card_default, {
                'BackgroundMethods': '6',
                'EventMatchForMixing': '0',
                'BinnedMixing': '0',
                'MixingFileList': str(donor_list),
                'MixingMaxOpenFiles': '3',
                'ApplyJetVeto': '0',
                'Run3JetID': '0',
                'Run3JetQA': '0',
            })

            # 3. Verify buildMixingLookupTable CLI opens each of the 10 donor files AT MOST ONCE
            _, build_counts = run_traced(
                ['./buildMixingLookupTable', str(donor_list), str(card_default)],
                'build_lookup.log'
            )
            self.assertTrue(sidecar_lookup.is_file())
            self.assertEqual(len(build_counts), 10, f"Expected 10 distinct donor files scanned, got {build_counts}")
            for fname, count in build_counts.items():
                self.assertEqual(count, 1, f"{fname} was opened {count} times during lookup build (expected 1)")

            # 4. Zero-jet signal file + prebuilt sidecar lookup table -> ZERO donor ROOT files opened!
            _, zero_jet_counts = run_traced(
                ['./eecAnalysis', str(fixtures / 'signal_zero_jets.root'), str(card_default),
                 str(tmpdir / 'out_zero.root'), '0', 'true', '1'],
                'zero_jets.log'
            )
            self.assertEqual(sum(zero_jet_counts.values()), 0,
                             f"Expected 0 donor ROOT file opens for 0-jet signal file, got {zero_jet_counts}")

            # 5. Signal file with qualifying jets + MixingMaxOpenFiles=3 across 4 simulated batch/CRAB shards:
            #    - Each job opens <= 3 distinct donor files, each opened AT MOST ONCE!
            #    - Different shard IDs spread load across different subsets of the 10-file donor pool!
            all_shards_seen = set()
            shard_sets = []
            for shard_id in (0, 3, 6, 9):
                _, job_counts = run_traced(
                    ['./eecAnalysis', str(fixtures / 'signal_jets.root'), str(card_default),
                     str(tmpdir / f'out_shard_{shard_id}.root'), '0', 'true', str(shard_id)],
                    f'shard_{shard_id}.log'
                )
                self.assertGreaterEqual(len(job_counts), 1, "Expected at least 1 donor file opened for matched jets")
                self.assertLessEqual(len(job_counts), 3,
                                     f"MixingMaxOpenFiles=3 violated: opened {len(job_counts)} donor files ({job_counts})")
                for fname, count in job_counts.items():
                    self.assertEqual(count, 1, f"Donor file {fname} reopened {count} times in a single job!")
                shard_sets.append(frozenset(job_counts.keys()))
                all_shards_seen.update(job_counts.keys())

            self.assertGreater(len(set(shard_sets)), 1, f"All shards opened identical donor files: {shard_sets}")
            self.assertGreater(len(all_shards_seen), 3, f"Sharding failed to spread across pool: {all_shards_seen}")

            # 6. Fallback mode (no sidecar .lookup.tsv) with MixingMaxOpenFiles=3:
            #    Only the 3 active shard files may be touched (never all 10 files in the list!)
            sidecar_backup = tmpdir / 'backup.lookup.tsv'
            sidecar_lookup.rename(sidecar_backup)
            _, fallback_counts = run_traced(
                ['./eecAnalysis', str(fixtures / 'signal_jets.root'), str(card_default),
                 str(tmpdir / 'out_fallback.root'), '0', 'true', '2'],
                'fallback.log'
            )
            self.assertLessEqual(len(fallback_counts), 3,
                                 f"Fallback mode touched more than MixingMaxOpenFiles=3 files: {fallback_counts}")
            for fname, count in fallback_counts.items():
                # At most 1 open for metadata scan + 1 open for trackTree prefetch = <= 2 opens
                self.assertLessEqual(count, 2, f"Fallback mode opened {fname} {count} times (expected <= 2)")

            # Restore sidecar lookup and test mismatch handling:
            # (a) Auto-discovered sidecar with mismatched ZVertexCut logs a note and falls back cleanly
            sidecar_backup.rename(sidecar_lookup)
            card_mismatch_sidecar = tmpdir / 'card_mismatch_sidecar.input'
            write_oo_mixing_card(card_default, card_mismatch_sidecar, {'ZVertexCut': '10.0'})
            res_sidecar, _ = run_traced(
                ['./eecAnalysis', str(fixtures / 'signal_zero_jets.root'), str(card_mismatch_sidecar),
                 str(tmpdir / 'out_mismatch_sidecar.root'), '0', 'true', '1'],
                'mismatch_sidecar.log'
            )
            self.assertIn('skipped (ZVertexCut mismatch)', res_sidecar.stdout)

            # (b) Explicit MixingLookupTable pointing to mismatched table MUST fail fast!
            card_mismatch_explicit = tmpdir / 'card_mismatch_explicit.input'
            write_oo_mixing_card(card_mismatch_sidecar, card_mismatch_explicit, {
                'MixingLookupTable': str(sidecar_lookup),
            })
            res_explicit, _ = run_traced(
                ['./eecAnalysis', str(fixtures / 'signal_zero_jets.root'), str(card_mismatch_explicit),
                 str(tmpdir / 'out_mismatch_explicit.root'), '0', 'true', '1'],
                'mismatch_explicit.log',
                check=False
            )
            self.assertNotEqual(res_explicit.returncode, 0,
                                "Explicitly configured mismatched MixingLookupTable must fail fast")


if __name__ == '__main__':
    unittest.main()
