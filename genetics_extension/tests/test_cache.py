import tempfile
import unittest
from pathlib import Path

from cache import run_stage


class CheckpointContracts(unittest.TestCase):
    def test_reuse_does_not_call_producer(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            code = root / 'fixture.py'
            code.write_text('# fixture\n')
            calls = []

            def produce(output):
                calls.append(True)
                (output / 'result.txt').write_text('validated\n')
                return {'fixture': True}

            first, manifest = run_stage(root, 'fixture', {'input': 'a'}, {}, [code], produce)
            second, reused = run_stage(root, 'fixture', {'input': 'a'}, {}, [code], produce)
            self.assertEqual(first, second)
            self.assertEqual(manifest, reused)
            self.assertEqual(len(calls), 1)

    def test_tampered_output_is_preserved_and_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            code = root / 'fixture.py'
            code.write_text('# fixture\n')

            def produce(output):
                (output / 'result.txt').write_text('original\n')
                return {}

            output, _ = run_stage(root, 'fixture', {'input': 'a'}, {}, [code], produce)
            (output / 'result.txt').write_text('changed\n')
            with self.assertRaisesRegex(RuntimeError, 'checksum mismatch'):
                run_stage(root, 'fixture', {'input': 'a'}, {}, [code], produce)
            self.assertEqual((output / 'result.txt').read_text(), 'changed\n')

    def test_interruption_is_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            code = root / 'fixture.py'
            code.write_text('# fixture\n')

            def fail(output):
                (output / 'partial.txt').write_text('incomplete\n')
                raise RuntimeError('fixture interruption')

            with self.assertRaisesRegex(RuntimeError, 'fixture interruption'):
                run_stage(root, 'fixture', {'input': 'a'}, {}, [code], fail)
            directories = list((root / 'results').iterdir())
            self.assertEqual(len(directories), 1)
            self.assertTrue((directories[0] / 'FAILED.json').is_file())
            with self.assertRaisesRegex(RuntimeError, 'incomplete checkpoint'):
                run_stage(root, 'fixture', {'input': 'a'}, {}, [code], fail)
            self.assertTrue((directories[0] / 'partial.txt').is_file())


if __name__ == '__main__':
    unittest.main()
