import inspect
import unittest
from os import path

import astair2.cli as run

current = path.abspath(path.dirname(__file__))


class AsTaiRunTest(unittest.TestCase):
    """"""

    def test_run(self):
        """"""
        self.assertEqual(inspect.cleandoc(run.cli.help), 'asTair2 (tools for processing cytosine modification sequencing data)')
        
if __name__ == '__main__':
    unittest.main()
