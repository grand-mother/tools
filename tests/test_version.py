# -*- coding: utf-8 -*-
"""
Unit tests for the tools.version module
"""

import unittest
import sys

import tools
from grand_pkg import git


try:
    import tools.version
except:
    # Skip version tests for non release builds
    pass
else:
    class VersionTest(unittest.TestCase):
        """Unit tests for the version module"""

        def test_hash(self):
            githash = git("rev-parse", "HEAD")
            self.assertEqual(githash.strip(), tools.version.__githash__)

        def test_version(self):
            self.assertIsNotNone(tools.version.__version__)


if __name__ == "__main__":
    unittest.main()
