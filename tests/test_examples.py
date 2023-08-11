"""
Tests the examples.
Notes: The paths are relative to '../'
"""
import logging
import os
from testbook import testbook

print(">>>>>>>>>>>>>>>>> current dir:"+os.getcwd())
if not os.getcwd().endswith("/examples"):
    os.chdir('./examples/')
    print(">>>>>>>>>>>>>>>>> current dir:"+os.getcwd())

@testbook("viz_image.ipynb", execute=True)
def test_viz_image(tb):
    pass # execute only because tests are present in the notebook itself
    return
