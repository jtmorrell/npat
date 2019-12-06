import pytest
### Run using command `pytest` in npat/test directory


def setup_function():
	print('Setup')


def test_add():
	assert 1+1==2

def test_mult():
	assert 1*1==2


def teardown_function():
	print('Teardown')