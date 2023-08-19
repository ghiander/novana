.PHONY: build install test coverage

build:
	python setup.py bdist_wheel
	python setup.py sdist

upload_test:
	twine upload -r testpypi dist/*

test:
	pytest -svv