.PHONY: build install test coverage

build:
	python setup.py bdist_wheel
	python setup.py sdist

upload_test:
	twine upload -r testpypi dist/*

upload:
	twine upload -r pypi dist/*

test:
	pytest -svv

install_local:
	pip install dist/novana*.whl