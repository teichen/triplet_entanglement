.PHONY: pip-install-dev clean test-dev

help:
	@echo "pip-install-dev - Install the Python dependencies for the software"
	@echo "clean - clean out the code"
	@echo "test-dev - master testing for the software"

pip-install-dev:
	pip install -r requirements.txt

clean:
	rm -fr build/

test-dev:
	pytest unit_tests/
