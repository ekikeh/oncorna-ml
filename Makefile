.PHONY: test lint check

test:
	python -m pytest

lint:
	ruff check .

check: lint test
