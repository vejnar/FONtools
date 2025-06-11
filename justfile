default:
    @just --list

_check-uv:
    @if [ ! -f "$(which uv || true)" ]; then echo "uv not found"; exit 1 ; fi

clean:
    @find . -type f -name "*.py[co]" -delete -o -type d -name __pycache__ -delete
    @rm -rf src/*.egg-info/ build/ dist/ .pytest_cache/ .ruff_cache/

publish: _check-uv clean
    @if [ -z "$(git describe --exact-match 2>/dev/null || true)" ]; then echo "Current commit is not a tag; abort task"; exit 1; fi
    uv build
    uv publish
