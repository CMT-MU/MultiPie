build:
	jupyter book build --all -v .

upload:	build
	ghp-import -n -p -f _build/html

clean:
	rm -rf _build
