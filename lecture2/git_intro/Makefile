all: reveal beamer pdf

reveal:
	pandoc --self-contained --variable theme="simple" -t revealjs -s github.md -o github.html

beamer:
	pandoc -t beamer -s github.md -o github_slides.pdf

pdf:
	pandoc -t beamer -s github.md -o github.pdf
