all: clean build test install

clean:
	@echo "Cleaning up"

build:
	perl Build.PL && ./Build installdeps && ./Build

test:
	./Build test

install: build
	./Build install
