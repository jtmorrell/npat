install :
	python setup.py install --user

test :
	cd test && python test.py

validate :
	echo "WARNING: Make validate has not been implemented yet"

all :
	make install

.PHONY : test