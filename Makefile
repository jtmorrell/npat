install :
	python setup.py install --user

test :
	echo "WARNING: Make test has not been implemented yet"

validate :
	echo "WARNING: Make validate has not been implemented yet"

all :
	make install

.PHONY : install