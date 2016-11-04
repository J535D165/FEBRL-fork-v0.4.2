# Shell script to run all Febrl testing scripts
#
# Peter Christen, 24/10/2007
# -----------------------------------------------------------------------------
# Start with basic modules before doing more complex modules

echo "Run auxiliaryTest.py"
python auxiliaryTest.py

echo "Run mymathTest.py"
python mymathTest.py

echo "Run encodeTest.py"
python encodeTest.py

echo "Run stringcmpTest.py"
python stringcmpTest.py

echo "Run phonenumTest.py"
python phonenumTest.py

echo "Run datasetTest.py"
python datasetTest.py

echo "Run lookupTest.py"
python lookupTest.py

echo "Run simplehmmTest.py"
python simplehmmTest.py

echo "Run indexingTest.py"
python indexingTest.py

echo "Run comparisonTest.py"
python comparisonTest.py

echo "Run classificationTest.py"
python classificationTest.py

echo "Run standardisationTest.py"
python standardisationTest.py


