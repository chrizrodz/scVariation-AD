git clone https://github.com/lanagarmire/deepimpute
cd deepimpute/
python3 setup.py
cd deepimpute/
cat ../../GrubmanData/scRNA_rawCounts.tsv | tr "\\t" "," > ../../GrubmanData/scRNA_rawCounts.csv
python3 deepimpute/deepImpute.py ../GrubmanData/scRNA_rawCounts.tsv 
python3 deepImpute.py ../../GrubmanData/scRNA_rawCounts.csv -o ../../GrubmanData/test_impute.csv --cell-axis columns
cd deepimpute/
