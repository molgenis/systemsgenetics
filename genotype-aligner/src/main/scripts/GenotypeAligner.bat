@echo off

set genotypeAlignerDir=%~dp0

java -Xmx1g -jar "%genotypeAlignerDir%GenotypeAligner.jar" %*