@echo off

cd /d %~dp0

GenotypeHarmonizer.bat --inputType PLINK_BED --input .\exampleData\hapmap3CeuChr20B37Mb6RandomStrand --update-id --outputType PLINK_BED --output .\exampleOutput\binaryPlinkExampleOut --refType VCF --ref .\exampleData\1000gCeuChr20Mb6

