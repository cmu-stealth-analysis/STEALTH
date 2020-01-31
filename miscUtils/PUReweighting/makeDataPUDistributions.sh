pileupCalc.py -i Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt --inputLumiJSON pileup_info_2016.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 dataPU_2016.root
mv -v dataPU_2016.root ../../getMCSystematics/data/dataPU_2016.root

pileupCalc.py -i Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt --inputLumiJSON pileup_info_2017.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 dataPU_2017.root
mv -v dataPU_2017.root ../../getMCSystematics/data/dataPU_2017.root

pileupCalc.py -i Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt --inputLumiJSON pileup_info_2018.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 dataPU_2018.root
mv -v dataPU_2018.root ../../getMCSystematics/data/dataPU_2018.root
