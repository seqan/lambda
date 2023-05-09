cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

## Index input files

set(BASEURL "https://raw.githubusercontent.com/h-2/lambda-testdata/0a52e3bed1daf8646a496b63721b912fef3364bd")

declare_datasource (FILE db_nucl.fasta.gz
                    URL ${BASEURL}/input_files/db_nucl.fasta.gz
                    URL_HASH SHA256=614f8d7863c40facb7fffe666ce04341368a43ea90b8125cb675907f315bb0a2)

declare_datasource (FILE db_nucl_bs.fasta.gz
                    URL ${BASEURL}/input_files/db_nucl_bs.fasta.gz
                    URL_HASH SHA256=160375ac5ff4426f1768981a0215495efd6de0b847e3226c5c7c112370a599a1)

declare_datasource (FILE db_prot.fasta.gz
                    URL ${BASEURL}/input_files/db_prot.fasta.gz
                    URL_HASH SHA256=2c27f09f77e1f8ec0fea0aa1cae75f6634bab852a843a48033fc0f5251d57626)

## Index output files

declare_datasource (FILE db_nucl_fm.fasta.gz.lba
                    URL ${BASEURL}/output_files/db_nucl_fm.fasta.gz.lba
                    URL_HASH SHA256=68de7b2822b3f47268efee8becd400697cb1581a7055d1b581af365dd27f7a89)

declare_datasource (FILE db_nucl_bifm.fasta.gz.lba
                    URL ${BASEURL}/output_files/db_nucl_bifm.fasta.gz.lba
                    URL_HASH SHA256=e812f75f1672b5dc4b8cb6078ce209ec809d21783424836ccbe0c9b78f37887c)

declare_datasource (FILE db_nucl_bs_fm.fasta.gz.lba
                    URL ${BASEURL}/output_files/db_nucl_bs_fm.fasta.gz.lba
                    URL_HASH SHA256=b17d4e29f424af54fa802344ae9d4af974dbfc16d0445741e08de6f3cfface0b)

declare_datasource (FILE db_nucl_bs_bifm.fasta.gz.lba
                    URL ${BASEURL}/output_files/db_nucl_bs_bifm.fasta.gz.lba
                    URL_HASH SHA256=6dc216251505f0a913dffb9dc97699898492d511934dbcd1c9961337a151a62b)

declare_datasource (FILE db_prot_fm.fasta.gz.lba
                    URL ${BASEURL}/output_files/db_prot_fm.fasta.gz.lba
                    URL_HASH SHA256=57fa7a4319c5eb4dded850a79beefa1c5eca6927c91541cd7d658fde0ab7d166)

declare_datasource (FILE db_prot_bifm.fasta.gz.lba
                    URL ${BASEURL}/output_files/db_prot_bifm.fasta.gz.lba
                    URL_HASH SHA256=76d0405235aff8145fe0d90399ffb7d3619ad2fedc37776771a1ca8a38b79a5d)

declare_datasource (FILE db_trans_fm.fasta.gz.lba
                    URL ${BASEURL}/output_files/db_trans_fm.fasta.gz.lba
                    URL_HASH SHA256=dc94a97d19a2ad224b995640bcacfc8fb6c7aa87ef134f3e61647c5b31286696)

declare_datasource (FILE db_trans_bifm.fasta.gz.lba
                    URL ${BASEURL}/output_files/db_trans_bifm.fasta.gz.lba
                    URL_HASH SHA256=3ffce0f1115888bce562460e7f4f1222af7052cd8a85a32db0813eac317d8be5)

## Query input files

declare_datasource (FILE queries_nucl.fasta.gz
                    URL ${BASEURL}/input_files/queries_nucl.fasta.gz
                    URL_HASH SHA256=7a8adcc7ee5d967992a0624b9fb704118223abd112534c2c6456104b30ddad88)

declare_datasource (FILE queries_nucl_bs.fasta.gz
                    URL ${BASEURL}/input_files/queries_nucl_bs.fasta.gz
                    URL_HASH SHA256=a358ace3e6f35bd379854fd8870a63a6da62b6bd3e72d95f11e5386398881f0a)

declare_datasource (FILE queries_prot.fasta.gz
                    URL ${BASEURL}/input_files/queries_prot.fasta.gz
                    URL_HASH SHA256=e21411f422c1dd844696c8ca8a08b93e8243d7f77307731c5f2ee8dcc60d0703)

## Search output files

declare_datasource (FILE output_blastn_fm.bam
                    URL ${BASEURL}/output_files/output_blastn_fm.bam
                    URL_HASH SHA256=89e4052abc650d8ab01e927a7e7783de78f8d7e8029f0c406c226e8802176260)

declare_datasource (FILE output_blastn_fm.m0
                    URL ${BASEURL}/output_files/output_blastn_fm.m0
                    URL_HASH SHA256=cac63b12c30ffef52393224ffbe01a17c89ba8b89cb8349acbbac869eb51273b)

declare_datasource (FILE output_blastn_fm.m8
                    URL ${BASEURL}/output_files/output_blastn_fm.m8
                    URL_HASH SHA256=2c139f4cf35d3ba5afdf3c0ff7f04c77904d73b7a66895c5d86fcc987b3cf64a)

declare_datasource (FILE output_blastn_fm.m9
                    URL ${BASEURL}/output_files/output_blastn_fm.m9
                    URL_HASH SHA256=5f7a343a5d19238332af16d08b38ea70aa7da5285a8874b361cffe763505a536)

declare_datasource (FILE output_blastn_fm.sam
                    URL ${BASEURL}/output_files/output_blastn_fm.sam
                    URL_HASH SHA256=61b1499cead7990c102f7d66ea7323baa1756ce8c832c3caf68d3134df3e4841)

declare_datasource (FILE output_blastn_bs_fm.bam
                    URL ${BASEURL}/output_files/output_blastn_bs_fm.bam
                    URL_HASH SHA256=7ffbf515acb743e2f20fe2210493b8de8af88c17ad5339867f1afb1b1ec11bb0)

declare_datasource (FILE output_blastn_bs_fm.m0
                    URL ${BASEURL}/output_files/output_blastn_bs_fm.m0
                    URL_HASH SHA256=8234dc029edac95d878ecea6dac7239584d8c5558f974912a13cc0e37988bc56)

declare_datasource (FILE output_blastn_bs_fm.m8
                    URL ${BASEURL}/output_files/output_blastn_bs_fm.m8
                    URL_HASH SHA256=145cf9f87566e820e05c9e84f94037860dbc0c809ae5148fda8b626ce172d37b)

declare_datasource (FILE output_blastn_bs_fm.m9
                    URL ${BASEURL}/output_files/output_blastn_bs_fm.m9
                    URL_HASH SHA256=329298632838e82c7fa5c0cae908d0c3c4fa429166d8a1224f0fabef164d5834)

declare_datasource (FILE output_blastn_bs_fm.sam
                    URL ${BASEURL}/output_files/output_blastn_bs_fm.sam
                    URL_HASH SHA256=5e975b66c5d034afd7a0adc438eb2951367377b2bf5e55723d0fa00ab2bc5c20)

declare_datasource (FILE output_blastp_fm.bam
                    URL ${BASEURL}/output_files/output_blastp_fm.bam
                    URL_HASH SHA256=a750d396a54ccd7b414f2567d4e96596e9cbd576bb6e2cb40ae7cdf0d049d546)

declare_datasource (FILE output_blastp_fm.m0
                    URL ${BASEURL}/output_files/output_blastp_fm.m0
                    URL_HASH SHA256=27c2da9595f3276255eaaaf3910c32809e9472dd4bb190fcea95ebfd14249eb7)

declare_datasource (FILE output_blastp_fm.m8
                    URL ${BASEURL}/output_files/output_blastp_fm.m8
                    URL_HASH SHA256=a3e2c3c0d305cd8ba2e8a5fe52bb1656e41d3091c24c2b80e7a5a5bb678f401b)

declare_datasource (FILE output_blastp_fm.m9
                    URL ${BASEURL}/output_files/output_blastp_fm.m9
                    URL_HASH SHA256=7ee3cc8cf7fb2ea784ec404ec97b9cede14875cfba8de597a16e7b12dc19a1f4)

declare_datasource (FILE output_blastp_fm.sam
                    URL ${BASEURL}/output_files/output_blastp_fm.sam
                    URL_HASH SHA256=57e2c5c7297a259441976afbc3f4e77c6ea8f92ec8b02e08a103937da8cd32a0)

declare_datasource (FILE output_blastx_fm.bam
                    URL ${BASEURL}/output_files/output_blastx_fm.bam
                    URL_HASH SHA256=db5daa713193dfcf693cc24e52e4e64631ec58dcaba7eea67901b916d4951843)

declare_datasource (FILE output_blastx_fm.m0
                    URL ${BASEURL}/output_files/output_blastx_fm.m0
                    URL_HASH SHA256=0db7f6f44ba0a4cf5052b1c3abf094e374c881ba79e3c9440156cc3f5e7f7c52)

declare_datasource (FILE output_blastx_fm.m8
                    URL ${BASEURL}/output_files/output_blastx_fm.m8
                    URL_HASH SHA256=bfa3b1d8c0aebea0f0cd2c99fedfe582020c0ac24ab5cd51ecd67649ec75f467)

declare_datasource (FILE output_blastx_fm.m9
                    URL ${BASEURL}/output_files/output_blastx_fm.m9
                    URL_HASH SHA256=cd60e7bf440cfc70362d1757842d9cceb5361cfa2e8cbf9223650617dc54dc47)

declare_datasource (FILE output_blastx_fm.sam
                    URL ${BASEURL}/output_files/output_blastx_fm.sam
                    URL_HASH SHA256=47a7dcc2645221645368acc23d891ca10d04e3f797e0962a7d2734a8415d5d01)

declare_datasource (FILE output_tblastn_fm.bam
                    URL ${BASEURL}/output_files/output_tblastn_fm.bam
                    URL_HASH SHA256=d5dc8d696932d4bc595388368afa60fc6281f5a74be6896e3b6301c048d602af)

declare_datasource (FILE output_tblastn_fm.m0
                    URL ${BASEURL}/output_files/output_tblastn_fm.m0
                    URL_HASH SHA256=dbe9949ba4c70f30adb9bf22362a09934c735df31b2d81bbb674373ac897ba90)

declare_datasource (FILE output_tblastn_fm.m8
                    URL ${BASEURL}/output_files/output_tblastn_fm.m8
                    URL_HASH SHA256=ec174aad5c9199673cd6b9123bf5e9eada348ca3bc68b3409e36b81cbf533197)

declare_datasource (FILE output_tblastn_fm.m9
                    URL ${BASEURL}/output_files/output_tblastn_fm.m9
                    URL_HASH SHA256=4498af84a412e273067c3396845b3c8d6a1ef9c72f1c3368207ef6d9ed89b517)

declare_datasource (FILE output_tblastn_fm.sam
                    URL ${BASEURL}/output_files/output_tblastn_fm.sam
                    URL_HASH SHA256=f674bcb028742f89b367ee2c1b59284a83ecbac0c29c793574f744887162552d)

declare_datasource (FILE output_tblastx_fm.bam
                    URL ${BASEURL}/output_files/output_tblastx_fm.bam
                    URL_HASH SHA256=eded153112ea8975a6fe736e74c931bfd212c06fa61c0197dafc7911e2058d82)

declare_datasource (FILE output_tblastx_fm.m0
                    URL ${BASEURL}/output_files/output_tblastx_fm.m0
                    URL_HASH SHA256=52b19c9a3545471c52a826e95610e67dcf0e6b7851f852a09a54f97b6eaa677b)

declare_datasource (FILE output_tblastx_fm.m8
                    URL ${BASEURL}/output_files/output_tblastx_fm.m8
                    URL_HASH SHA256=a2215049149c581dd34c3c4d28fbb90d654cd59b399b2c7e001af08d858b82c6)

declare_datasource (FILE output_tblastx_fm.m9
                    URL ${BASEURL}/output_files/output_tblastx_fm.m9
                    URL_HASH SHA256=385a42654fa9ea76ee6dc851702bfad487aa75ff39fab424ec0eec84a2069238)

declare_datasource (FILE output_tblastx_fm.sam
                    URL ${BASEURL}/output_files/output_tblastx_fm.sam
                    URL_HASH SHA256=3b2046496fc651fa1987b666b899c6b6d2fbc62bc6af623f575a7c7e2ef3b32a)

declare_datasource (FILE output_blastn_fm_fast.sam
                    URL ${BASEURL}/output_files/output_blastn_fm_fast.sam
                    URL_HASH SHA256=de57c8ee12447df2d635b1b2e9943bca9d4089e3cba433bb9e78ce4ec7db729f)

declare_datasource (FILE output_blastn_fm_fast.m8
                    URL ${BASEURL}/output_files/output_blastn_fm_fast.m8
                    URL_HASH SHA256=1a7d1a07b8673e7fb6ee0033d27aa705be5e366febb95e8d465947f0680f3521)

declare_datasource (FILE output_blastn_bs_fm_fast.sam
                    URL ${BASEURL}/output_files/output_blastn_bs_fm_fast.sam
                    URL_HASH SHA256=5e975b66c5d034afd7a0adc438eb2951367377b2bf5e55723d0fa00ab2bc5c20)

declare_datasource (FILE output_blastn_bs_fm_fast.m8
                    URL ${BASEURL}/output_files/output_blastn_bs_fm_fast.m8
                    URL_HASH SHA256=145cf9f87566e820e05c9e84f94037860dbc0c809ae5148fda8b626ce172d37b)

declare_datasource (FILE output_blastp_fm_fast.sam
                    URL ${BASEURL}/output_files/output_blastp_fm_fast.sam
                    URL_HASH SHA256=410ff4b7316dd8353cfdc42b68eb01bbdda189b6bbfbbe3fc8333eb7d80328ba)

declare_datasource (FILE output_blastp_fm_fast.m8
                    URL ${BASEURL}/output_files/output_blastp_fm_fast.m8
                    URL_HASH SHA256=47fce11a96fbcce92f62bfacb1545de432e4d97f1d5b74c4d1e6f3fd9418898e)

declare_datasource (FILE output_blastx_fm_fast.sam
                    URL ${BASEURL}/output_files/output_blastx_fm_fast.sam
                    URL_HASH SHA256=512373692eb80bf426f027dac5a4355c6f94589cc60bd97d16e9964bc5d8689c)

declare_datasource (FILE output_blastx_fm_fast.m8
                    URL ${BASEURL}/output_files/output_blastx_fm_fast.m8
                    URL_HASH SHA256=b039b8fabbcfb67bdd2d92c99c4db161f820bc9ea1886ac088327c0b3803f3e3)

declare_datasource (FILE output_tblastn_fm_fast.sam
                    URL ${BASEURL}/output_files/output_tblastn_fm_fast.sam
                    URL_HASH SHA256=d74c48bb71ed0cff236d3f9b14ae694155d1e798f309a84113a3748fba24605d)

declare_datasource (FILE output_tblastn_fm_fast.m8
                    URL ${BASEURL}/output_files/output_tblastn_fm_fast.m8
                    URL_HASH SHA256=ae3f8e5f2966579c09688667d56f98b48a5ac5cd513d408cd2de09f7090aa1d9)

declare_datasource (FILE output_tblastx_fm_fast.sam
                    URL ${BASEURL}/output_files/output_tblastx_fm_fast.sam
                    URL_HASH SHA256=ffae5e644c963294766a41879fb990c11ec308e56df123649ea2ce9d677eaa93)

declare_datasource (FILE output_tblastx_fm_fast.m8
                    URL ${BASEURL}/output_files/output_tblastx_fm_fast.m8
                    URL_HASH SHA256=686450a6a1ff546244a75bfb05925963bbce79ecf9a4c67a3566d968134cd0e7)

declare_datasource (FILE output_blastn_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_blastn_fm_sensitive.sam
                    URL_HASH SHA256=34e6942316c33300004f2d28cc8644dd4bd6cc73f18e70d3a0a357f74519e2ea)

declare_datasource (FILE output_blastn_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_blastn_fm_sensitive.m8
                    URL_HASH SHA256=16b36d8cab68aa08ca6130be8b91b4c9f48411eaf6edcf4d97b87bde5fe76134)

declare_datasource (FILE output_blastn_bs_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_blastn_bs_fm_sensitive.sam
                    URL_HASH SHA256=9369c82e236d879928eb5b6731fea04042193a4df7d92fd7fd95cb82ed368d8b)

declare_datasource (FILE output_blastn_bs_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_blastn_bs_fm_sensitive.m8
                    URL_HASH SHA256=09642471eb94906c1706deb597b58b9bdea5c1924622f5668cc7d9d55c88ba27)

declare_datasource (FILE output_blastp_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_blastp_fm_sensitive.sam
                    URL_HASH SHA256=582fe1d78cc03eb94e70b78ab8275d7c34a597bb7a896277b5a93573c489dfa8)

declare_datasource (FILE output_blastp_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_blastp_fm_sensitive.m8
                    URL_HASH SHA256=0fe2075c016e847b9ffb2b817e0c59326228d2fbeaef86af5542ad187f4ba283)

declare_datasource (FILE output_blastx_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_blastx_fm_sensitive.sam
                    URL_HASH SHA256=0009d2e5cb5019d8c9867b898f9cd2593306069c6b0bbad33b97914cef0d741b)

declare_datasource (FILE output_blastx_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_blastx_fm_sensitive.m8
                    URL_HASH SHA256=93d341c0cf30a9f42221008fad3ff46c1d213f3c89ce1168871d259b02824260)

declare_datasource (FILE output_tblastn_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_tblastn_fm_sensitive.sam
                    URL_HASH SHA256=f1eb59c90b051798ef0c18f4000c0354fe8c97f83c384f614a19d1d294f7aab6)

declare_datasource (FILE output_tblastn_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_tblastn_fm_sensitive.m8
                    URL_HASH SHA256=292ec38ced08c68eca49dd0cd08bc0d89d9c011cddb7ccba445c9c61708cc51b)

declare_datasource (FILE output_tblastx_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_tblastx_fm_sensitive.sam
                    URL_HASH SHA256=727edd5f9c1794ba61af90ccbd3008f1c736383031591a3c76537a55333fcfa3)

declare_datasource (FILE output_tblastx_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_tblastx_fm_sensitive.m8
                    URL_HASH SHA256=fb0109e002e6ef9805a3aa270ec0455e3c67ceed34a0a297fb6d463472086005)
