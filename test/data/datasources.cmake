cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

## Index input files

declare_datasource (FILE db_nucl.fasta.gz
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/input_files/db_nucl.fasta.gz
                    URL_HASH SHA256=614f8d7863c40facb7fffe666ce04341368a43ea90b8125cb675907f315bb0a2)

declare_datasource (FILE db_nucl_bs.fasta.gz
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/input_files/db_nucl_bs.fasta.gz
                    URL_HASH SHA256=160375ac5ff4426f1768981a0215495efd6de0b847e3226c5c7c112370a599a1)

declare_datasource (FILE db_prot.fasta.gz
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/input_files/db_prot.fasta.gz
                    URL_HASH SHA256=2c27f09f77e1f8ec0fea0aa1cae75f6634bab852a843a48033fc0f5251d57626)

## Index output files

declare_datasource (FILE db_nucl_fm.fasta.gz.lba
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/db_nucl_fm.fasta.gz.lba
                    URL_HASH SHA256=7a139ea6d9da493ae80a1b0de0060f018cc408ce5b130496320b7f18b6f43eda)

declare_datasource (FILE db_nucl_bifm.fasta.gz.lba
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/db_nucl_bifm.fasta.gz.lba
                    URL_HASH SHA256=4bc3894485ac11210149cb29cbc7a04ec3e8f01d1cc55d50b7c15e61d1a0c115)

declare_datasource (FILE db_nucl_bs_fm.fasta.gz.lba
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/db_nucl_bs_fm.fasta.gz.lba
                    URL_HASH SHA256=e7429e16fa79114228a42e367843719f8bf56374c32537a5c25e5c9a062f3faf)

declare_datasource (FILE db_nucl_bs_bifm.fasta.gz.lba
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/db_nucl_bs_bifm.fasta.gz.lba
                    URL_HASH SHA256=dae02a72f5dd40ccb67f8c4fc8f84b04d148074f84e5e088fff0e61e7c2cda99)

declare_datasource (FILE db_prot_fm.fasta.gz.lba
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/db_prot_fm.fasta.gz.lba
                    URL_HASH SHA256=cc9a563f5a2a2e7fe27b9d19761e8326371625dadf2509289b7a52df65b31da3)

declare_datasource (FILE db_prot_bifm.fasta.gz.lba
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/db_prot_bifm.fasta.gz.lba
                    URL_HASH SHA256=0b2cdbea634481efb840c901ed768bc9bc2811ac87f16cca751c327094c50c75)

declare_datasource (FILE db_trans_fm.fasta.gz.lba
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/db_trans_fm.fasta.gz.lba
                    URL_HASH SHA256=49c0b6784df6de0d6e064c65d034ad3ced961283135f5cc14b86b1048059a0c9)

declare_datasource (FILE db_trans_bifm.fasta.gz.lba
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/db_trans_bifm.fasta.gz.lba
                    URL_HASH SHA256=172f3fb7c529bcf1eff086238efa1c66ce9a609c7abf6146105fa243c3cfcc4d)

## Query input files

declare_datasource (FILE queries_nucl.fasta.gz
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/input_files/queries_nucl.fasta.gz
                    URL_HASH SHA256=7a8adcc7ee5d967992a0624b9fb704118223abd112534c2c6456104b30ddad88)

declare_datasource (FILE queries_nucl_bs.fasta.gz
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/input_files/queries_nucl_bs.fasta.gz
                    URL_HASH SHA256=a358ace3e6f35bd379854fd8870a63a6da62b6bd3e72d95f11e5386398881f0a)

declare_datasource (FILE queries_prot.fasta.gz
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/input_files/queries_prot.fasta.gz
                    URL_HASH SHA256=e21411f422c1dd844696c8ca8a08b93e8243d7f77307731c5f2ee8dcc60d0703)

## Search output files

declare_datasource (FILE output_blastn_fm.bam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_fm.bam
                    URL_HASH SHA256=89e4052abc650d8ab01e927a7e7783de78f8d7e8029f0c406c226e8802176260)

declare_datasource (FILE output_blastn_fm.m0
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_fm.m0
                    URL_HASH SHA256=cac63b12c30ffef52393224ffbe01a17c89ba8b89cb8349acbbac869eb51273b)

declare_datasource (FILE output_blastn_fm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_fm.m8
                    URL_HASH SHA256=2c139f4cf35d3ba5afdf3c0ff7f04c77904d73b7a66895c5d86fcc987b3cf64a)

declare_datasource (FILE output_blastn_fm.m9
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_fm.m9
                    URL_HASH SHA256=5f7a343a5d19238332af16d08b38ea70aa7da5285a8874b361cffe763505a536)

declare_datasource (FILE output_blastn_fm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_fm.sam
                    URL_HASH SHA256=61b1499cead7990c102f7d66ea7323baa1756ce8c832c3caf68d3134df3e4841)

declare_datasource (FILE output_blastn_bs_fm.bam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_bs_fm.bam
                    URL_HASH SHA256=36c07091621de03f8448998adb8dafed041c07ef094a3d14bf2f1fc1f459bb25)

declare_datasource (FILE output_blastn_bs_fm.m0
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_bs_fm.m0
                    URL_HASH SHA256=a565f3bdb3c2617e6a90189fb9c833e619ce2f729248326388a6781b72a5f87c)

declare_datasource (FILE output_blastn_bs_fm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_bs_fm.m8
                    URL_HASH SHA256=37f9f4f93e12c780200511c94418b9e4280b3b8ce3c7008389e89dc1dca9f54b)

declare_datasource (FILE output_blastn_bs_fm.m9
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_bs_fm.m9
                    URL_HASH SHA256=0a6eac342cc4790d3399555dfc2cef4d5804ec058139b26f5d89341309e1c4fc)

declare_datasource (FILE output_blastn_bs_fm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_bs_fm.sam
                    URL_HASH SHA256=0c6fccf892618751187d1195207716abc5b11be3dfbbaf689c3f6696288f3113)

declare_datasource (FILE output_blastp_fm.bam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastp_fm.bam
                    URL_HASH SHA256=a750d396a54ccd7b414f2567d4e96596e9cbd576bb6e2cb40ae7cdf0d049d546)

declare_datasource (FILE output_blastp_fm.m0
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastp_fm.m0
                    URL_HASH SHA256=27c2da9595f3276255eaaaf3910c32809e9472dd4bb190fcea95ebfd14249eb7)

declare_datasource (FILE output_blastp_fm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastp_fm.m8
                    URL_HASH SHA256=a3e2c3c0d305cd8ba2e8a5fe52bb1656e41d3091c24c2b80e7a5a5bb678f401b)

declare_datasource (FILE output_blastp_fm.m9
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastp_fm.m9
                    URL_HASH SHA256=7ee3cc8cf7fb2ea784ec404ec97b9cede14875cfba8de597a16e7b12dc19a1f4)

declare_datasource (FILE output_blastp_fm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastp_fm.sam
                    URL_HASH SHA256=57e2c5c7297a259441976afbc3f4e77c6ea8f92ec8b02e08a103937da8cd32a0)

declare_datasource (FILE output_blastx_fm.bam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastx_fm.bam
                    URL_HASH SHA256=db5daa713193dfcf693cc24e52e4e64631ec58dcaba7eea67901b916d4951843)

declare_datasource (FILE output_blastx_fm.m0
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastx_fm.m0
                    URL_HASH SHA256=0db7f6f44ba0a4cf5052b1c3abf094e374c881ba79e3c9440156cc3f5e7f7c52)

declare_datasource (FILE output_blastx_fm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastx_fm.m8
                    URL_HASH SHA256=bfa3b1d8c0aebea0f0cd2c99fedfe582020c0ac24ab5cd51ecd67649ec75f467)

declare_datasource (FILE output_blastx_fm.m9
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastx_fm.m9
                    URL_HASH SHA256=cd60e7bf440cfc70362d1757842d9cceb5361cfa2e8cbf9223650617dc54dc47)

declare_datasource (FILE output_blastx_fm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastx_fm.sam
                    URL_HASH SHA256=47a7dcc2645221645368acc23d891ca10d04e3f797e0962a7d2734a8415d5d01)

declare_datasource (FILE output_tblastn_fm.bam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastn_fm.bam
                    URL_HASH SHA256=d5dc8d696932d4bc595388368afa60fc6281f5a74be6896e3b6301c048d602af)

declare_datasource (FILE output_tblastn_fm.m0
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastn_fm.m0
                    URL_HASH SHA256=dbe9949ba4c70f30adb9bf22362a09934c735df31b2d81bbb674373ac897ba90)

declare_datasource (FILE output_tblastn_fm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastn_fm.m8
                    URL_HASH SHA256=ec174aad5c9199673cd6b9123bf5e9eada348ca3bc68b3409e36b81cbf533197)

declare_datasource (FILE output_tblastn_fm.m9
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastn_fm.m9
                    URL_HASH SHA256=4498af84a412e273067c3396845b3c8d6a1ef9c72f1c3368207ef6d9ed89b517)

declare_datasource (FILE output_tblastn_fm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastn_fm.sam
                    URL_HASH SHA256=f674bcb028742f89b367ee2c1b59284a83ecbac0c29c793574f744887162552d)

declare_datasource (FILE output_tblastx_fm.bam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastx_fm.bam
                    URL_HASH SHA256=eded153112ea8975a6fe736e74c931bfd212c06fa61c0197dafc7911e2058d82)

declare_datasource (FILE output_tblastx_fm.m0
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastx_fm.m0
                    URL_HASH SHA256=52b19c9a3545471c52a826e95610e67dcf0e6b7851f852a09a54f97b6eaa677b)

declare_datasource (FILE output_tblastx_fm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastx_fm.m8
                    URL_HASH SHA256=a2215049149c581dd34c3c4d28fbb90d654cd59b399b2c7e001af08d858b82c6)

declare_datasource (FILE output_tblastx_fm.m9
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastx_fm.m9
                    URL_HASH SHA256=385a42654fa9ea76ee6dc851702bfad487aa75ff39fab424ec0eec84a2069238)

declare_datasource (FILE output_tblastx_fm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastx_fm.sam
                    URL_HASH SHA256=3b2046496fc651fa1987b666b899c6b6d2fbc62bc6af623f575a7c7e2ef3b32a)

declare_datasource (FILE output_blastn_fm_fast.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_fm_fast.sam
                    URL_HASH SHA256=de57c8ee12447df2d635b1b2e9943bca9d4089e3cba433bb9e78ce4ec7db729f)

declare_datasource (FILE output_blastn_fm_fast.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_fm_fast.m8
                    URL_HASH SHA256=1a7d1a07b8673e7fb6ee0033d27aa705be5e366febb95e8d465947f0680f3521)

declare_datasource (FILE output_blastn_bs_fm_fast.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_bs_fm_fast.sam
                    URL_HASH SHA256=0c6fccf892618751187d1195207716abc5b11be3dfbbaf689c3f6696288f3113)

declare_datasource (FILE output_blastn_bs_fm_fast.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_bs_fm_fast.m8
                    URL_HASH SHA256=37f9f4f93e12c780200511c94418b9e4280b3b8ce3c7008389e89dc1dca9f54b)

declare_datasource (FILE output_blastp_fm_fast.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastp_fm_fast.sam
                    URL_HASH SHA256=410ff4b7316dd8353cfdc42b68eb01bbdda189b6bbfbbe3fc8333eb7d80328ba)

declare_datasource (FILE output_blastp_fm_fast.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastp_fm_fast.m8
                    URL_HASH SHA256=47fce11a96fbcce92f62bfacb1545de432e4d97f1d5b74c4d1e6f3fd9418898e)

declare_datasource (FILE output_blastx_fm_fast.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastx_fm_fast.sam
                    URL_HASH SHA256=512373692eb80bf426f027dac5a4355c6f94589cc60bd97d16e9964bc5d8689c)

declare_datasource (FILE output_blastx_fm_fast.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastx_fm_fast.m8
                    URL_HASH SHA256=b039b8fabbcfb67bdd2d92c99c4db161f820bc9ea1886ac088327c0b3803f3e3)

declare_datasource (FILE output_tblastn_fm_fast.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastn_fm_fast.sam
                    URL_HASH SHA256=d74c48bb71ed0cff236d3f9b14ae694155d1e798f309a84113a3748fba24605d)

declare_datasource (FILE output_tblastn_fm_fast.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastn_fm_fast.m8
                    URL_HASH SHA256=ae3f8e5f2966579c09688667d56f98b48a5ac5cd513d408cd2de09f7090aa1d9)

declare_datasource (FILE output_tblastx_fm_fast.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastx_fm_fast.sam
                    URL_HASH SHA256=ffae5e644c963294766a41879fb990c11ec308e56df123649ea2ce9d677eaa93)

declare_datasource (FILE output_tblastx_fm_fast.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastx_fm_fast.m8
                    URL_HASH SHA256=686450a6a1ff546244a75bfb05925963bbce79ecf9a4c67a3566d968134cd0e7)

declare_datasource (FILE output_blastn_fm_sensitive.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_fm_sensitive.sam
                    URL_HASH SHA256=34e6942316c33300004f2d28cc8644dd4bd6cc73f18e70d3a0a357f74519e2ea)

declare_datasource (FILE output_blastn_fm_sensitive.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_fm_sensitive.m8
                    URL_HASH SHA256=16b36d8cab68aa08ca6130be8b91b4c9f48411eaf6edcf4d97b87bde5fe76134)

declare_datasource (FILE output_blastn_bs_fm_sensitive.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_bs_fm_sensitive.sam
                    URL_HASH SHA256=c26dfe05d85a7dae1ea13949d8e93a874d817620973e763ca482e135b54485c1)

declare_datasource (FILE output_blastn_bs_fm_sensitive.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastn_bs_fm_sensitive.m8
                    URL_HASH SHA256=3d47f9c25332d116a9d1ef8ed438beda2b3ef47691adcb63ff0f9e6347698535)

declare_datasource (FILE output_blastp_fm_sensitive.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastp_fm_sensitive.sam
                    URL_HASH SHA256=582fe1d78cc03eb94e70b78ab8275d7c34a597bb7a896277b5a93573c489dfa8)

declare_datasource (FILE output_blastp_fm_sensitive.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastp_fm_sensitive.m8
                    URL_HASH SHA256=0fe2075c016e847b9ffb2b817e0c59326228d2fbeaef86af5542ad187f4ba283)

declare_datasource (FILE output_blastx_fm_sensitive.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastx_fm_sensitive.sam
                    URL_HASH SHA256=0009d2e5cb5019d8c9867b898f9cd2593306069c6b0bbad33b97914cef0d741b)

declare_datasource (FILE output_blastx_fm_sensitive.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_blastx_fm_sensitive.m8
                    URL_HASH SHA256=93d341c0cf30a9f42221008fad3ff46c1d213f3c89ce1168871d259b02824260)

declare_datasource (FILE output_tblastn_fm_sensitive.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastn_fm_sensitive.sam
                    URL_HASH SHA256=f1eb59c90b051798ef0c18f4000c0354fe8c97f83c384f614a19d1d294f7aab6)

declare_datasource (FILE output_tblastn_fm_sensitive.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastn_fm_sensitive.m8
                    URL_HASH SHA256=292ec38ced08c68eca49dd0cd08bc0d89d9c011cddb7ccba445c9c61708cc51b)

declare_datasource (FILE output_tblastx_fm_sensitive.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastx_fm_sensitive.sam
                    URL_HASH SHA256=727edd5f9c1794ba61af90ccbd3008f1c736383031591a3c76537a55333fcfa3)

declare_datasource (FILE output_tblastx_fm_sensitive.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/b2e625245dfcb9119f0eea01c96cd5dd59591b54/output_files/output_tblastx_fm_sensitive.m8
                    URL_HASH SHA256=fb0109e002e6ef9805a3aa270ec0455e3c67ceed34a0a297fb6d463472086005)
