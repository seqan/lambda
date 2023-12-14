cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

## Index input files

set(BASEURL "https://raw.githubusercontent.com/h-2/lambda-testdata/832453c8721094af511a2d7041dbb107a9ecc912")

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

declare_datasource (FILE output_blastn_bs_fm_fast.m8
                    URL ${BASEURL}/output_files/output_blastn_bs_fm_fast.m8
                    URL_HASH SHA256=732b439bded780e0b4b68478e6e89368934cab43d261c2f41e48d90ab505a01e)
declare_datasource (FILE output_blastn_bs_fm_fast.sam
                    URL ${BASEURL}/output_files/output_blastn_bs_fm_fast.sam
                    URL_HASH SHA256=2247d450da5ad0fa76f630e373bbf70d0fde79a736637d19f8f25d985178fcfa)
declare_datasource (FILE output_blastn_bs_fm_pairs_default.m8
                    URL ${BASEURL}/output_files/output_blastn_bs_fm_pairs_default.m8
                    URL_HASH SHA256=11c0e206f0c48bbe87bb464c175bd6bea3ce25c02769d33961beeadba15fcc98)
declare_datasource (FILE output_blastn_bs_fm_pairs_default.sam
                    URL ${BASEURL}/output_files/output_blastn_bs_fm_pairs_default.sam
                    URL_HASH SHA256=55e2f168ab0a3fd2d1d55e8344473ab67c75c509f8e42c3fa3d0cb5d655141cf)
declare_datasource (FILE output_blastn_bs_fm_pairs_sensitive.m8
                    URL ${BASEURL}/output_files/output_blastn_bs_fm_pairs_sensitive.m8
                    URL_HASH SHA256=1bd31e06f5c20c273b5079a8f483ae574daa6c8e38b50cdf45bacad18f940901)
declare_datasource (FILE output_blastn_bs_fm_pairs_sensitive.sam
                    URL ${BASEURL}/output_files/output_blastn_bs_fm_pairs_sensitive.sam
                    URL_HASH SHA256=ff854acd2bba55142275ace681c9270d2e287e3be7123f7ce69598cde2eb98ca)
declare_datasource (FILE output_blastn_bs_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_blastn_bs_fm_sensitive.m8
                    URL_HASH SHA256=5836c6941f04119823ec9e23ebf7cfec1c9b1e59d8cc486211b0aeb1d34ca0ca)
declare_datasource (FILE output_blastn_bs_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_blastn_bs_fm_sensitive.sam
                    URL_HASH SHA256=bee4b36f98572c6a4d6cc4f9c6cacf40b627b57dc8c5964d82921c9e7bd884c8)
declare_datasource (FILE output_blastn_bs_fm.bam
                    URL ${BASEURL}/output_files/output_blastn_bs_fm.bam
                    URL_HASH SHA256=cb8204f66fe135f39a01a495d553a4d2ebb7f190aab1cc92fc78e158d8ebb6ae)
declare_datasource (FILE output_blastn_bs_fm.m0
                    URL ${BASEURL}/output_files/output_blastn_bs_fm.m0
                    URL_HASH SHA256=20227f3db16bf6dd5ecb5c22d831eb42e976f664df3678a49670369366ea3e00)
declare_datasource (FILE output_blastn_bs_fm.m8
                    URL ${BASEURL}/output_files/output_blastn_bs_fm.m8
                    URL_HASH SHA256=732b439bded780e0b4b68478e6e89368934cab43d261c2f41e48d90ab505a01e)
declare_datasource (FILE output_blastn_bs_fm.m9
                    URL ${BASEURL}/output_files/output_blastn_bs_fm.m9
                    URL_HASH SHA256=0bdcd2dae96480227f4da0dd431118ba0564c923f1fe9ea1f3eaf1cf143af44f)
declare_datasource (FILE output_blastn_bs_fm.sam
                    URL ${BASEURL}/output_files/output_blastn_bs_fm.sam
                    URL_HASH SHA256=2247d450da5ad0fa76f630e373bbf70d0fde79a736637d19f8f25d985178fcfa)
declare_datasource (FILE output_blastn_fm_fast.m8
                    URL ${BASEURL}/output_files/output_blastn_fm_fast.m8
                    URL_HASH SHA256=8486e576c6ce0ed0f605f20f794ec2a7804cdb302c8261ab975fdd3b8f8f02b0)
declare_datasource (FILE output_blastn_fm_fast.sam
                    URL ${BASEURL}/output_files/output_blastn_fm_fast.sam
                    URL_HASH SHA256=105b00b92787338dbd3b2bed3322d452b75427f66d9cdc6c71028361afccc623)
declare_datasource (FILE output_blastn_fm_pairs_default.m8
                    URL ${BASEURL}/output_files/output_blastn_fm_pairs_default.m8
                    URL_HASH SHA256=f3ce1c515a7e48d1594ef7de7ec3af03fed352b111dfb387304012dda5d33b53)
declare_datasource (FILE output_blastn_fm_pairs_default.sam
                    URL ${BASEURL}/output_files/output_blastn_fm_pairs_default.sam
                    URL_HASH SHA256=f9edd30876b9a75d99d4bf52b71f6da22038de43ba33f711805044ad8d1d5794)
declare_datasource (FILE output_blastn_fm_pairs_sensitive.m8
                    URL ${BASEURL}/output_files/output_blastn_fm_pairs_sensitive.m8
                    URL_HASH SHA256=367660b88239aaa5ed32dbce7702e3c22eb5adb7c6f1c6f043a8e881c930e872)
declare_datasource (FILE output_blastn_fm_pairs_sensitive.sam
                    URL ${BASEURL}/output_files/output_blastn_fm_pairs_sensitive.sam
                    URL_HASH SHA256=e65eab36c18315d830d1abcf6b052657db113797e7a3cd7db7fc095b34a00b54)
declare_datasource (FILE output_blastn_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_blastn_fm_sensitive.m8
                    URL_HASH SHA256=e173961e41df452e67dff97e3eeee8dfaec6cdfa3b2250d2c6772c06bb2801de)
declare_datasource (FILE output_blastn_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_blastn_fm_sensitive.sam
                    URL_HASH SHA256=9737c9d2fa1f82f851a614fc84d91298e6b2c33b1809ea787708077c3e47c803)
declare_datasource (FILE output_blastn_fm.bam
                    URL ${BASEURL}/output_files/output_blastn_fm.bam
                    URL_HASH SHA256=989ecf8caddc2bbfe84c1c47d74fe88752027effc19e02b8ad30308022835996)
declare_datasource (FILE output_blastn_fm.m0
                    URL ${BASEURL}/output_files/output_blastn_fm.m0
                    URL_HASH SHA256=b533756dc14ee58fd350f5fc4f3ed8dac38dba8100185f3e1a23800fb06cba01)
declare_datasource (FILE output_blastn_fm.m8
                    URL ${BASEURL}/output_files/output_blastn_fm.m8
                    URL_HASH SHA256=18b7a0feb4b5e76a44be7ec4ae26ceb1be05fa257df0725a6d231f841de01d54)
declare_datasource (FILE output_blastn_fm.m9
                    URL ${BASEURL}/output_files/output_blastn_fm.m9
                    URL_HASH SHA256=f20ecfabbaf801ec7532b956931b36cbcaa4de5f9a946e139597334dfb50c132)
declare_datasource (FILE output_blastn_fm.sam
                    URL ${BASEURL}/output_files/output_blastn_fm.sam
                    URL_HASH SHA256=ca8202a68e22d59e731d023c021244e27e664077d8aa68662041346ac828c110)
declare_datasource (FILE output_blastp_fm_fast.m8
                    URL ${BASEURL}/output_files/output_blastp_fm_fast.m8
                    URL_HASH SHA256=753acc7afccbcacc8e585e7a83b7b295548ea536b23800ca90a9566bcb7a48cb)
declare_datasource (FILE output_blastp_fm_fast.sam
                    URL ${BASEURL}/output_files/output_blastp_fm_fast.sam
                    URL_HASH SHA256=36bca3f2ee56d2234ed44f8351ad83b5fe267c9e75a9eea5c13d34565f20f2ac)
declare_datasource (FILE output_blastp_fm_pairs_default.m8
                    URL ${BASEURL}/output_files/output_blastp_fm_pairs_default.m8
                    URL_HASH SHA256=8dd2211fbdf3dcbdbd5a53ce02f66d4cbe37ffe5543ba61d3cc324f5c6b698a8)
declare_datasource (FILE output_blastp_fm_pairs_default.sam
                    URL ${BASEURL}/output_files/output_blastp_fm_pairs_default.sam
                    URL_HASH SHA256=f1f874880d6ba5280f1fcb2c0e5ff2746ed343a7e9c16e49e4d0db4b48e32041)
declare_datasource (FILE output_blastp_fm_pairs_sensitive.m8
                    URL ${BASEURL}/output_files/output_blastp_fm_pairs_sensitive.m8
                    URL_HASH SHA256=9c9bb4c98315f1c3fee2116bc565b0eb9abc061af09282680275118a0c884fcf)
declare_datasource (FILE output_blastp_fm_pairs_sensitive.sam
                    URL ${BASEURL}/output_files/output_blastp_fm_pairs_sensitive.sam
                    URL_HASH SHA256=e6d733403afc9333429a8617ce91a326861b9c82cb676a819d8c281c7e1fdcea)
declare_datasource (FILE output_blastp_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_blastp_fm_sensitive.m8
                    URL_HASH SHA256=25b8fa3de4eda51b016baa90ade64f5ff4f409b37b32b16226e46ca350673823)
declare_datasource (FILE output_blastp_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_blastp_fm_sensitive.sam
                    URL_HASH SHA256=83fc52cfe429c8ac4d790f69f32856f1426ed003c128d6dc4c8767e977e5e9c1)
declare_datasource (FILE output_blastp_fm.bam
                    URL ${BASEURL}/output_files/output_blastp_fm.bam
                    URL_HASH SHA256=322932054405c802211f9fd1f9e86a96c8b19ea85d2eaa8cbd18f4bcfb1d5131)
declare_datasource (FILE output_blastp_fm.m0
                    URL ${BASEURL}/output_files/output_blastp_fm.m0
                    URL_HASH SHA256=56e29329f35f1868d6c5a2ab0c88e7ea1446d433f5d7bb70308ceb97599fcd5f)
declare_datasource (FILE output_blastp_fm.m8
                    URL ${BASEURL}/output_files/output_blastp_fm.m8
                    URL_HASH SHA256=99f520bb55f5c1b371ae5ba41b881c21360a84e7fcb8ff097e01036beb801d4c)
declare_datasource (FILE output_blastp_fm.m9
                    URL ${BASEURL}/output_files/output_blastp_fm.m9
                    URL_HASH SHA256=3d7b519b26419361c37f5512b28220f3004620aafa2aae62f0e13bbee3ee92f2)
declare_datasource (FILE output_blastp_fm.sam
                    URL ${BASEURL}/output_files/output_blastp_fm.sam
                    URL_HASH SHA256=d7267490e2c6c1ebf838fb7f08a5f370bf0f2b23024ca76db6bfed440851fc91)
declare_datasource (FILE output_blastx_fm_fast.m8
                    URL ${BASEURL}/output_files/output_blastx_fm_fast.m8
                    URL_HASH SHA256=5554e8dcd4de9ef02818e18ad8fb185594762a425aded7929a0b2a574137f46e)
declare_datasource (FILE output_blastx_fm_fast.sam
                    URL ${BASEURL}/output_files/output_blastx_fm_fast.sam
                    URL_HASH SHA256=557394a3c75c2da49d3d15c34f10f1eca2e4a934102cfbe92081df0f704b1828)
declare_datasource (FILE output_blastx_fm_pairs_default.m8
                    URL ${BASEURL}/output_files/output_blastx_fm_pairs_default.m8
                    URL_HASH SHA256=9479d7bb3ca8e53b1d425b70867f0d4aeecd73026596997724f4b1d3458ea16f)
declare_datasource (FILE output_blastx_fm_pairs_default.sam
                    URL ${BASEURL}/output_files/output_blastx_fm_pairs_default.sam
                    URL_HASH SHA256=44c8a52fa3c76ca943de5dc47e8a20a3f09f70f3c513096f1bc1101b822699e4)
declare_datasource (FILE output_blastx_fm_pairs_sensitive.m8
                    URL ${BASEURL}/output_files/output_blastx_fm_pairs_sensitive.m8
                    URL_HASH SHA256=abdd532001a794d3abbe3411e842271eae503684bca26c259c0170c2d6da0a61)
declare_datasource (FILE output_blastx_fm_pairs_sensitive.sam
                    URL ${BASEURL}/output_files/output_blastx_fm_pairs_sensitive.sam
                    URL_HASH SHA256=49b24eb192c7dbcee68662ba789ecd772a7160229e296653f7960e1ab639f1cf)
declare_datasource (FILE output_blastx_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_blastx_fm_sensitive.m8
                    URL_HASH SHA256=59f22a6158b48463b408b329db9472412838455dfbd74a6df7ab8ed4a0240914)
declare_datasource (FILE output_blastx_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_blastx_fm_sensitive.sam
                    URL_HASH SHA256=1402fe5c5161be0df142b279ef293f6edb6685ace4f85c2a6e8b5ff43e41a8b8)
declare_datasource (FILE output_blastx_fm.bam
                    URL ${BASEURL}/output_files/output_blastx_fm.bam
                    URL_HASH SHA256=7a4cad866e89ef84839bba6f9eece7f20528758075fc7faac2419c9406be84cc)
declare_datasource (FILE output_blastx_fm.m0
                    URL ${BASEURL}/output_files/output_blastx_fm.m0
                    URL_HASH SHA256=ad5d8ed1fa3ea34edc9e3f5f0c256a7be83db7910fb45a65b514bb506a6b91fe)
declare_datasource (FILE output_blastx_fm.m8
                    URL ${BASEURL}/output_files/output_blastx_fm.m8
                    URL_HASH SHA256=d3bb3a5f101e9deab0b7605b83d3c2ce03270545923a65d10b0c7afb52c324b7)
declare_datasource (FILE output_blastx_fm.m9
                    URL ${BASEURL}/output_files/output_blastx_fm.m9
                    URL_HASH SHA256=daabecc9b0ef1dd9f27ac7a27baf00327043b73ee94f3dc0a398a3ef2faebc02)
declare_datasource (FILE output_blastx_fm.sam
                    URL ${BASEURL}/output_files/output_blastx_fm.sam
                    URL_HASH SHA256=78f55f6dbff1205e09efa00b65498cb0de5e1d907be394cc8b73e75be7069a3e)
declare_datasource (FILE output_tblastn_fm_fast.m8
                    URL ${BASEURL}/output_files/output_tblastn_fm_fast.m8
                    URL_HASH SHA256=e7e34710b246347dfdc12619b064da1ec0f1d248316901209d89672e92d9bbb2)
declare_datasource (FILE output_tblastn_fm_fast.sam
                    URL ${BASEURL}/output_files/output_tblastn_fm_fast.sam
                    URL_HASH SHA256=e416ddf4d79ea72f8240d291c17c5e41fc2a2d7b4560f1c9ad600b498f9e4d1a)
declare_datasource (FILE output_tblastn_fm_pairs_default.m8
                    URL ${BASEURL}/output_files/output_tblastn_fm_pairs_default.m8
                    URL_HASH SHA256=9a8416ca080499af6c91e63d6c6ee5315a5440177dc8fc167c490d3d0df05c7c)
declare_datasource (FILE output_tblastn_fm_pairs_default.sam
                    URL ${BASEURL}/output_files/output_tblastn_fm_pairs_default.sam
                    URL_HASH SHA256=06d8a2fefc4bfb3fc39605f06bf28b96c442db25c228f59e24fc5c51c68597dd)
declare_datasource (FILE output_tblastn_fm_pairs_sensitive.m8
                    URL ${BASEURL}/output_files/output_tblastn_fm_pairs_sensitive.m8
                    URL_HASH SHA256=877bb8dbed8a69a065fa82dd3af3f68e393db5ccda575ab62f1ac19c69cb08ee)
declare_datasource (FILE output_tblastn_fm_pairs_sensitive.sam
                    URL ${BASEURL}/output_files/output_tblastn_fm_pairs_sensitive.sam
                    URL_HASH SHA256=5cc60d9032f4f8f7de1643d1f5181534288af78969b77051640efc0a7e1ab852)
declare_datasource (FILE output_tblastn_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_tblastn_fm_sensitive.m8
                    URL_HASH SHA256=c07e11add4b48d1d910c76f95a007de452c16a8fe8ea83f695a1c75ee4152db1)
declare_datasource (FILE output_tblastn_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_tblastn_fm_sensitive.sam
                    URL_HASH SHA256=b9b54e58f158746e1d69fc24d841ba637c859849bc9fe50e2c75cf3bc634a3e8)
declare_datasource (FILE output_tblastn_fm.bam
                    URL ${BASEURL}/output_files/output_tblastn_fm.bam
                    URL_HASH SHA256=0d646463d379eeadca4af6705081854034d2e07bce870dee1165f11dc50f9336)
declare_datasource (FILE output_tblastn_fm.m0
                    URL ${BASEURL}/output_files/output_tblastn_fm.m0
                    URL_HASH SHA256=d3922d4db0245d3c797b4ca5a982e73c3331be7b7f3282754e31daf2820691aa)
declare_datasource (FILE output_tblastn_fm.m8
                    URL ${BASEURL}/output_files/output_tblastn_fm.m8
                    URL_HASH SHA256=316ff6d392d1553bf3115f41e14e5cfced72f6b22459e487b0def80c37123b44)
declare_datasource (FILE output_tblastn_fm.m9
                    URL ${BASEURL}/output_files/output_tblastn_fm.m9
                    URL_HASH SHA256=720895bca93a35543fc58f2dbc490ba9ce531195fe71d0dbbc0ed41e6392002f)
declare_datasource (FILE output_tblastn_fm.sam
                    URL ${BASEURL}/output_files/output_tblastn_fm.sam
                    URL_HASH SHA256=b2ca3d794d03154ed956d095e531858d71aac8d81b4fdee5ce18ed585300e034)
declare_datasource (FILE output_tblastx_fm_fast.m8
                    URL ${BASEURL}/output_files/output_tblastx_fm_fast.m8
                    URL_HASH SHA256=dea0901a96f8dafd9e3f7518223257ac6ad0b8e4bf2ec42fe05701707ed72b4a)
declare_datasource (FILE output_tblastx_fm_fast.sam
                    URL ${BASEURL}/output_files/output_tblastx_fm_fast.sam
                    URL_HASH SHA256=3c23be9c52c1f05acef72acbfc2ea316e18a9e3fa9a9ad1f59903118783c4129)
declare_datasource (FILE output_tblastx_fm_pairs_default.m8
                    URL ${BASEURL}/output_files/output_tblastx_fm_pairs_default.m8
                    URL_HASH SHA256=8539cce5bb7ca18974eccf060b361fc09271449add941b21a99b040b001c4683)
declare_datasource (FILE output_tblastx_fm_pairs_default.sam
                    URL ${BASEURL}/output_files/output_tblastx_fm_pairs_default.sam
                    URL_HASH SHA256=0ebc191fac3f01ab6d72dd29729e49d0b06da3329281209165abd6a0260af02e)
declare_datasource (FILE output_tblastx_fm_pairs_sensitive.m8
                    URL ${BASEURL}/output_files/output_tblastx_fm_pairs_sensitive.m8
                    URL_HASH SHA256=59d760626164cdb2564686ce30c36d71c881f657eebca8b2a06ebe3b9a70d50d)
declare_datasource (FILE output_tblastx_fm_pairs_sensitive.sam
                    URL ${BASEURL}/output_files/output_tblastx_fm_pairs_sensitive.sam
                    URL_HASH SHA256=5bf3a0fe8ae736b072fd99df764c2039b2dcd28c4486b21c3796b5712d05c1fc)
declare_datasource (FILE output_tblastx_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_tblastx_fm_sensitive.m8
                    URL_HASH SHA256=a590cff9711dbbad654a7317c081af5ffacda9c934ece3bf7b59eec223c2cf1f)
declare_datasource (FILE output_tblastx_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_tblastx_fm_sensitive.sam
                    URL_HASH SHA256=8b400bf75ba657dd29983fcf0e5f0f8de2243d8b68287af8b884aab61524ba2e)
declare_datasource (FILE output_tblastx_fm.bam
                    URL ${BASEURL}/output_files/output_tblastx_fm.bam
                    URL_HASH SHA256=ab3d60dc3beaa6d5ac78aaa7c6bd212297c9e7ab982fc9c9f6d3c69c6b819b87)
declare_datasource (FILE output_tblastx_fm.m0
                    URL ${BASEURL}/output_files/output_tblastx_fm.m0
                    URL_HASH SHA256=fcb43c8084305372ebaf931736bc1dc61a1589021aa3d36ff1d0610e7a423cbb)
declare_datasource (FILE output_tblastx_fm.m8
                    URL ${BASEURL}/output_files/output_tblastx_fm.m8
                    URL_HASH SHA256=d040e8debb4e5cff6443e676e77449a65f4ffdb83924660f339d640b121964f1)
declare_datasource (FILE output_tblastx_fm.m9
                    URL ${BASEURL}/output_files/output_tblastx_fm.m9
                    URL_HASH SHA256=970a2a628e1816243645e6eabc26ee8f20f2f2ae447e5d18b432fb125d25e867)
declare_datasource (FILE output_tblastx_fm.sam
                    URL ${BASEURL}/output_files/output_tblastx_fm.sam
                    URL_HASH SHA256=6f218b39e8785e1bebbe8be486237e03e06d9fccd41a5e57e575539cfe4eab6e)
