## Taking GLF as input
../bin/polymutt -p test.ped -d test.dat -g test.gif -c 0.9 --minDepth 150 --maxDepth 200 --nthreads 4  --out_vcf test.out.vcf

## Taking VCF as input
../bin/polymutt -p test.ped -d test.dat --in_vcf testvcf.in.vcf --out_vcf testvcf.out.vcf

## Mixture of related and unrelated individuals
../bin/polymutt -p test.mix.ped -d test.dat -g test.gif  --out_vcf test.mix.out.vcf

## Calling de novo mutations
../bin/polymutt -p test.ped -d test.dat -g test.gif --nthreads 4  --out_vcf test.denovo.out.vcf --denovo --rate_denovo 1.5e-07
