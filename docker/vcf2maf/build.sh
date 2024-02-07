# Build image
#docker build -f Dockerfile --tag=ccbr_vcf2maf:v102.0.0 .
docker build -t=ccbr_vcf2maf:v102.0.0 .

# Tag image with version and reset latest
#docker tag ccbr_vcf2maf:v102.0.0 skchronicles/ccbr_vcf2maf:v102.0.0
#docker tag ccbr_vcf2maf:v102.0.0 skchronicles/ccbr_vcf2maf
#docker tag ccbr_vcf2maf:v102.0.0 nciccbr/ccbr_vcf2maf:v102.0.0
#docker tag ccbr_vcf2maf:v102.0.0 nciccbr/ccbr_vcf2maf

#docker tag ccbr_vcf2maf:v106.0.0 nciccbr/ccbr_vcf2maf
docker tag ccbr_vcf2maf:v102.0.0 dnousome/ccbr_vcf2maf:v102.0.0
docker tag ccbr_vcf2maf:v102.0.0 dnousome/ccbr_vcf2maf

# Push image to DockerHub
#docker push skchronicles/ccbr_vcf2maf:v102.0.0
#docker push skchronicles/ccbr_vcf2maf:latest
#docker push nciccbr/ccbr_vcf2maf:v102.0.0
#docker push nciccbr/ccbr_vcf2maf:latest

docker push dnousome/ccbr_vcf2maf:v102.0.0
docker push dnousome/ccbr_vcf2maf:latest
