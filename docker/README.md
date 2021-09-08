# homelette Docker Image

Since a full installation of all dependencies for `homelette` currently requires a personal license for MODELLER, we are only able to share a docker image that contains everything **expect** the MODELLER license key. We have provided here scripts that, given a valid MODELLER license key, assemble a **local** image with all dependencies available. **Please be aware that this image contains your MODELLER license key and should not be shared on DockerHub or other repositries!**

A license key for MODELLER for academic use can be easily aquired [here](https://salilab.org/modeller/registration.html). 

## Constructing local homelette image

The local, full installation of homelette is based on the `homelette_template:latest` image shared on DockerHub. 

A bash script (`construct_homelette_image.sh`) has been provided. This script pulls the latest version of the `homelette_template` image from DockerHub and then attempts to construct the local homelette image with the given MODELLER license key:

```bash
./construct_homelette_image.sh "YOUR MODELLER KEY HERE"
```


## Using the local homelette image

After constructing the local `homelette:latest` image, you can use it as every other docker image. In order to make access somewhat more comfortable, we have developed a small bash script that makes access to the local image in different modes easier: `homelette.sh`

For a full documentation, please run `./homelette.sh -h`

---

For more details, please visit our [documentation](https://homelette.readthedocs.io/).
