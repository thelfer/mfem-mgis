# Installation guide

The installation guide is available at
<https://thelfer.github.io/mfem-mgis/installation_guide/installation_guide.html>.

In short, `mfem-mgis` is packaged in [Spack](https://spack.io/), which also
installs all its dependencies. These commands set up Spack:

```sh
git clone --depth=2 --branch=v1.2.2 https://github.com/spack/spack.git
source spack/share/spack/setup-env.sh
spack repo update builtin --branch develop
```

The last command selects the `develop` branch of the Spack packages, because
their releases do not provide `mfem-mgis` yet.

The `@master` version of `mfem-mgis`, its development version, is currently
recommended, as it includes many fixes:

```sh
spack install mfem-mgis@master
```

The last release (1.0.4) is installed by:

```sh
spack install mfem-mgis
```
