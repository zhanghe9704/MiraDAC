# JLL runtime tests

Upload the `MiraDAC.v<version>.<platform>.tar.gz` files from BinaryBuilder as
artifacts of an Actions run **in this repository** (one artifact containing all
tarballs, or separate artifacts), for example with a producer workflow step:

```yaml
- uses: actions/upload-artifact@v4
  with:
    name: miradac-tarballs
    path: products/MiraDAC.v*.tar.gz
```

Adjust `path` to the recipe's actual tarball output directory. On the
repository's **Actions → JLL runtime tests → Run workflow** page, select the
branch containing this workflow and
enter that artifact-producing run's numeric ID as `run_id`. The ID is the number
at the end of the producer run's Actions URL. The workflow downloads its
artifacts, tests each available platform on Julia 1.10 and 1, and logs `SKIP`
for missing platforms.

The branch under test must include MiraDAC's `MiraDAC_jll` dependency and
no `libmiradac` preference in `julia/MiraDAC/LocalPreferences.toml`. This
workflow creates a temporary local JLL pointing at the extracted binary,
so it can run before `MiraDAC_jll` is registered. It does not publish a JLL.
