# dbdahl.r-universe.dev

R-universe registry and public package staging for
[dbdahl.r-universe.dev](https://dbdahl.r-universe.dev).

## Branch layout

- **`main`** — registry (`packages.json`) and tooling only
- **`pkg/<package>`** — published source snapshot for each package

## Installing packages

```r
install.packages(
  "splinclust",
  repos = c("https://dbdahl.r-universe.dev", "https://cloud.r-project.org")
)
```

For repeated use, add to `~/.Rprofile`:

```r
options(repos = c(
  dbdahl = "https://dbdahl.r-universe.dev",
  CRAN   = "https://cloud.r-project.org"
))
```

Then `install.packages("splinclust")` works normally.

Always include CRAN alongside the universe URL so that dependencies resolve
correctly.
