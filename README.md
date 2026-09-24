# dbdahl.r-universe.dev

R-universe registry for
[dbdahl.r-universe.dev](https://dbdahl.r-universe.dev).

## How packages get here

`packages.json` points at each package's own public GitHub repository.
R-universe builds and checks the default branch automatically after each push.

## Installing packages

```r
install.packages(
  "salso",
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

Then `install.packages("salso")` works normally.

Always include CRAN alongside the universe URL so that dependencies resolve
correctly.

## Latest build results

[Latest build results](https://github.com/r-universe/dbdahl/actions)
