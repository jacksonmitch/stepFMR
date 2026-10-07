
.onAttach <- function(libname, pkgname) {
  logo <- paste0(
"         __             ________  _______
   _____/ /____  ____  / ____/  |/  / __ \\
  / ___/ __/ _ \\/ __ \\/ /_  / /|_/ / /_/ /
 (__  ) /_/  __/ /_/ / __/ / /  / / _, _/
/____/\\__/\\___/ .___/_/   /_/  /_/_/ |_|
             /_/         version ", utils::packageVersion("stepFMR"))
  packageStartupMessage(
    "==========================================="
  )
  packageStartupMessage(logo)
  packageStartupMessage(
    "==========================================="
  )
}
