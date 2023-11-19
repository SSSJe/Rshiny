## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----basic-demo-code, eval=FALSE----------------------------------------------
#  library(shiny)
#  library(bs4Dash)
#  
#  shinyApp(
#      ui = dashboardPage(
#          header = dashboardHeader(
#              title = dashboardBrand(
#                  title = "My dashboard",
#                  color = "primary",
#                  href = "https://adminlte.io/themes/v3",
#                  image = "https://adminlte.io/themes/v3/dist/img/AdminLTELogo.png"
#              )
#          ),
#          sidebar = dashboardSidebar(),
#          body = dashboardBody(
#              lapply(getAdminLTEColors(), function(color) {
#                  box(status = color)
#              })
#          ),
#          controlbar = dashboardControlbar(),
#          title = "DashboardPage"
#      ),
#      server = function(input, output) { }
#  )

## ---- echo=FALSE--------------------------------------------------------------
bs4Dash:::app_container("https://dgranjon.shinyapps.io/bs4Dash-scrollToTop/", TRUE)

