## ----setup, include=FALSE-----------------------------------------------------
library(shiny)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----skinSelector-code, eval=FALSE--------------------------------------------
#  library(shiny)
#  library(bs4Dash)
#  shinyApp(
#    ui = dashboardPage(
#      header = dashboardHeader(),
#      sidebar = dashboardSidebar(
#        sidebarMenu(
#          menuItem(
#            text = "Item 1"
#          ),
#          menuItem(
#            text = "Item 2"
#          )
#        )
#      ),
#      body = dashboardBody(),
#      controlbar = dashboardControlbar(
#        collapsed = FALSE,
#        div(class = "p-3", skinSelector()),
#        pinned = TRUE
#      ),
#      title = "Skin Selector"
#    ),
#    server = function(input, output) { }
#  )

## ---- echo=FALSE--------------------------------------------------------------
# App output
bs4Dash:::app_container("https://dgranjon.shinyapps.io/bs4Dash-skinSelector/", TRUE)

## ----fresh-code, eval=FALSE---------------------------------------------------
#  library(fresh)
#  # create the theme with a cyberpunk color palette
#  theme <- create_theme(
#    bs4dash_vars(
#      navbar_light_color = "#bec5cb",
#      navbar_light_active_color = "#FFF",
#      navbar_light_hover_color = "#FFF"
#    ),
#    bs4dash_yiq(
#      contrasted_threshold = 10,
#      text_dark = "#FFF",
#      text_light = "#272c30"
#    ),
#    bs4dash_layout(
#      main_bg = "#353c42"
#    ),
#    bs4dash_sidebar_light(
#      bg = "#272c30",
#      color = "#bec5cb",
#      hover_color = "#FFF",
#      submenu_bg = "#272c30",
#      submenu_color = "#FFF",
#      submenu_hover_color = "#FFF"
#    ),
#    bs4dash_status(
#      primary = "#5E81AC", danger = "#BF616A", light = "#272c30"
#    ),
#    bs4dash_color(
#      gray_900 = "#FFF", white = "#272c30"
#    )
#  )
#  
#  # create tribble for box global config
#  box_config <- tibble::tribble(
#    ~background, ~labelStatus,
#    "danger", "warning",
#    "purple", "success",
#    "success", "primary",
#    "warning", "danger",
#    "fuchsia", "info"
#  )
#  
#  # box factory function
#  box_factory <- function(background, labelStatus) {
#    box(
#      title = "Cyberpunk Box",
#      collapsible = TRUE,
#      background = background,
#      height = "200px",
#      label = boxLabel(1, labelStatus)
#    )
#  }
#  
#  # pmap magic
#  boxes <- purrr::pmap(box_config, box_factory)
#  
#  shinyApp(
#    ui = dashboardPage(
#      freshTheme = theme,
#      header = dashboardHeader(
#        leftUi = dropdownMenu(
#          type = "messages",
#          badgeStatus = "success",
#          messageItem(
#            from = "Support Team",
#            message = "This is the content of a message.",
#            time = "5 mins"
#          ),
#          messageItem(
#            from = "Support Team",
#            message = "This is the content of another message.",
#            time = "2 hours"
#          )
#        )
#      ),
#      sidebar = dashboardSidebar(),
#      body = dashboardBody(boxes),
#      controlbar = dashboardControlbar(),
#      title = "Fresh theming"
#    ),
#    server = function(input, output) { }
#  )

## ---- echo=FALSE--------------------------------------------------------------
# App output
bs4Dash:::app_container("https://dgranjon.shinyapps.io/bs4Dash-fresh-theming/")

