## ----boxAPI-code, eval=FALSE--------------------------------------------------
#  library(shiny)
#  library(bs4Dash)
#  ui <- dashboardPage(
#    title = "Box API",
#    dashboardHeader(),
#    dashboardSidebar(),
#    dashboardBody(
#      tags$style("body { background-color: ghostwhite}"),
#      fluidRow(
#        actionButton("toggle_box", "Toggle Box"),
#        actionButton("remove_box", "Remove Box", class = "bg-danger"),
#        actionButton("restore_box", "Restore Box", class = "bg-success"),
#        actionButton("update_box", "Update Box", class = "bg-primary")
#      ),
#      br(),
#      box(
#        title = textOutput("box_state"),
#        "Box body",
#        id = "mybox",
#        collapsible = TRUE,
#        closable = TRUE,
#        plotOutput("plot")
#      )
#    )
#  )
#  
#  server <- function(input, output, session) {
#    output$plot <- renderPlot({
#      req(!input$mybox$collapsed)
#      plot(rnorm(200))
#    })
#  
#    output$box_state <- renderText({
#      state <- if (input$mybox$collapsed) "collapsed" else "uncollapsed"
#      paste("My box is", state)
#    })
#  
#    observeEvent(input$toggle_box, {
#      updateBox("mybox", action = "toggle")
#    })
#  
#    observeEvent(input$remove_box, {
#      updateBox("mybox", action = "remove")
#    })
#  
#    observeEvent(input$restore_box, {
#      updateBox("mybox", action = "restore")
#    })
#  
#    observeEvent(input$update_box, {
#      updateBox(
#        "mybox",
#        action = "update",
#        options = list(
#          title = h2("New title", dashboardBadge(1, color = "primary")),
#          status = "danger",
#          solidHeader = TRUE,
#          width = 4
#        )
#      )
#    })
#  
#    observeEvent(input$mybox$visible, {
#      collapsed <- if (input$mybox$collapsed) "collapsed" else "uncollapsed"
#      visible <- if (input$mybox$visible) "visible" else "hidden"
#      message <- paste("My box is", collapsed, "and", visible)
#      toast(
#        title = message,
#        options = list(
#          autohide = TRUE,
#          class = "bg-pink",
#          position = "topRight"
#        )
#      )
#    })
#  }
#  
#  shinyApp(ui, server)

## ---- echo=FALSE--------------------------------------------------------------
bs4Dash:::app_container("https://dgranjon.shinyapps.io/bs4Dash-box-api/", deps = TRUE)

## ----boxTools, echo=FALSE, fig.cap='Box Tools. From left to right: boxLabel, boxDropdown, collapsible and closable buttons, boxSidebar trigger.', fig.align = 'center', out.width='50%'----
knitr::include_graphics("figures/boxTools.png")

## ----boxSidebar-code, eval=FALSE----------------------------------------------
#  shinyApp(
#    ui = dashboardPage(
#      header = dashboardHeader(),
#      body = dashboardBody(
#        box(
#          title = "Update box sidebar",
#          closable = TRUE,
#          width = 12,
#          height = "500px",
#          solidHeader = FALSE,
#          collapsible = TRUE,
#          actionButton("update", "Toggle card sidebar"),
#          sidebar = boxSidebar(
#            id = "mycardsidebar",
#            sliderInput(
#              "obs",
#              "Number of observations:",
#              min = 0,
#              max = 1000,
#              value = 500
#            )
#          ),
#          plotOutput("distPlot")
#        )
#      ),
#      sidebar = dashboardSidebar()
#    ),
#    server = function(input, output, session) {
#      observe(print(input$mycardsidebar))
#  
#      output$distPlot <- renderPlot({
#        hist(rnorm(input$obs))
#      })
#  
#      observeEvent(input$update, {
#        updateBoxSidebar("mycardsidebar")
#      })
#    }
#  )

## ---- echo=FALSE--------------------------------------------------------------
bs4Dash:::app_container("https://dgranjon.shinyapps.io/bs4Dash-box-sidebar/")

## ----boxDropdown-code, eval=FALSE---------------------------------------------
#  shinyApp(
#    ui = dashboardPage(
#      dashboardHeader(),
#      dashboardSidebar(),
#      dashboardBody(
#        box(
#          title = "Closable Box with dropdown",
#          closable = TRUE,
#          width = 12,
#          status = "warning",
#          solidHeader = FALSE,
#          collapsible = TRUE,
#          dropdownMenu = boxDropdown(
#            boxDropdownItem("Click me", id = "dropdownItem", icon = icon("heart")),
#            boxDropdownItem("item 2", href = "https://www.google.com/"),
#            dropdownDivider(),
#            boxDropdownItem("item 3", icon = icon("th"))
#          ),
#          "My box"
#        )
#      )
#    ),
#    server = function(input, output) {
#      observeEvent(input$dropdownItem, {
#        toast(
#          title = "Hi!",
#          options = list(
#            autohide = TRUE,
#            class = "bg-pink",
#            position = "topRight"
#          )
#        )
#      })
#    }
#  )

## ---- echo=FALSE--------------------------------------------------------------
bs4Dash:::app_container("https://dgranjon.shinyapps.io/bs4Dash-box-dropdown/")

## ----user-description, eval=FALSE---------------------------------------------
#  userDescription(
#    title = "Nadia Carmichael",
#    subtitle = "lead Developer",
#    type = 2,
#    image = "https://adminlte.io/themes/AdminLTE/dist/img/user7-128x128.jpg",
#  )

## ----userBox-code, eval=FALSE-------------------------------------------------
#  shinyApp(
#    ui = dashboardPage(
#      dashboardHeader(),
#      dashboardSidebar(),
#      dashboardBody(
#        actionButton("update_box", "Update"),
#        userBox(
#          id = "userbox",
#          title = userDescription(
#            title = "Nadia Carmichael",
#            subtitle = "lead Developer",
#            type = 2,
#            image = "https://adminlte.io/themes/AdminLTE/dist/img/user7-128x128.jpg",
#          ),
#          status = "primary",
#          gradient = TRUE,
#          background = "primary",
#          boxToolSize = "xl",
#          "Some text here!",
#          footer = "The footer here!"
#        )
#      ),
#      title = "userBox"
#    ),
#    server = function(input, output) {
#      observeEvent(input$update_box, {
#        updateBox(
#          "userbox",
#          action = "update",
#          options = list(
#            title = userDescription(
#              title = "Jean Box",
#              subtitle = "Developer",
#              type = 1,
#              image = "https://adminlte.io/themes/AdminLTE/dist/img/user3-128x128.jpg",
#            ),
#            status = "danger",
#            background = NULL,
#            solidHeader = FALSE,
#            width = 4
#          )
#        )
#      })
#    }
#  )

## ---- echo=FALSE--------------------------------------------------------------
bs4Dash:::app_container("https://dgranjon.shinyapps.io/bs4Dash-userBox/")

## ----user-block, eval=FALSE---------------------------------------------------
#  userBlock(
#    image = "https://adminlte.io/themes/AdminLTE/dist/img/user4-128x128.jpg",
#    title = "Social Box",
#    subtitle = "example-01.05.2018"
#  )

## ---- eval=FALSE--------------------------------------------------------------
#  shinyApp(
#    ui = dashboardPage(
#      dashboardHeader(),
#      dashboardSidebar(),
#      dashboardBody(
#        socialBox(
#          id = "socialbox",
#          title = userBlock(
#            image = "https://adminlte.io/themes/AdminLTE/dist/img/user4-128x128.jpg",
#            title = "Social Box",
#            subtitle = "example-01.05.2018"
#          ),
#          actionButton("update_box", "Refresh"),
#          "Some text here!",
#          br(), br(),
#          tabsetPanel(
#            tabPanel(
#              "News",
#              attachmentBlock(
#                image = "https://www.sammobile.com/wp-content/uploads/2017/11/Camel.png",
#                title = "Test",
#                href = "https://google.com",
#                "This is the content"
#              )
#            ),
#            tabPanel(
#              "Messages",
#              userMessages(
#                width = 12,
#                status = "danger",
#                userMessage(
#                  author = "Alexander Pierce",
#                  date = "20 Jan 2:00 pm",
#                  image = "https://adminlte.io/themes/AdminLTE/dist/img/user1-128x128.jpg",
#                  type = "received",
#                  "Is this template really for free? That's unbelievable!"
#                ),
#                userMessage(
#                  author = "Sarah Bullock",
#                  date = "23 Jan 2:05 pm",
#                  image = "https://adminlte.io/themes/AdminLTE/dist/img/user3-128x128.jpg",
#                  type = "sent",
#                  "You better believe it!"
#                )
#              )
#            )
#          ),
#          lapply(X = 1:10, FUN = function(i) {
#            boxComment(
#              image = "https://adminlte.io/themes/AdminLTE/dist/img/user3-128x128.jpg",
#              title = paste("Comment", i),
#              date = "01.05.2018",
#              paste0("The ", i, "-th comment")
#            )
#          }),
#          footer = "The footer here!"
#        )
#      ),
#      title = "Social Box"
#    ),
#    server = function(input, output) {
#      observeEvent(input$update_box, {
#        updateBox(
#          "socialbox",
#          action = "update",
#          options = list(
#            title = userBlock(
#              image = "https://adminlte.io/themes/AdminLTE/dist/img/user3-128x128.jpg",
#              title = "Social Box updated",
#              subtitle = "today"
#            )
#          )
#        )
#      })
#    }
#  )

## ---- echo=FALSE--------------------------------------------------------------
bs4Dash:::app_container("https://dgranjon.shinyapps.io/bs4Dash-socialBox/")

