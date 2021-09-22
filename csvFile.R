
# Module UI function
csvFileInput <- function(id, label = "CSV file") {
  # Create a namespace function using the provided id
  ns <- NS(id)

  tagList(
    fileInput(ns("file"), label),
    checkboxInput(ns("heading"), "Heading in first row"),
    selectInput(ns("sep"), "Values separated by", c(
      "None" = "",
      "Comma" = ",",
      "Semicolon" = ";",
      "tabs"="\t"
    ))
  )
}


csvFile <- function(input, output, session, stringsAsFactors) {
  # The selected file, if any
  userFile <- reactive({
    # If no file is selected, don't do anything
    validate(need(input$file, message = FALSE))
    input$file
  })

  # The user's data, parsed into a data frame
  dataframe <- reactive({
    read.csv(userFile()$datapath,
             header = input$heading,
             stringsAsFactors = stringsAsFactors,
             sep=input$sep)
  })

  # We can run observers in here if we want to
  observe({
    msg <- sprintf("File %s was uploaded", userFile()$name)
    cat(msg, "\n")
  })

  # Return the reactive that yields the data frame
  return(dataframe)
}
