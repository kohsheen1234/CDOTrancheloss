
library(shiny)

lossdistr <- function(x, defthresh=-2, rho=0.1, lgd=0.4) {
  qnormv <- ifelse(x/lgd < 0.999, qnorm(x/lgd), 3.1)
  sqrt((1-rho)/rho)*exp(-(sqrt(1-rho)*qnormv - defthresh)^2/(2*rho) + qnormv^2/2)/lgd
}  # end lossdistr

# Define Vasicek cumulative loss distribution
# (with error handling for x)
cumlossdistr <- function(x, defthresh=(-2), rho=0.2, lgd=0.4) {
  qnormv <- ifelse(x/lgd < 0.999, qnorm(x/lgd), 3.1)
  pnorm((sqrt(1-rho)*qnormv - defthresh)/sqrt(rho))
} # end cumlossdistr


# Define UI for shiny app
ui <- fluidPage(
  titlePanel("CDO Tranche Expected Losses"),
  fluidRow(
    column(width=4, sliderInput("rho", label="Correlation:",
                                min=0.0, max=0.9, value=0.1, step=0.01)),
    column(width=4,sliderInput("defprob", label = "Default probability:",
                               min = 0, max = 1, value = 0.2, step = 0.01)),
    column(width=4, sliderInput("lgd", label="Loss severity:",
                                min=0.0, max=0.9, value=0.4, step=0.01)),
    column(width=4, sliderInput("attachp", label="CDO attachment point:",
                                min=0.0, max=0.5, value=0.15, step=0.01)),
    column(width=4,sliderInput("detachp", label = "CDO detachment point:",
                               min = 0, max = 1, value = 0.2, step = 0.01)),
  ),
  shiny::plotOutput("loss_plot")
)

server <- function(input, output) {
  
  # Calculate expected tranche loss
  tranchel <- reactive({
    attachp <- input$attachp
    defprob <- input$defprob
    rho <- input$rho
    lgd <- input$lgd
    exploss <- lgd*defprob
    defthresh <- qnorm(defprob)
    detachp <- input$detachp
    
    # Integrate the loss distribution between the attachment and detachment points
    tranche_loss <- round(integrate(function(x, attachp) (x-attachp)*lossdistr(x, defthresh=defthresh, rho=rho, lgd=lgd),
                                    low=attachp, up=detachp, attachp=attachp)$value/(detachp-attachp) + # Loss in excess of detachp
                            (1-cumlossdistr(x=detachp, defthresh=defthresh, rho=rho, lgd=lgd)),digits = 5)
    
  })
  
  output$loss_plot <- renderPlot({
    # Generate loss values to plot
    loss_vals <- seq(0, 3*input$lgd*input$defprob, length.out=1000)
    loss_densities <- sapply(loss_vals, function(x) lossdistr(x, defthresh=qnorm(input$defprob), rho=input$rho, lgd=input$lgd))
    
    attachp <- input$attachp
    defprob <- input$defprob
    rho <- input$rho
    lgd <- input$lgd
    exploss <- lgd*defprob
    defthresh <- qnorm(defprob)
    detachp <- input$detachp
    
    exploss <- lgd*defprob
    defthresh <- qnorm(defprob)
    
    # Calculate max x-axis range
    xmax <- max(3*exploss, detachp)
    ymax <- max(sapply(seq(fr=0.01, to=lgd/2, length.out=10), lossdistr, defthresh=defthresh, rho=rho, lgd=lgd))
    
    
    curve(expr=lossdistr(x, defthresh=defthresh, rho=rho, lgd=lgd),
          cex.main=1.5, cex.lab=1.5, cex.axis=1.5, 
          type="l", xlim=c(0, xmax), 
          xlab="Percentage loss", ylab="Density", lwd=3,
          col="orange", main="Distribution of Losses")
    
    # Add vertical lines
    abline(v = input$attachp, lty = 2, col = "blue", lwd = 2, 
           name = "CDO Attachment Point")
    text(x=input$attachp-0.001, y=ymax/2, labels="CDO Attachment Point",
         lwd=2, srt=90, pos=3, cex=1)
    abline(v = input$detachp, lty = 2, col = "blue", lwd = 2, 
           name = "CDO Detachment Point")
    text(x=input$detachp-0.001, y=ymax/2, labels="CDO Detachment Point",
         lwd=2, srt=90, pos=3, cex=1)
    abline(v = input$lgd*input$defprob, lty = 2, col = "red", lwd = 2, 
           name = "Expected Loss")
    text(x=input$lgd*input$defprob-0.001, y=ymax/2, labels="expected loss",
         lwd=2, srt=90, pos=3, cex=1)
    
    
    # Calculate tranche shading for CVaR
    var_max <- detachp
    varv <- seq(attachp, var_max, length=100)
    densv <- sapply(varv, lossdistr, defthresh=defthresh, rho=rho, lgd=lgd)
    # Draw shaded polygon
    polygon(c(attachp, varv, var_max),
            c(-1, densv, -1), col="red", border=NA, density=10)
    # text(x=0.045, y=0, labels="CVaR", lwd=2, pos=3)
    
    # Add text
    text(input$lgd*input$defprob-0.001, max(loss_densities), 
         "expected loss", pos = 4, size = 6, col = "red")
    
    # Text with tranche attachment
    text(xmax, ymax, 
         lab=paste0(
           "Default probability = ", format(100*input$defprob, digits=3), "%", "\n",
           "Loss severity = ", format(100*input$lgd, digits=3), "%", "\n",
           "Correlation = ", format(100*input$rho, digits=3), "%", "\n",
           "Tranche attachment = ", format(100*input$attachp, digits=3), "%", "\n",
           "Tranche detachment = ", format(100*input$detachp, digits=3), "%", "\n",
           "Expected Loss = ", format(100*tranchel(), digits=5), "%", "\n"),
         adj=c(1, 1), cex=1.5, lwd=2)
  })
}

# Create Shiny object
shinyApp(ui = ui, server = server)

