library(splines)
library(lattice)


data("mtcars")

workingModel <- lm(mpg ~ factor(gear) + bs(wt, knots = 5) + hp, data = mtcars)

bs(mtcars$wt, knots = 4)

bs(cbind(mtcars$wt,mtcars$hp), knots = 4)

model <- lm(mpg ~ ns(wt, knots = 3) + ns(hp, knots = 3) - 1, mtcars)
model <- lm(mpg ~ ns(wt, knots = 2):ns(hp, knots = 2), mtcars)
model.matrix(model)
# grid over wt and hp

myGrid <- expand.grid(wt = seq(1.6, 5.3, length.out = 20),
                      hp = seq(52, 335, length.out = 20))

myGrid$mpgPred <- predict(model, newdata = myGrid)


wireframe(mpgPred ~ wt * hp, data=myGrid, shade=TRUE, scales=list(arrows=FALSE))


x<- unique(myGrid$wt)
y <- unique(myGrid$hp)
z <- matrix(myGrid$mpgPred, ncol = length(x), nrow = length(y))

plot_ly() %>% add_surface(x = ~x, y = ~y, z = ~z)

library("scatterplot3d") 

scatterplot3d(x = mtcars$wt, y=mtcars$hp, z=mtcars$mpg)

library(plotly)
fig <- plot_ly(mtcars, x = ~wt, y = ~hp, z = ~mpg)
fig <- fig %>% add_markers()

fig


nim_approx_tri(myModel$xC, myModel$xH, myModel$xD,
               myModel$grid, 
               myModel$yAlarm,
               myModel$smoothC[10], myModel$smoothH[10], myModel$smoothD[10])
