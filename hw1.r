library(pryr)
library(stats4)

enclosing_env("glm")
is_active_binding(add.fit)
mem_change(enclosing_env("glm"))
names_c()
print
parenvs(abilitydata)

y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
nLL <- function(lambda) -sum(dpois(y, lambda, log = TRUE))
fit <- mle(nLL, start = list(lambda = 5), nobs = length(y))

otype(ability.lgt)
otype(res.dev)

loadedNamespaces()
parents <- function(env = globalenv()) {
  Name = paste0(if(isNamespace(env)) "namespace:" else "",
                environmentName(env))
  if(identical(env, emptyenv()))
    Name
  else
    c(Name, parents(parent.env(env)))
}

ftype(parents)
ftype(glm)
ftype(hist)

otype(parents)
otype(hist)

find_uses("package:base", "summary")
is_s3_generic("glm")
is_s3_method("glm")

object_size(add.fit) #add.fit is an additive glm model that was for STAT526
