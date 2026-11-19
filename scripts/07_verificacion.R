# ============================================================================
# 07_verificacion.R — Verificación estadística de soporte a la sustentación
# ----------------------------------------------------------------------------
# Responde, con cifras reproducibles, las preguntas típicas de una defensa:
#   A. Factores de expansión y estimador del total
#   B. Modelo de regresión: coeficientes, colinealidad y aporte de variables
#   C. Validación de hiperparámetros con validación cruzada espacial
#   D. Comparación con modelos simples (árbol de decisión y regresión logística)
#   E. Análisis de patrones de puntos (muestra y activos Climate TRACE)
#   F. Supuestos del modelo de diseño de experimentos
# Uso: Rscript scripts/07_verificacion.R   (no modifica ningún dato)
# ============================================================================
options(width = 220)
p <- function(...) cat(..., "\n", sep = "")
suppressPackageStartupMessages({
  library(dplyr); library(sf); library(ranger); library(rpart)
  library(pROC); library(spdep); library(car); library(FNN)
})
sf_use_s2(FALSE)
g <- readRDS("datos/procesados/grilla_co2.rds")
m <- readRDS("datos/procesados/muestra_final.rds")

# ---------------------------------------------------------------------------
# A. FACTORES DE EXPANSIÓN
# ---------------------------------------------------------------------------
p("=========== A. FACTORES DE EXPANSION ===========")
agg <- m |> st_drop_geometry() |> group_by(estrato_final) |>
  summarise(n = n(), N_h = first(N_h), n_h = first(n_h), peso = first(peso),
            y_medio = mean(co2_2023), aporte_mt = sum(co2_2023 * peso) / 1e6,
            .groups = "drop")
print(as.data.frame(agg))
p("HT total (Mt): ", round(sum(m$co2_2023 * m$peso) / 1e6, 2))
p("suma de pesos: ", round(sum(m$peso)), " | N del marco: ", nrow(g))
p("peso: min ", round(min(m$peso), 2), " | max ", round(max(m$peso), 2),
  " | medio ", round(mean(m$peso), 2))

# ---------------------------------------------------------------------------
# B. MODELO DE REGRESIÓN
# ---------------------------------------------------------------------------
p("")
p("=========== B. MODELO DE REGRESION (OLS) ===========")
a <- g |> st_drop_geometry() |>
  mutate(log_dens = log_dens_2023, log_pob = log1p(poblacion_2020),
         log_popdens = log1p(pop_dens_2020), log_dist = log1p(dist_fncer_km),
         region = factor(region), urbano = factor(urbano)) |>
  select(celda_edgar, x_9377, y_9377, log_dens, log_pob, log_popdens, urbano,
         region, log_dist, n_proy, cap_mw, area_km2, lat, lon) |> na.omit()
fit <- lm(log_dens ~ log_pob + log_popdens + urbano + region + log_dist +
            n_proy + cap_mw + area_km2 + lat + lon, data = a)
s <- summary(fit)
print(round(s$coefficients, 5))
p("R2 en muestra: ", round(s$r.squared, 4), " | R2 ajustado: ",
  round(s$adj.r.squared, 4), " | n: ", nrow(a))
p("--- VIF (colinealidad) ---")
print(round(tryCatch(car::vif(fit), error = function(e) NA), 2))
num <- c("log_pob", "log_popdens", "log_dist", "n_proy", "cap_mw", "area_km2", "lat", "lon")
b <- coef(fit)[num]; sx <- sapply(a[num], sd); sy <- sd(a$log_dens)
p("--- beta estandarizados ---"); print(round(b * sx / sy, 3))
av <- anova(fit); ss <- av[["Sum Sq"]]
p("--- eta2 secuencial (tipo I) ---")
print(round(ss[-length(ss)] / (ss[-length(ss)] + ss[length(ss)]), 3))
p("--- correlacion entre predictores ---")
print(round(cor(a[, c("log_pob", "log_popdens", "log_dist", "area_km2", "lat", "lon")]), 2))

# ---------------------------------------------------------------------------
# C. HIPERPARÁMETROS
# ---------------------------------------------------------------------------
p("")
p("=========== C. BUSQUEDA DE HIPERPARAMETROS (CV espacial, 5 bloques) ===========")
X <- model.matrix(~ log_pob + log_popdens + urbano + region + log_dist + n_proy +
                    cap_mw + area_km2 + lat + lon, data = a)[, -1]
y <- a$log_dens
d_reg <- data.frame(y = y, X, check.names = FALSE)
d_x <- data.frame(X, check.names = FALSE)
set.seed(2026)
kmb <- kmeans(a[, c("x_9377", "y_9377")], centers = 5, nstart = 20)$cluster
rmse <- function(o, pr) sqrt(mean((o - pr)^2))
r2 <- function(o, pr) 1 - sum((o - pr)^2) / sum((o - mean(o))^2)
grid <- expand.grid(mtry = c(2, 3, 4, 6, 8), mns = c(1, 5, 20))
cv1 <- function(mtry, mns) {
  pr <- numeric(length(y))
  for (k in 1:5) {
    tr <- kmb != k; te <- kmb == k
    f <- ranger(x = X[tr, , drop = FALSE], y = y[tr], num.trees = 300, mtry = mtry,
                min.node.size = mns, num.threads = 4, seed = 2026)
    pr[te] <- predict(f, data = X[te, , drop = FALSE])$predictions
  }
  c(RMSE = rmse(y, pr), R2 = r2(y, pr))
}
res <- cbind(grid, round(t(mapply(cv1, grid$mtry, grid$mns)), 4))
print(res[order(res$RMSE), ], row.names = FALSE)
p("Mejor: mtry = ", res$mtry[which.min(res$RMSE)], ", min.node.size = ",
  res$mns[which.min(res$RMSE)], " -> RMSE = ", round(min(res$RMSE), 4))

# ---------------------------------------------------------------------------
# D. MODELOS SIMPLES
# ---------------------------------------------------------------------------
p("")
p("=========== D1. REGRESION: OLS vs ARBOL vs RF ===========")
pr_ols <- pr_tree <- pr_rf <- numeric(length(y))
for (k in 1:5) {
  tr <- kmb != k; te <- kmb == k
  pr_ols[te] <- predict(lm(y ~ ., data = d_reg[tr, ]), d_reg[te, ])
  pr_tree[te] <- predict(rpart(y ~ ., data = d_reg[tr, ], method = "anova", cp = 0.01),
                         d_reg[te, ])
  fr <- ranger(x = X[tr, , drop = FALSE], y = y[tr], num.trees = 300, mtry = 4,
               num.threads = 4, seed = 2026)
  pr_rf[te] <- predict(fr, data = X[te, , drop = FALSE])$predictions
}
print(data.frame(modelo = c("OLS", "Arbol (rpart)", "Random Forest"),
                 RMSE = round(c(rmse(y, pr_ols), rmse(y, pr_tree), rmse(y, pr_rf)), 4),
                 R2_agrupado = round(c(r2(y, pr_ols), r2(y, pr_tree), r2(y, pr_rf)), 4)),
      row.names = FALSE)

p("")
p("=========== D2. CLASIFICACION: alta emision (sobre la mediana) ===========")
alta <- as.integer(a$log_dens > median(a$log_dens))
pr_l <- pr_r <- pr_t <- numeric(nrow(a))
for (k in 1:5) {
  tr <- kmb != k; te <- kmb == k
  gl <- glm(alta ~ ., data = data.frame(alta = alta[tr], d_x[tr, ]), family = binomial)
  pr_l[te] <- predict(gl, newdata = data.frame(d_x[te, ]), type = "response")
  fr <- ranger(x = X[tr, , drop = FALSE], y = factor(alta[tr]), num.trees = 300,
               probability = TRUE, num.threads = 4, seed = 2026)
  pr_r[te] <- predict(fr, data = X[te, , drop = FALSE])$predictions[, 2]
  tt <- rpart(factor(alta) ~ ., data = data.frame(alta = factor(alta[tr]), d_x[tr, ]),
              cp = 0.01)
  pr_t[te] <- predict(tt, newdata = data.frame(d_x[te, ]), type = "prob")[, 2]
}
auc <- function(pr) as.numeric(pROC::auc(pROC::roc(alta, pr, quiet = TRUE)))
p("AUC logistica: ", round(auc(pr_l), 4), " | AUC arbol: ", round(auc(pr_t), 4),
  " | AUC Random Forest: ", round(auc(pr_r), 4))
for (nm in c("logistica", "arbol", "RF")) {
  pr <- switch(nm, logistica = pr_l, arbol = pr_t, RF = pr_r)
  p("  exactitud ", nm, ": ",
    round(100 * mean(as.integer(pr > 0.5) == alta), 1), " %")
}
gl <- glm(alta ~ ., data = data.frame(alta = alta, d_x), family = binomial)
p("--- coeficientes logisticos (log-odds, p) ---")
print(round(summary(gl)$coefficients[, c(1, 4)], 4))
p("--- arbol de regresion (estructura) ---")
print(rpart(y ~ ., data = d_reg, method = "anova", cp = 0.01))

# ---------------------------------------------------------------------------
# E. PATRONES DE PUNTOS
# ---------------------------------------------------------------------------
p("")
p("=========== E. PATRONES DE PUNTOS ===========")
mc <- st_centroid(st_geometry(st_transform(m, 9377)))
xy_m <- st_coordinates(mc)
D <- as.matrix(dist(xy_m)) / 1000; diag(D) <- NA
nn <- apply(D, 1, min, na.rm = TRUE)
area <- 1137740
p("Muestra: n = ", length(nn), " | vecino mediano = ", round(median(nn), 2),
  " km | media = ", round(mean(nn), 2), " | p90 = ", round(quantile(nn, .9), 2),
  " | max = ", round(max(nn), 2))
p("  esperado aleatorio = ", round(0.5 / sqrt(400 / area), 2), " km | R Clark-Evans = ",
  round(mean(nn) / (0.5 / sqrt(400 / area)), 3))
p("  departamentos cubiertos: ", length(unique(m$depto)), " de 33 | municipios: ",
  length(unique(m$mun_id)), " de ", length(unique(g$mun_id)))
d2 <- as.matrix(st_distance(st_centroid(st_geometry(st_transform(g, 9377))), mc)) / 1000
dmin <- apply(d2, 1, min)
for (x in c(25, 50, 100)) p("  marco a <= ", x, " km de una celda muestreada: ",
                            round(100 * mean(dmin <= x), 1), " %")

ct <- read.csv("datos/climate_trace_COL_2023.csv", stringsAsFactors = FALSE)
ct <- ct[!is.na(ct$lon) & !is.na(ct$lat), ]
xy_todos <- as.matrix(ct[, c("lon", "lat")]); xy <- unique(xy_todos)
p("Climate TRACE: filas = ", nrow(xy_todos), " | ubicaciones unicas = ", nrow(xy))
dnn <- get.knnx(xy, xy, k = 2)$nn.dist[, 2] * 111
lam <- nrow(xy) / area
p("  vecino mediano = ", round(median(dnn), 2), " km | esperado = ",
  round(0.5 / sqrt(lam), 2), " km | R Clark-Evans = ", round(mean(dnn) / (0.5 / sqrt(lam)), 3))
p("  % de ubicaciones a <5 km de otra: ", round(100 * mean(dnn < 5), 1))
p("  vecino mediano por sector (km):")
print(round(sort(tapply(dnn, ct$sector[!duplicated(xy_todos)], median)), 2))

# ---------------------------------------------------------------------------
# F. SUPUESTOS DEL DOE
# ---------------------------------------------------------------------------
p("")
p("=========== F. SUPUESTOS DEL DOE ===========")
do <- readRDS("datos/procesados/resultados_doe.rds")
print(do$diagnostico)
p("--- robustez (resumen) ---")
print(do$robustez$tabla[, c("prueba", "efecto", "estadistico", "gl", "p", "eta2_parc")],
      row.names = FALSE)
p("¿Se mantienen las conclusiones? ", do$robustez$mantiene_conclusiones)
