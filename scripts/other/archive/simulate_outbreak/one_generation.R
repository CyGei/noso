source(here::here("scripts", "helpers.R"))

gamma <- c(15, 0.42)
f <- c(0.2, 0.8)
pi <- pi_formula(gamma, f)
R0 <- c(100, 100)
group <- c("A", "B")

true_table <- matrix(c(31, 7, 13, 11), nrow = 2, byrow = TRUE,
                     dimnames = list(c("A", "B"), c("A", "B"))) %>%
  prop.table() %>%
  round(2)



# Simulate one generation
#Group A
idx <- which(group == "A")
A_A = pi[idx] * R0[idx]
A_B = R0[idx] - A_A

#Group B
idx <- which(group == "B")
B_B = pi[idx] * R0[idx]
B_A = R0[idx] - B_B

#table
table_start <- matrix(c(A_A, A_B, B_A, B_B), nrow = 2, byrow = TRUE,
                      dimnames = list(c("A", "B"), c("A", "B"))) %>%
  prop.table() %>%
  round(2)

#Mo
obs = diag(table_start) / colSums(table_start)
exp = colSums(table_start) / sum(table_start)
obs / exp


# Relative Sucecptibility
S <- c(2, 1)
idx <- which(group == "A")
A_A = (pi[idx] * R0[idx]) * S[idx]
idx <- which(group == "B")
B_B = (pi[idx] * R0[idx]) * S[idx]

#table
table <- matrix(c(A_A, A_B, B_A, B_B), nrow = 2, byrow = TRUE,
                dimnames = list(c("A", "B"), c("A", "B"))) %>%
  prop.table() %>%
  round(2)

newR0 <- rowSums(table)

# Renormalize R0
normR0 <- R0 * (R0/newR0)
normR0 <- round(normR0)

#Final model
A_A = pi[1] * normR0[1]
A_B = normR0[1] - A_A
A_A = A_A * S[1]
B_B = pi[2] * normR0[2]
B_A = normR0[2] - B_B
B_B = B_B * S[2]

#table
table <- matrix(c(A_A, A_B, B_A, B_B), nrow = 2, byrow = TRUE,
                dimnames = list(c("A", "B"), c("A", "B"))) %>%
  prop.table() %>%
  round(2)
#Mo
obs = diag(table) / colSums(table)
exp = colSums(table) / sum(table)
obs / exp


# tets --------------------------------------------------------------------


gamma <- c(30, 0.42)
f <- c(0.2, 0.8)
pi <- get_pi(gamma, f)
R0 <- c(100, 100)
group <- c("A", "B")

# Simulate one generation
#Group A
idx <- which(group == "A")

A_A = pi[idx] * R0[idx]
A_B = R0[idx] - A_A

#Group B
idx <- which(group == "B")
B_B = pi[idx] * R0[idx]
B_A = R0[idx] - B_B

#table
table_start <- matrix(c(A_A, A_B, B_A, B_B), nrow = 2, byrow = TRUE)
colnames(table_start) <- c("A", "B")
rownames(table_start) <- c("A", "B")
round(table_start)
round(prop.table(table_start), 2)
