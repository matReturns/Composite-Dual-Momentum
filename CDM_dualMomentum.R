



# ============================================================
# Composite Dual Momentum Replication
# Gary Antonacci-style module-based dual momentum
# ============================================================

stratStats <- function(rets, digits = 4) {
  
  require(PerformanceAnalytics)
  
  rets <- na.omit(rets)
  
  stats <- rbind(
    "Annualized Return" = Return.annualized(rets),
    "Annualized Std Dev" = StdDev.annualized(rets),
    "Annualized Sharpe (Rf=0%)" = SharpeRatio.annualized(rets, Rf = 0),
    "Worst Drawdown" = maxDrawdown(rets),
    "Calmar Ratio" = CalmarRatio(rets)
  )
  
  return(round(stats, digits))
}



run_composite_dual_momentum <- function(dataStartDate = "2006-01-01",
                                        analysisStartDate = "2008-07-01",
                                        endDate = Sys.Date(),
                                        modules = list(
                                          Equities = c("SPY", "EFA"),
                                          Credit = c("HYG", "BIV"),
                                          RealEstate = c("VNQ", "REM"),
                                          EconomicStress = c("GLD", "TLT")
                                        ),
                                        proxyAsset = "BIL",
                                        lookbackMonths = 12,
                                        verbose = TRUE) {
  
  require(quantmod)
  require(PerformanceAnalytics)
  require(xts)
  require(zoo)
  
  # -----------------------------
  # Check module structure
  # -----------------------------
  if (!all(sapply(modules, length) == 2)) {
    stop("Each module must contain exactly two assets.")
  }
  
  riskyAssets <- unique(unlist(modules))
  allSymbols <- unique(c(riskyAssets, proxyAsset))
  nModules <- length(modules)
  
  # -----------------------------
  # Download adjusted prices
  # -----------------------------
  priceList <- list()
  
  for (sym in allSymbols) {
    
    if (verbose) {
      message("Downloading: ", sym)
    }
    
    tmp <- getSymbols(
      sym,
      from = dataStartDate,
      to = endDate,
      auto.assign = FALSE,
      warnings = FALSE
    )
    
    px <- Ad(tmp)
    colnames(px) <- sym
    priceList[[sym]] <- px
  }
  
  dailyPrices <- do.call(merge, priceList)
  dailyPrices <- na.omit(dailyPrices[, allSymbols])
  
  # -----------------------------
  # Convert to month-end prices and returns
  # -----------------------------
  monthlyPrices <- apply.monthly(dailyPrices, last)
  monthlyPrices <- na.omit(monthlyPrices)
  
  monthlyReturns <- na.omit(Return.calculate(monthlyPrices))
  monthlyReturns <- monthlyReturns[, allSymbols]
  
  if (verbose) {
    message("Monthly return data begins: ", as.character(first(index(monthlyReturns))))
    message("Monthly return data ends:   ", as.character(last(index(monthlyReturns))))
  }
  
  # -----------------------------
  # Helper: cumulative return
  # -----------------------------
  cumulative_return <- function(x) {
    prod(1 + as.numeric(x), na.rm = FALSE) - 1
  }
  
  # -----------------------------
  # 12-month momentum
  # -----------------------------
  momentum <- rollapply(
    monthlyReturns,
    width = lookbackMonths,
    FUN = cumulative_return,
    by.column = TRUE,
    align = "right",
    fill = NA
  )
  
  momentum <- momentum[, allSymbols]
  
  # -----------------------------
  # Preallocate weights and logs
  # -----------------------------
  nRows <- NROW(monthlyReturns)
  nCols <- length(allSymbols)
  
  weightsMat <- matrix(
    NA_real_,
    nrow = nRows,
    ncol = nCols
  )
  
  colnames(weightsMat) <- allSymbols
  
  choiceRecords <- vector("list", nRows * nModules)
  recordCounter <- 1
  
  # -----------------------------
  # Main signal loop
  # -----------------------------
  for (i in seq_len(nRows)) {
    
    currentMomentum <- momentum[i, ]
    
    # Skip months before the 12-month signal is valid
    if (any(is.na(currentMomentum))) {
      next
    }
    
    w <- rep(0, nCols)
    names(w) <- allSymbols
    
    for (moduleName in names(modules)) {
      
      pair <- modules[[moduleName]]
      assetA <- pair[1]
      assetB <- pair[2]
      
      momA <- as.numeric(currentMomentum[, assetA])
      momB <- as.numeric(currentMomentum[, assetB])
      proxyMom <- as.numeric(currentMomentum[, proxyAsset])
      
      # Relative momentum: choose stronger asset inside the module
      if (momA >= momB) {
        relativeWinner <- assetA
        winnerMomentum <- momA
      } else {
        relativeWinner <- assetB
        winnerMomentum <- momB
      }
      
      # Absolute momentum: compare winner to proxy asset
      if (winnerMomentum > proxyMom) {
        finalAsset <- relativeWinner
        usesProxy <- FALSE
      } else {
        finalAsset <- proxyAsset
        usesProxy <- TRUE
      }
      
      # Each module gets 25% of portfolio
      w[finalAsset] <- w[finalAsset] + 1 / nModules
      
      choiceRecords[[recordCounter]] <- data.frame(
        date = as.Date(index(monthlyReturns)[i]),
        module = moduleName,
        assetA = assetA,
        assetB = assetB,
        relativeWinner = relativeWinner,
        finalAsset = finalAsset,
        winnerMomentum = winnerMomentum,
        proxyMomentum = proxyMom,
        usesProxy = usesProxy,
        stringsAsFactors = FALSE
      )
      
      recordCounter <- recordCounter + 1
    }
    
    weightsMat[i, ] <- w
  }
  
  # -----------------------------
  # Convert weights to xts and lag by one month
  # -----------------------------
  weights <- xts(
    weightsMat,
    order.by = index(monthlyReturns)
  )
  
  # Important: lag weights to avoid look-ahead bias.
  weightsLag <- lag(weights, k = 1)
  
  # -----------------------------
  # Strategy returns
  # -----------------------------
  cdmReturns <- xts(
    rowSums(weightsLag * monthlyReturns[, allSymbols], na.rm = FALSE),
    order.by = index(monthlyReturns)
  )
  
  colnames(cdmReturns) <- "CDM"
  
  # -----------------------------
  # Equal-weight benchmark
  # Monthly rebalanced equal weight across the same eight risky assets
  # -----------------------------
  equalWeightReturns <- xts(
    rowMeans(monthlyReturns[, riskyAssets], na.rm = FALSE),
    order.by = index(monthlyReturns)
  )
  
  colnames(equalWeightReturns) <- "EqualWeight"
  
  # -----------------------------
  # Optional proxy return reference
  # -----------------------------
  proxyReturns <- monthlyReturns[, proxyAsset]
  colnames(proxyReturns) <- paste0("Proxy_", proxyAsset)
  
  # -----------------------------
  # Combine and apply analysis start date
  # -----------------------------
  returns <- merge(
    cdmReturns,
    equalWeightReturns
  )
  
  returns <- returns[paste0(analysisStartDate, "/")]
  returns <- na.omit(returns)
  
  referenceReturns <- merge(
    cdmReturns,
    equalWeightReturns,
    proxyReturns
  )
  
  referenceReturns <- referenceReturns[paste0(analysisStartDate, "/")]
  referenceReturns <- na.omit(referenceReturns)
  
  # -----------------------------
  # Clean choice log
  # -----------------------------
  choiceRecords <- choiceRecords[!sapply(choiceRecords, is.null)]
  choiceLog <- do.call(rbind, choiceRecords)
  
  choiceLog <- subset(
    choiceLog,
    date >= as.Date(first(index(returns))) &
      date <= as.Date(last(index(returns)))
  )
  
  # -----------------------------
  # Annual returns
  # -----------------------------
  annualReturns <- apply.yearly(returns, Return.cumulative)
  
  # -----------------------------
  # Exposure stats
  # Based on lagged live portfolio weights
  # -----------------------------
  liveWeights <- weightsLag[index(returns), ]
  proxyWeight <- liveWeights[, proxyAsset]
  
  exposureStats <- data.frame(
    avgProxyWeight = mean(as.numeric(proxyWeight), na.rm = TRUE),
    avgRiskyWeight = 1 - mean(as.numeric(proxyWeight), na.rm = TRUE),
    pctNoProxy = mean(as.numeric(proxyWeight) == 0, na.rm = TRUE),
    pctSomeProxy = mean(as.numeric(proxyWeight) > 0 & as.numeric(proxyWeight) < 1, na.rm = TRUE),
    pctAllProxy = mean(as.numeric(proxyWeight) == 1, na.rm = TRUE)
  )
  
  exposureStats <- round(exposureStats, 4)
  
  # -----------------------------
  # Module-level selection stats
  # -----------------------------
  moduleChoiceStats <- list()
  
  for (moduleName in names(modules)) {
    
    tmp <- subset(choiceLog, module == moduleName)
    
    moduleChoiceStats[[moduleName]] <- list(
      finalAssetCounts = table(tmp$finalAsset),
      finalAssetWeights = round(prop.table(table(tmp$finalAsset)), 4),
      proxyUseRate = round(mean(tmp$usesProxy), 4)
    )
  }
  
  # -----------------------------
  # Overall asset selection stats
  # -----------------------------
  finalAssetCounts <- table(choiceLog$finalAsset)
  finalAssetWeights <- round(prop.table(finalAssetCounts), 4)
  
  # -----------------------------
  # Return result object
  # -----------------------------
  return(list(
    returns = returns,
    referenceReturns = referenceReturns,
    summary = stratStats(returns),
    annualReturns = round(annualReturns, 4),
    exposureStats = exposureStats,
    monthlyPrices = monthlyPrices,
    monthlyReturns = monthlyReturns,
    momentum = momentum,
    weights = weights,
    weightsLag = weightsLag,
    choiceLog = choiceLog,
    moduleChoiceStats = moduleChoiceStats,
    finalAssetCounts = finalAssetCounts,
    finalAssetWeights = finalAssetWeights,
    settings = list(
      dataStartDate = dataStartDate,
      analysisStartDate = as.character(first(index(returns))),
      endDate = as.character(last(index(returns))),
      modules = modules,
      proxyAsset = proxyAsset,
      lookbackMonths = lookbackMonths,
      benchmark = "Monthly rebalanced equal weight portfolio using the same eight risky assets"
    )
  ))
}



cdmTest <- run_composite_dual_momentum(
  dataStartDate = "2006-01-01",
  analysisStartDate = "2007-07-01",
  endDate = "2026-07-09",
  modules = list(
    Equities = c("SPY", "EFA"),
    Credit = c("HYG", "BIV"),
    RealEstate = c("VNQ", "REM"),
    EconomicStress = c("GLD", "TLT")
  ),
  proxyAsset = "BIL",
  lookbackMonths = 12,
  verbose = TRUE
)




cdmTest$summary
cdmTest$exposureStats
cdmTest$annualReturns
cdmTest$finalAssetCounts
cdmTest$finalAssetWeights
cdmTest$moduleChoiceStats
cdmTest$settings




charts.PerformanceSummary(
  cdmTest$returns,
  main = "Composite Dual Momentum vs. Equal Weight Benchmark",
  wealth.index = T,
  colorset = c("darkgreen", "darkorange")
)



chart.CumReturns(
  cdmTest$returns,
  wealth.index = TRUE,
  main = "Composite Dual Momentum vs. Equal Weight Benchmark"
)



chart.Drawdown(
  cdmTest$returns,
  main = "Composite Dual Momentum Drawdowns vs. Equal Weight Benchmark",
  colorset = c("darkgreen", "darkorange"),
  legend.loc = "bottomright"
)



table.CalendarReturns(cdmTest$returns)



head(cdmTest$choiceLog, 20)
tail(cdmTest$choiceLog, 20)



cdmTest$moduleChoiceStats$Equities
cdmTest$moduleChoiceStats$Credit
cdmTest$moduleChoiceStats$RealEstate
cdmTest$moduleChoiceStats$EconomicStress



cdmTest$finalAssetCounts
cdmTest$finalAssetWeights















