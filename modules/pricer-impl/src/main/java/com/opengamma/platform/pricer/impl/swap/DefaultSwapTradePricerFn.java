/**
 * Copyright (C) 2014 - present by OpenGamma Inc. and the OpenGamma group of companies
 *
 * Please see distribution for license.
 */
package com.opengamma.platform.pricer.impl.swap;

import com.opengamma.basics.currency.MultiCurrencyAmount;
import com.opengamma.collect.ArgChecker;
import com.opengamma.platform.finance.swap.Swap;
import com.opengamma.platform.finance.swap.SwapProduct;
import com.opengamma.platform.finance.swap.SwapTrade;
import com.opengamma.platform.pricer.PricingEnvironment;
import com.opengamma.platform.pricer.TradePricerFn;
import com.opengamma.platform.pricer.swap.SwapProductPricerFn;

/**
 * Pricer implementation for swap trades.
 * <p>
 * The swap trade is priced by examining the embedded swap.
 */
public class DefaultSwapTradePricerFn
    implements TradePricerFn<SwapTrade> {

  /**
   * Default implementation.
   */
  public static final DefaultSwapTradePricerFn DEFAULT = new DefaultSwapTradePricerFn(
      ExpandingSwapProductPricerFn.DEFAULT);

  /**
   * Pricer for {@link SwapProduct}.
   */
  private final SwapProductPricerFn<? super Swap> swapPricerFn;

  /**
   * Creates an instance.
   * 
   * @param swapPricerFn  the swap pricer
   */
  public DefaultSwapTradePricerFn(
      SwapProductPricerFn<? super Swap> swapPricerFn) {
    this.swapPricerFn = ArgChecker.notNull(swapPricerFn, "swapPricerFn");
  }

  //-------------------------------------------------------------------------
  @Override
  public MultiCurrencyAmount presentValue(PricingEnvironment env, SwapTrade trade) {
    return swapPricerFn.presentValue(env, trade.getSwap());
  }

}