/**
 * Copyright (C) 2015 - present by OpenGamma Inc. and the OpenGamma group of companies
 * 
 * Please see distribution for license.
 */
package com.opengamma.strata.examples.report;

import com.opengamma.strata.product.TradeList;

/**
 * Parameter converter for {@link TradeList}.
 */
public class TradeListParameterConverter
    extends JodaBeanParameterConverter<TradeList> {

  @Override
  protected Class<TradeList> getExpectedType() {
    return TradeList.class;
  }

}
