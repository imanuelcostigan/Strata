/**
 * Copyright (C) 2015 - present by OpenGamma Inc. and the OpenGamma group of companies
 *
 * Please see distribution for license.
 */
package com.opengamma.strata.pricer.curve;

import static java.time.temporal.ChronoUnit.MONTHS;

import java.time.YearMonth;

import com.opengamma.analytics.math.curve.InterpolatedDoublesCurve;

/**
 * Implementation of a PriceIndexCurve where the curve is represented by an interpolated curve.
 */
public class PriceIndexInterpolatedSeasonalityCurve implements PriceIndexCurve {
  
  /** The doubles curve underlying the price index curve. 
   * The X dimension on the curve represent the number of months between. **/
  private final InterpolatedDoublesCurve curve;
  /** The reference month to compute the number of months from. **/
  private final YearMonth valuationMonth;
  /** Describes the seasonal adjustments. The array has a dimension of 12, one element for each month, 
   * starting from January. The adjustments are multiplicative. For each month, the price index is the one obtained
   * from the interpolated part of the curve multiplied by the seasonal adjustment. */
  private final double[] seasonality;

  /**
   * Constructor.
   * @param curve  the underlying curve
   * @param seasonality  the seasonal adjustment for each month, starting with January.
   * @param referenceMonth  the reference month
   */
  public PriceIndexInterpolatedSeasonalityCurve(InterpolatedDoublesCurve curve, double[] seasonality, 
      YearMonth referenceMonth) {
    this.curve = curve;
    this.seasonality = seasonality;
    this.valuationMonth = referenceMonth;
  }

  @Override
  public String getName() {
    return curve.getName();
  }

  @Override
  public double getPriceIndex(YearMonth month) {
    int monthInt = month.getMonthValue();
    double nbMonth = valuationMonth.until(month, MONTHS);
    double indexInterpolated = curve.getYValue(nbMonth);
    double adjustment = seasonality[monthInt];
    return indexInterpolated * adjustment;
  }

  @Override
  public Double[] getPriceIndexParameterSensitivity(YearMonth month) {
    double nbMonth = valuationMonth.until(month, MONTHS); // TODO: review - multiply by seasonality
    return curve.getYValueParameterSensitivity(nbMonth);
  }

  @Override
  public int getNumberOfParameters() {
    return curve.size();
  }

  @Override
  public PriceIndexCurve shiftCurve(double[] shifts) {
    double[] x = curve.getXDataAsPrimitive();
    double[] y = curve.getYDataAsPrimitive();
    double[] yShifted = new double[y.length];
    for(int i=0; i<y.length; i++) {
      yShifted[i] = y[i] + shifts[i];
    }
    InterpolatedDoublesCurve curveShifted = new InterpolatedDoublesCurve(x, yShifted, curve.getInterpolator(), true);
    return new PriceIndexInterpolatedSeasonalityCurve(curveShifted, seasonality, valuationMonth);
  }

}
