/**
 * Copyright (C) 2016 - present by OpenGamma Inc. and the OpenGamma group of companies
 *
 * Please see distribution for license.
 */
package com.opengamma.strata.pricer.capfloor;

import static com.opengamma.strata.basics.date.DayCounts.ACT_ACT_ISDA;
import static com.opengamma.strata.basics.index.IborIndices.USD_LIBOR_3M;
import static org.testng.Assert.assertEquals;

import java.time.Period;
import java.util.List;

import org.testng.annotations.Test;

import com.google.common.collect.ImmutableList;
import com.opengamma.strata.collect.array.DoubleArray;
import com.opengamma.strata.collect.array.DoubleMatrix;
import com.opengamma.strata.collect.tuple.Pair;
import com.opengamma.strata.market.ValueType;
import com.opengamma.strata.market.curve.interpolator.CurveExtrapolators;
import com.opengamma.strata.market.curve.interpolator.CurveInterpolators;
import com.opengamma.strata.market.surface.ConstantSurface;
import com.opengamma.strata.market.surface.Surfaces;
import com.opengamma.strata.pricer.model.SabrVolatilityFormula;
import com.opengamma.strata.pricer.option.RawOptionData;
import com.opengamma.strata.product.capfloor.ResolvedIborCapFloorLeg;

/**
 * Test {@link SabrIborCapletFloorletVolatilityCalibrator}.
 */
@Test
public class SabrIborCapletFloorletVolatilityCalibratorTest
    extends CapletStrippingSetup {

  private static final SabrIborCapletFloorletVolatilityCalibrator CALIBRATOR =
      SabrIborCapletFloorletVolatilityCalibrator.DEFAULT;
  private static final double TOL = 1.0e-3;

  public void test_recovery_black() {
    // choose nodes close to expiries of caps - 0.25y before end dates
    DoubleArray alphaKnots = DoubleArray.of(0.75, 1.75, 2.75, 4.75, 6.75, 9.75);
    DoubleArray rhoKnots = DoubleArray.of(0.75, 2.75, 4.75);
    DoubleArray nuKnots = DoubleArray.of(0.75, 1.75, 2.75, 4.75, 6.75, 9.75);
    SabrIborCapletFloorletVolatilityCalibrationDefinition definition =
        SabrIborCapletFloorletVolatilityCalibrationDefinition.of(
            IborCapletFloorletVolatilitiesName.of("test"), USD_LIBOR_3M, ACT_ACT_ISDA, 0.7,
            alphaKnots, rhoKnots, nuKnots,
            CurveInterpolators.DOUBLE_QUADRATIC, CurveExtrapolators.FLAT, CurveExtrapolators.FLAT,
            SabrVolatilityFormula.hagan());
    ImmutableList<Period> maturities = createBlackMaturities();
    DoubleArray strikes = createBlackStrikes();
    DoubleMatrix volData = createFullBlackDataMatrix();
    DoubleMatrix error = DoubleMatrix.filled(volData.rowCount(), volData.columnCount(), 0.01);
    RawOptionData data = RawOptionData.of(
        maturities, strikes, ValueType.STRIKE, volData, error, ValueType.BLACK_VOLATILITY);
    IborCapletFloorletVolatilityCalibrationResult res = CALIBRATOR.calibrate(definition, CALIBRATION_TIME, data, RATES_PROVIDER);
    SabrParametersIborCapletFloorletVolatilities resVols =
        (SabrParametersIborCapletFloorletVolatilities) res.getVolatilities();
    for (int i = 0; i < NUM_BLACK_STRIKES; ++i) {
      Pair<List<ResolvedIborCapFloorLeg>, List<Double>> capsAndVols = getCapsBlackVols(i);
      List<ResolvedIborCapFloorLeg> caps = capsAndVols.getFirst();
      List<Double> vols = capsAndVols.getSecond();
      int nCaps = caps.size();
      for (int j = 0; j < nCaps; ++j) {
        ConstantSurface volSurface = ConstantSurface.of(
            Surfaces.blackVolatilityByExpiryStrike("test", ACT_ACT_ISDA), vols.get(j));
        BlackIborCapletFloorletExpiryStrikeVolatilities constVol = BlackIborCapletFloorletExpiryStrikeVolatilities.of(
            USD_LIBOR_3M, CALIBRATION_TIME, volSurface);
        double priceOrg = LEG_PRICER_BLACK.presentValue(caps.get(j), RATES_PROVIDER, constVol).getAmount();
        double priceCalib = LEG_PRICER_SABR.presentValue(caps.get(j), RATES_PROVIDER, resVols).getAmount();
        assertEquals(priceOrg, priceCalib, Math.max(priceOrg, 1d) * TOL * 3d);
      }
    }
    assertEquals(resVols.getIndex(), USD_LIBOR_3M);
    assertEquals(resVols.getName(), definition.getName());
    assertEquals(resVols.getValuationDateTime(), CALIBRATION_TIME);
    assertEquals(resVols.getParameterCount(), alphaKnots.size() + rhoKnots.size() + nuKnots.size() + 2); // beta, shift counted
    assertEquals(resVols.getParameters().getShiftCurve(), definition.getShiftCurve());
    assertEquals(resVols.getParameters().getBetaCurve(), definition.getBetaCurve().get());
  }

  public void test_recovery_black_shift() {
    // choose nodes close to expiries of caps - 0.25y before end dates
    DoubleArray alphaKnots = DoubleArray.of(0.75, 1.75, 2.75, 4.75, 6.75, 9.75);
    DoubleArray rhoKnots = DoubleArray.of(0.75, 2.75, 4.75);
    DoubleArray nuKnots = DoubleArray.of(0.75, 1.75, 2.75, 4.75, 6.75, 9.75);
    SabrIborCapletFloorletVolatilityCalibrationDefinition definition =
        SabrIborCapletFloorletVolatilityCalibrationDefinition.of(
            IborCapletFloorletVolatilitiesName.of("test"), USD_LIBOR_3M, ACT_ACT_ISDA, 0.7, 0.05, alphaKnots, rhoKnots, nuKnots,
            CurveInterpolators.PCHIP, CurveExtrapolators.FLAT, CurveExtrapolators.FLAT, SabrVolatilityFormula.hagan(),
            DoubleArray.of(0.03, 0.7, -0.2, 0.9));
    ImmutableList<Period> maturities = createBlackMaturities();
    DoubleArray strikes = createBlackStrikes();
    DoubleMatrix volData =createFullBlackDataMatrix();
    DoubleMatrix error = DoubleMatrix.filled(volData.rowCount(), volData.columnCount(), 0.1);
    RawOptionData data = RawOptionData.of(
        maturities, strikes, ValueType.STRIKE, volData, error, ValueType.BLACK_VOLATILITY);
    IborCapletFloorletVolatilityCalibrationResult res = CALIBRATOR.calibrate(definition, CALIBRATION_TIME, data, RATES_PROVIDER);
    SabrIborCapletFloorletVolatilities resVols =
        (SabrIborCapletFloorletVolatilities) res.getVolatilities();
    for (int i = 0; i < NUM_BLACK_STRIKES; ++i) {
      Pair<List<ResolvedIborCapFloorLeg>, List<Double>> capsAndVols = getCapsBlackVols(i);
      List<ResolvedIborCapFloorLeg> caps = capsAndVols.getFirst();
      List<Double> vols = capsAndVols.getSecond();
      int nCaps = caps.size();
      for (int j = 0; j < nCaps; ++j) {
        ConstantSurface volSurface = ConstantSurface.of(
            Surfaces.blackVolatilityByExpiryStrike("test", ACT_ACT_ISDA), vols.get(j));
        BlackIborCapletFloorletExpiryStrikeVolatilities constVol = BlackIborCapletFloorletExpiryStrikeVolatilities.of(
            USD_LIBOR_3M, CALIBRATION_TIME, volSurface);
        double priceOrg = LEG_PRICER_BLACK.presentValue(caps.get(j), RATES_PROVIDER, constVol).getAmount();
        double priceCalib = LEG_PRICER_SABR.presentValue(caps.get(j), RATES_PROVIDER, resVols).getAmount();
        assertEquals(priceOrg, priceCalib, Math.max(priceOrg, 1d) * TOL * 5d);
      }
    }
  }

  public void recovery_test_flat() {
    // choose nodes close to expiries of caps - 0.25y before end dates
    DoubleArray alphaKnots = DoubleArray.of(0.75, 1.75, 2.75, 4.75, 6.75, 9.75);
    DoubleArray rhoKnots = DoubleArray.of(0.75, 2.75, 4.75);
    DoubleArray nuKnots = DoubleArray.of(0.75, 1.75, 2.75, 4.75, 6.75, 9.75);
    SabrIborCapletFloorletVolatilityCalibrationDefinition definition =
        SabrIborCapletFloorletVolatilityCalibrationDefinition.of(
            IborCapletFloorletVolatilitiesName.of("test"), USD_LIBOR_3M, ACT_ACT_ISDA, 0.95,
            alphaKnots, rhoKnots, nuKnots,
            CurveInterpolators.LINEAR, CurveExtrapolators.FLAT, CurveExtrapolators.FLAT, SabrVolatilityFormula.hagan(),
            DoubleArray.of(0.5, 0.95, 0.0, 0.05));
    DoubleArray strikes = createBlackStrikes();
    RawOptionData data = RawOptionData.of(
        createBlackMaturities(), strikes, ValueType.STRIKE, createFullFlatBlackDataMatrix(), ValueType.BLACK_VOLATILITY);
    IborCapletFloorletVolatilityCalibrationResult res = CALIBRATOR.calibrate(definition, CALIBRATION_TIME, data, RATES_PROVIDER);
    SabrIborCapletFloorletVolatilities resVol =
        (SabrIborCapletFloorletVolatilities) res.getVolatilities();
    for (int i = 0; i < NUM_BLACK_STRIKES; ++i) {
      Pair<List<ResolvedIborCapFloorLeg>, List<Double>> capsAndVols = getCapsFlatBlackVols(i);
      List<ResolvedIborCapFloorLeg> caps = capsAndVols.getFirst();
      List<Double> vols = capsAndVols.getSecond();
      int nCaps = caps.size();
      for (int j = 0; j < nCaps; ++j) {
        ConstantSurface volSurface = ConstantSurface.of(
            Surfaces.blackVolatilityByExpiryStrike("test", ACT_ACT_ISDA), vols.get(j));
        BlackIborCapletFloorletExpiryStrikeVolatilities constVol = BlackIborCapletFloorletExpiryStrikeVolatilities.of(
            USD_LIBOR_3M, CALIBRATION_TIME, volSurface);
        double priceOrg = LEG_PRICER_BLACK.presentValue(caps.get(j), RATES_PROVIDER, constVol).getAmount();
        double priceCalib = LEG_PRICER_SABR.presentValue(caps.get(j), RATES_PROVIDER, resVol).getAmount();
        assertEquals(priceOrg, priceCalib, Math.max(priceOrg, 1d) * TOL);
      }
    }
//    print(res, strikes, 10d);
  }

  public void test_recovery_normal() {
    // choose nodes close to expiries of caps - 0.25y before end dates
    DoubleArray alphaKnots = DoubleArray.of(0.75, 1.75, 2.75, 4.75, 6.75, 9.75);
    DoubleArray rhoKnots = DoubleArray.of(0.75, 2.75, 4.75);
    DoubleArray nuKnots = DoubleArray.of(0.75, 1.75, 2.75, 4.75, 6.75, 9.75);

    SabrIborCapletFloorletVolatilityCalibrationDefinition definition =
        SabrIborCapletFloorletVolatilityCalibrationDefinition.of(
            IborCapletFloorletVolatilitiesName.of("test"), USD_LIBOR_3M, ACT_ACT_ISDA, 0.7,
            alphaKnots, rhoKnots, nuKnots,
            CurveInterpolators.DOUBLE_QUADRATIC, CurveExtrapolators.FLAT, CurveExtrapolators.FLAT,
            SabrVolatilityFormula.hagan());
    ImmutableList<Period> maturities = createNormalEquivMaturities();
    DoubleArray strikes = createNormalEquivStrikes();
    DoubleMatrix volData = createFullNormalEquivDataMatrix();
    DoubleMatrix error = DoubleMatrix.filled(volData.rowCount(), volData.columnCount(), 0.01);
    RawOptionData data = RawOptionData.of(
        maturities, strikes, ValueType.STRIKE, volData, error, ValueType.NORMAL_VOLATILITY);
    IborCapletFloorletVolatilityCalibrationResult res = CALIBRATOR.calibrate(definition, CALIBRATION_TIME, data, RATES_PROVIDER);
    SabrIborCapletFloorletVolatilities resVols =
        (SabrIborCapletFloorletVolatilities) res.getVolatilities();
    for (int i = 0; i < strikes.size(); ++i) {
      Pair<List<ResolvedIborCapFloorLeg>, List<Double>> capsAndVols = getCapsNormalEquivVols(i);
      List<ResolvedIborCapFloorLeg> caps = capsAndVols.getFirst();
      List<Double> vols = capsAndVols.getSecond();
      int nCaps = caps.size();
      for (int j = 0; j < nCaps; ++j) {
        ConstantSurface volSurface = ConstantSurface.of(
            Surfaces.normalVolatilityByExpiryStrike("test", ACT_ACT_ISDA), vols.get(j));
        NormalIborCapletFloorletExpiryStrikeVolatilities constVol = NormalIborCapletFloorletExpiryStrikeVolatilities.of(
            USD_LIBOR_3M, CALIBRATION_TIME, volSurface);
        double priceOrg = LEG_PRICER_NORMAL.presentValue(caps.get(j), RATES_PROVIDER, constVol).getAmount();
        double priceCalib = LEG_PRICER_SABR.presentValue(caps.get(j), RATES_PROVIDER, resVols).getAmount();
//        System.out.println(priceOrg + "\t" + priceCalib);
        assertEquals(priceOrg, priceCalib, Math.max(priceOrg, 1d) * TOL * 3d);
      }
    }

//    print(res, strikes, 10d);
//    assertTrue(res.getChiSquare() > 0d);
//    assertEquals(resVols.getIndex(), USD_LIBOR_3M);
//    assertEquals(resVols.getName(), definition.getName());
//    assertEquals(resVols.getValuationDateTime(), CALIBRATION_TIME);
//    assertEquals(resVols.getParameters().getShiftCurve(), definition.getShiftCurve());
//    assertEquals(resVols.getParameters().getBetaCurve(), definition.getBetaCurve().get());
  }

}
