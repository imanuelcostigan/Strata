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

@Test
public class SabrIborCapletFloorletVolatilityCalibratorTest
    extends CapletStrippingSetup {

  private static final SabrIborCapletFloorletVolatilityCalibrator CALIBRATOR =
      SabrIborCapletFloorletVolatilityCalibrator.DEFAULT;
  private static final BlackIborCapFloorLegPricer LEG_PRICER_BLACK = BlackIborCapFloorLegPricer.DEFAULT;
  private static final NormalIborCapFloorLegPricer LEG_PRICER_NORMAL = NormalIborCapFloorLegPricer.DEFAULT;
  private static final SabrIborCapFloorLegPricer LEG_PRICER_SABR = SabrIborCapFloorLegPricer.DEFAULT;
  private static final double TOL = 1.0e-3;

  public void test_recovery_black() {

    DoubleArray alphaKnots = DoubleArray.of(1, 2, 3, 5, 7, 10);
    DoubleArray betaKnots = DoubleArray.of(1);
    DoubleArray rhoKnots = DoubleArray.of(1, 3, 7);
    DoubleArray nuKnots = DoubleArray.of(1, 2, 3, 5, 7, 10);
    double error = 0.001;

    SabrIborCapletFloorletCalibrationDefinition definition = SabrIborCapletFloorletCalibrationDefinition.of(
        IborCapletFloorletVolatilitiesName.of("test"), USD_LIBOR_3M, ACT_ACT_ISDA,
        ImmutableList.of(alphaKnots, betaKnots, rhoKnots, nuKnots),
        CurveInterpolators.DOUBLE_QUADRATIC, CurveExtrapolators.FLAT, CurveExtrapolators.FLAT,
        SabrVolatilityFormula.hagan());
    ImmutableList<Period> maturities = createBlackMaturities();
    DoubleArray strikes = createBlackStrikes();
    RawOptionData data = RawOptionData.of(
        maturities, strikes, ValueType.STRIKE, createFullBlackDataMatrix(), ValueType.BLACK_VOLATILITY,
        DoubleMatrix.filled(maturities.size(), strikes.size(), error));
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
        System.out.println(priceOrg + "\t" + priceCalib);
        assertEquals(priceOrg, priceCalib, Math.max(priceOrg, 1d) * TOL * 5d);
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
