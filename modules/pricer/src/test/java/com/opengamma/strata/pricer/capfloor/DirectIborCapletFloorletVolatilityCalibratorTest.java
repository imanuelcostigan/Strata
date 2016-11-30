package com.opengamma.strata.pricer.capfloor;

import static com.opengamma.strata.basics.date.DayCounts.ACT_ACT_ISDA;
import static com.opengamma.strata.basics.index.IborIndices.USD_LIBOR_3M;
import static org.testng.Assert.assertEquals;

import java.util.List;

import org.testng.annotations.Test;

import com.opengamma.strata.collect.tuple.Pair;
import com.opengamma.strata.market.ValueType;
import com.opengamma.strata.market.curve.interpolator.CurveInterpolators;
import com.opengamma.strata.market.surface.ConstantSurface;
import com.opengamma.strata.market.surface.Surfaces;
import com.opengamma.strata.market.surface.interpolator.GridSurfaceInterpolator;
import com.opengamma.strata.pricer.option.RawOptionData;
import com.opengamma.strata.product.capfloor.ResolvedIborCapFloorLeg;

@Test
public class DirectIborCapletFloorletVolatilityCalibratorTest
    extends CapletStrippingSetup {

  private static final DirectIborCapletFloorletVolatilityCalibrator CALIBRATOR =
      DirectIborCapletFloorletVolatilityCalibrator.DEFAULT;
  private static final BlackIborCapFloorLegPricer LEG_PRICER_BLACK = BlackIborCapFloorLegPricer.DEFAULT;
  private static final NormalIborCapFloorLegPricer LEG_PRICER_NORMAL = NormalIborCapFloorLegPricer.DEFAULT;
  private static final double TOL = 1.0e-4;

  public void test_recovery_black() {

    double lambdaT = 0.07;
    double lambdaK = 0.07;
    double error = 1.0e-5;

    DirectIborCapletFloorletDefinition definition = DirectIborCapletFloorletDefinition.of(
        IborCapletFloorletVolatilitiesName.of("test"), USD_LIBOR_3M, ACT_ACT_ISDA, lambdaT, lambdaK, error,
        GridSurfaceInterpolator.of(CurveInterpolators.LINEAR, CurveInterpolators.LINEAR));
    RawOptionData data = RawOptionData.of(
        createBlackMaturities(), createBlackStrikes(), ValueType.STRIKE, createFullBlackDataMatrix(), ValueType.BLACK_VOLATILITY);
    IborCapletFloorletVolatilityCalibrationResult res = CALIBRATOR.calibrate(definition, CALIBRATION_TIME, data, RATES_PROVIDER);
    BlackIborCapletFloorletExpiryStrikeVolatilities resVols =
        (BlackIborCapletFloorletExpiryStrikeVolatilities) res.getVolatilities();
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
        double priceCalib = LEG_PRICER_BLACK.presentValue(caps.get(j), RATES_PROVIDER, resVols).getAmount();
        assertEquals(priceOrg, priceCalib, Math.max(priceOrg, 1d) * TOL * 5d);
      }
    }

//    print(res, createBlackStrikes(), 10d);
//    assertTrue(res.getChiSquare() > 0d);
//    assertEquals(resVols.getIndex(), USD_LIBOR_3M);
//    assertEquals(resVols.getName(), definition.getName());
//    assertEquals(resVols.getValuationDateTime(), CALIBRATION_TIME);
//    assertEquals(resVols.getParameters().getShiftCurve(), definition.getShiftCurve());
//    assertEquals(resVols.getParameters().getBetaCurve(), definition.getBetaCurve().get());
  }

}
