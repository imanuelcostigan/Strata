package com.opengamma.strata.pricer.capfloor;

import static com.opengamma.strata.basics.date.DayCounts.ACT_ACT_ISDA;
import static com.opengamma.strata.basics.index.IborIndices.USD_LIBOR_3M;
import static org.testng.Assert.assertEquals;

import java.util.List;

import org.testng.annotations.Test;

import com.opengamma.strata.collect.array.DoubleArray;
import com.opengamma.strata.collect.tuple.Pair;
import com.opengamma.strata.market.ValueType;
import com.opengamma.strata.market.curve.interpolator.CurveInterpolators;
import com.opengamma.strata.market.surface.ConstantSurface;
import com.opengamma.strata.market.surface.Surfaces;
import com.opengamma.strata.pricer.option.RawOptionData;
import com.opengamma.strata.product.capfloor.ResolvedIborCapFloorLeg;

@Test
public class IborCapletFloorletVolatilityBootstrapperTest extends CapletStrippingSetup {

  // TODO test normal vol

  private static final IborCapletFloorletVolatilityBootstrapper CALIBRATOR = IborCapletFloorletVolatilityBootstrapper.DEFAULT;
  private static final BlackIborCapFloorLegPricer LEG_PRICER_BLACK = BlackIborCapFloorLegPricer.DEFAULT;
  private static final NormalIborCapFloorLegPricer LEG_PRICER_NORMAL = NormalIborCapFloorLegPricer.DEFAULT;
  private static final double TOL = 1.0e-10;

  public void recovery_test() {
    IborCapletFloorletBootstrapDefinition definition =
        IborCapletFloorletBootstrapDefinition.of(IborCapletFloorletVolatilitiesName.of("test"), USD_LIBOR_3M, ACT_ACT_ISDA,
            CurveInterpolators.LINEAR, CurveInterpolators.DOUBLE_QUADRATIC);

    DoubleArray timesTotal = DoubleArray.EMPTY;
    DoubleArray strikesTotal = DoubleArray.EMPTY;
    DoubleArray volsTotal = DoubleArray.EMPTY;
    for (int i = 0; i < NUM_BLACK_STRIKES; ++i) {
      Pair<List<ResolvedIborCapFloorLeg>, List<Double>> capsAndVols = getCapsBlackVols(i);
      List<ResolvedIborCapFloorLeg> caps = capsAndVols.getFirst();
      List<Double> vols = capsAndVols.getSecond();
      BlackIborCapletFloorletExpiryStrikeVolatilities res =
          CALIBRATOR.calibrate(definition, CALIBRATION_TIME, caps, vols, RATES_PROVIDER);
      int nCaps = caps.size();

      for (int j = 0; j < nCaps; ++j) {
        ConstantSurface volSurface = ConstantSurface.of(
            Surfaces.blackVolatilityByExpiryStrike("test", ACT_ACT_ISDA), vols.get(j));
        BlackIborCapletFloorletExpiryStrikeVolatilities constVol = BlackIborCapletFloorletExpiryStrikeVolatilities.of(
            USD_LIBOR_3M, CALIBRATION_TIME, volSurface);
        double priceOrg = LEG_PRICER_BLACK.presentValue(caps.get(j), RATES_PROVIDER, constVol).getAmount();
        double priceCalib = LEG_PRICER_BLACK.presentValue(caps.get(j), RATES_PROVIDER, res).getAmount();
        assertEquals(priceOrg, priceCalib, TOL);
      }
    }
  }

  public void recovery_test2() {
    IborCapletFloorletBootstrapDefinition definition =
        IborCapletFloorletBootstrapDefinition.of(IborCapletFloorletVolatilitiesName.of("test"), USD_LIBOR_3M, ACT_ACT_ISDA,
            CurveInterpolators.LINEAR, CurveInterpolators.LINEAR);
    RawOptionData data =
        RawOptionData.of(createBlackMaturities(), createBlackStrikes(), ValueType.STRIKE, createFullBlackDataMatrix(),
        ValueType.BLACK_VOLATILITY);
    BlackIborCapletFloorletExpiryStrikeVolatilities res =
        (BlackIborCapletFloorletExpiryStrikeVolatilities) CALIBRATOR.calibrate(definition, CALIBRATION_TIME, data,
            RATES_PROVIDER);
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
        double priceCalib = LEG_PRICER_BLACK.presentValue(caps.get(j), RATES_PROVIDER, res).getAmount();
        assertEquals(priceOrg, priceCalib, Math.max(priceOrg, 1d) * TOL);
      }
    }
//
//    final int nSamples = 51;
//    final int nStrikeSamples = 51;
//    System.out.print("\n");
//    DoubleArray strikes = crateBlackStrikes();
//    for (int i = 0; i < nStrikeSamples; i++) {
//      System.out.print("\t" + (strikes.get(0) + (strikes.get(strikes.size() - 1) - strikes.get(0)) * i) / (nStrikeSamples - 1));
//    }
//    System.out.print("\n");
//    for (int index = 0; index < nSamples; index++) {
//      final double t = 0.25 + index * 10.0 / (nSamples - 1);
//      System.out.print(t);
//      for (int i = 0; i < nStrikeSamples; i++) {
//        double strike = (strikes.get(0) + (strikes.get(strikes.size() - 1) - strikes.get(0)) * i) / (nStrikeSamples - 1);
//        System.out.print("\t" + res.getSurface().zValue(t, strike));
//      }
//      System.out.print("\n");
//    }
  }

  public void recovery_test2Flat() {
    IborCapletFloorletBootstrapDefinition definition =
        IborCapletFloorletBootstrapDefinition.of(IborCapletFloorletVolatilitiesName.of("test"), USD_LIBOR_3M, ACT_ACT_ISDA,
            CurveInterpolators.LINEAR, CurveInterpolators.LINEAR);
    RawOptionData data =
        RawOptionData.of(createBlackMaturities(), createBlackStrikes(), ValueType.STRIKE, createFullFlatBlackDataMatrix(),
            ValueType.BLACK_VOLATILITY);
    BlackIborCapletFloorletExpiryStrikeVolatilities res =
        (BlackIborCapletFloorletExpiryStrikeVolatilities) CALIBRATOR.calibrate(definition, CALIBRATION_TIME, data,
            RATES_PROVIDER);
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
        double priceCalib = LEG_PRICER_BLACK.presentValue(caps.get(j), RATES_PROVIDER, res).getAmount();
        assertEquals(priceOrg, priceCalib, Math.max(priceOrg, 1d) * TOL);
//        System.out.println(priceOrg);
      }
    }

//    final int nSamples = 51;
//    final int nStrikeSamples = 51;
//    System.out.print("\n");
//    DoubleArray strikes = crateBlackStrikes();
//    for (int i = 0; i < nStrikeSamples; i++) {
//      System.out.print("\t" + (strikes.get(0) + (strikes.get(strikes.size() - 1) - strikes.get(0)) * i) / (nStrikeSamples - 1));
//    }
//    System.out.print("\n");
//    for (int index = 0; index < nSamples; index++) {
//      final double t = 0.25 + index * 10.0 / (nSamples - 1);
//      System.out.print(t);
//      for (int i = 0; i < nStrikeSamples; i++) {
//        double strike = (strikes.get(0) + (strikes.get(strikes.size() - 1) - strikes.get(0)) * i) / (nStrikeSamples - 1);
//        System.out.print("\t" + res.getSurface().zValue(t, strike));
//      }
//      System.out.print("\n");
//    }
  }

  public void recovery_test3() {
    IborCapletFloorletBootstrapDefinition definition =
        IborCapletFloorletBootstrapDefinition.of(IborCapletFloorletVolatilitiesName.of("test"), USD_LIBOR_3M, ACT_ACT_ISDA,
            CurveInterpolators.LINEAR, CurveInterpolators.DOUBLE_QUADRATIC);
    RawOptionData data =
        RawOptionData.of(createNormalMaturities(), createNormalStrikes(), ValueType.STRIKE, createFullNormalDataMatrix(),
        ValueType.NORMAL_VOLATILITY);
    NormalIborCapletFloorletExpiryStrikeVolatilities res =
        (NormalIborCapletFloorletExpiryStrikeVolatilities) CALIBRATOR.calibrate(definition, CALIBRATION_TIME, data,
            RATES_PROVIDER);
    for (int i = 0; i < NUM_NORMAL_STRIKES; ++i) {
      Pair<List<ResolvedIborCapFloorLeg>, List<Double>> capsAndVols = getCapsNormalVols(i);
      List<ResolvedIborCapFloorLeg> caps = capsAndVols.getFirst();
      List<Double> vols = capsAndVols.getSecond();
      int nCaps = caps.size();
      for (int j = 0; j < nCaps; ++j) {
        ConstantSurface volSurface = ConstantSurface.of(
            Surfaces.normalVolatilityByExpiryStrike("test", ACT_ACT_ISDA), vols.get(j));
        NormalIborCapletFloorletExpiryStrikeVolatilities constVol = NormalIborCapletFloorletExpiryStrikeVolatilities.of(
            USD_LIBOR_3M, CALIBRATION_TIME, volSurface);
        double priceOrg = LEG_PRICER_NORMAL.presentValue(caps.get(j), RATES_PROVIDER, constVol).getAmount();
        double priceCalib = LEG_PRICER_NORMAL.presentValue(caps.get(j), RATES_PROVIDER, res).getAmount();
        assertEquals(priceOrg, priceCalib, Math.max(priceOrg, 1d) * TOL * 10d);
//        System.out.println(priceOrg);
      }
    }

//    final int nSamples = 51;
//    final int nStrikeSamples = 51;
//    System.out.print("\n");
//    DoubleArray strikes = createNormalStrikes();
//    for (int i = 0; i < nStrikeSamples; i++) {
//      System.out.print("\t" + (strikes.get(0) + (strikes.get(strikes.size() - 1) - strikes.get(0)) * i) / (nStrikeSamples - 1));
//    }
//    System.out.print("\n");
//    for (int index = 0; index < nSamples; index++) {
//      final double t = 0.25 + index * 20.0 / (nSamples - 1);
//      double forward = FWD_CURVE.yValue(t);
//      System.out.print(t);
//      for (int i = 0; i < nStrikeSamples; i++) {
//        double strike = (strikes.get(0) + (strikes.get(strikes.size() - 1) - strikes.get(0)) * i) / (nStrikeSamples - 1);
//        System.out.print("\t" + res.volatility(t, strike, forward));
//      }
//      System.out.print("\n");
//    }

  }

  public void recovery_test4() {
    IborCapletFloorletBootstrapDefinition definition =
        IborCapletFloorletBootstrapDefinition.of(IborCapletFloorletVolatilitiesName.of("test"), USD_LIBOR_3M, ACT_ACT_ISDA,
            CurveInterpolators.LINEAR, CurveInterpolators.DOUBLE_QUADRATIC);
    RawOptionData data =
        RawOptionData.of(createNormalEquivMaturities(), createNormalEquivStrikes(), ValueType.STRIKE,
            createFullNormalEquivDataMatrix(),
            ValueType.NORMAL_VOLATILITY);
    NormalIborCapletFloorletExpiryStrikeVolatilities res =
        (NormalIborCapletFloorletExpiryStrikeVolatilities) CALIBRATOR.calibrate(definition, CALIBRATION_TIME, data,
            RATES_PROVIDER);
    for (int i = 0; i < NUM_BLACK_STRIKES; ++i) {
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
        double priceCalib = LEG_PRICER_NORMAL.presentValue(caps.get(j), RATES_PROVIDER, res).getAmount();
        assertEquals(priceOrg, priceCalib, Math.max(priceOrg, 1d) * TOL * 100d);
//        System.out.println(priceOrg);
      }
    }
//
//    final int nSamples = 51;
//    final int nStrikeSamples = 51;
//    System.out.print("\n");
//    DoubleArray strikes = createNormalEquivStrikes();
//    for (int i = 0; i < nStrikeSamples; i++) {
//      System.out.print("\t" + (strikes.get(0) + (strikes.get(strikes.size() - 1) - strikes.get(0)) * i) / (nStrikeSamples - 1));
//    }
//    System.out.print("\n");
//    for (int index = 0; index < nSamples; index++) {
//      final double t = 0.25 + index * 10.0 / (nSamples - 1);
//      double forward = FWD_CURVE.yValue(t);
//      System.out.print(t);
//      for (int i = 0; i < nStrikeSamples; i++) {
//        double strike = (strikes.get(0) + (strikes.get(strikes.size() - 1) - strikes.get(0)) * i) / (nStrikeSamples - 1);
//        System.out.print("\t" + res.volatility(t, strike, forward));
//      }
//      System.out.print("\n");
//    }

  }
}
