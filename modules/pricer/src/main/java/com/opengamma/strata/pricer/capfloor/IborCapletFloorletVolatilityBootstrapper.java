package com.opengamma.strata.pricer.capfloor;

import java.time.LocalDate;
import java.time.Period;
import java.time.ZonedDateTime;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

import com.opengamma.strata.basics.ReferenceData;
import com.opengamma.strata.basics.date.DayCount;
import com.opengamma.strata.basics.index.IborIndex;
import com.opengamma.strata.collect.ArgChecker;
import com.opengamma.strata.collect.array.DoubleArray;
import com.opengamma.strata.market.surface.ConstantSurface;
import com.opengamma.strata.market.surface.InterpolatedNodalSurface;
import com.opengamma.strata.market.surface.Surface;
import com.opengamma.strata.market.surface.SurfaceMetadata;
import com.opengamma.strata.market.surface.Surfaces;
import com.opengamma.strata.pricer.impl.option.GenericImpliedVolatiltySolver;
import com.opengamma.strata.pricer.option.RawOptionData;
import com.opengamma.strata.pricer.rate.RatesProvider;
import com.opengamma.strata.product.capfloor.ResolvedIborCapFloorLeg;

public class IborCapletFloorletVolatilityBootstrapper extends IborCapletFloorletVolatilityCalibrator {

  public static final IborCapletFloorletVolatilityBootstrapper DEFAULT =
      new IborCapletFloorletVolatilityBootstrapper(VolatilityIborCapFloorLegPricer.DEFAULT, ReferenceData.standard());

  public static IborCapletFloorletVolatilityBootstrapper of(VolatilityIborCapFloorLegPricer pricer, ReferenceData referenceData) {
    return new IborCapletFloorletVolatilityBootstrapper(pricer, referenceData);
  }

  private IborCapletFloorletVolatilityBootstrapper(VolatilityIborCapFloorLegPricer pricer, ReferenceData referenceData) {
    super(pricer, referenceData);
  }

//  public BlackIborCapletFloorletVolatilities calibrate(
//      IborCapletFloorletDefinition definition,
//      ZonedDateTime calibrationDateTime,
//      RawOptionData capFloorData,
//      RatesProvider ratesProvider) {
//
//    IborIndex index = definition.getIndex();
//    LocalDate baseDate = index.getEffectiveDateOffset().adjust(calibrationDateTime.toLocalDate(), referenceData);
//    LocalDate startDate = baseDate.plus(index.getTenor());
//
//    DoubleArray timesTotal = DoubleArray.EMPTY;
//    DoubleArray strikesTotal = DoubleArray.EMPTY;
//    DoubleArray volsTotal = DoubleArray.EMPTY;
//
//    DoubleArray strikes = capFloorData.getStrikes();
//    int nStrikes = strikes.size();
//    List<Period> expiries = capFloorData.getExpiries();
//    int nExpiries = expiries.size();
//    for (int i = 0; i < nStrikes; ++i) {
//      List<Double> marketValues = new ArrayList<>();
//      List<ResolvedIborCapFloorLeg> capFloors = new ArrayList<>();
//      for (int j = 0; j < nExpiries; ++j) {
//        if (Double.isFinite(capFloorData.getData().get(j, i))) {
//          marketValues.add(capFloorData.getData().get(j, i));
//          LocalDate endDate = baseDate.plus(expiries.get(j));
//          IborCapFloorLeg capFloor = IborCapFloorLeg.builder()
//              .calculation(IborRateCalculation.of(index))
//              .capSchedule(ValueSchedule.of(strikes.get(i)))
//              .currency(index.getCurrency())
//              .notional(ValueSchedule.ALWAYS_1)
//              .paymentSchedule(
//                  PeriodicSchedule.of(startDate, endDate, Frequency.of(index.getTenor().getPeriod()),
//                      definition.getBusinessDayAdjustment(), StubConvention.NONE, RollConventions.NONE))
//              .payReceive(PayReceive.RECEIVE)
//              .build();
//          capFloors.add(capFloor.resolve(referenceData));
//        }
//      }
//      BlackIborCapletFloorletExpiryStrikeVolatilities volsForStrike =
//          calibrate(definition, calibrationDateTime, capFloors, marketValues, ratesProvider);
//      InterpolatedNodalSurface nodalSurface = (InterpolatedNodalSurface) volsForStrike.getSurface();
//      timesTotal = timesTotal.concat(nodalSurface.getXValues());
//      strikesTotal = strikesTotal.concat(nodalSurface.getYValues());
//      volsTotal = volsTotal.concat(nodalSurface.getZValues());
//    }
//
//    SurfaceInterpolator interp = GridSurfaceInterpolator.of(CurveInterpolators.LINEAR, CurveExtrapolators.FLAT,
//        CurveInterpolators.LINEAR, CurveExtrapolators.FLAT); // TODO from definition
//    SurfaceMetadata metadata = Surfaces.blackVolatilityByExpiryStrike(definition.getName().getName(), definition.getDayCount());
//    InterpolatedNodalSurface surface = InterpolatedNodalSurface.ofUnsorted(
//        metadata, timesTotal, strikesTotal, volsTotal, interp);
//    return BlackIborCapletFloorletExpiryStrikeVolatilities.of(definition.getIndex(), calibrationDateTime, surface);
//  }


  @Override
  public IborCapletFloorletVolatilities calibrate(
      IborCapletFloorletDefinition definition,
      ZonedDateTime calibrationDateTime,
      RawOptionData capFloorData,
      RatesProvider ratesProvider) {

    ArgChecker.isTrue(definition instanceof IborCapletFloorletBootstrapDefinition);
    IborCapletFloorletBootstrapDefinition bootstrapDefinition = (IborCapletFloorletBootstrapDefinition) definition;
    IborIndex index = bootstrapDefinition.getIndex();
    LocalDate calibrationDate = calibrationDateTime.toLocalDate();
    LocalDate baseDate = index.getEffectiveDateOffset().adjust(calibrationDate, referenceData);
    LocalDate startDate = baseDate.plus(index.getTenor());
    Function<Surface, IborCapletFloorletVolatilities> volatilitiesFunction = volatilitiesFunction(
        bootstrapDefinition, calibrationDateTime, capFloorData);
    SurfaceMetadata metaData = bootstrapDefinition.createMetadata(capFloorData);
    List<Period> expiries = capFloorData.getExpiries();
    int nExpiries = expiries.size();

    List<Double> timeList = new ArrayList<>();
    List<Double> strikeList = new ArrayList<>();
    List<Double> volList = new ArrayList<>();
    List<ResolvedIborCapFloorLeg> capList = new ArrayList<>();
    List<Double> priceList = new ArrayList<>();
    int[] startIndex = new int[nExpiries + 1];
    for (int i = 0; i < nExpiries; ++i) {
      LocalDate endDate = baseDate.plus(expiries.get(i));
      DoubleArray volatilityData = capFloorData.getData().row(i);
      reduceRawData(
          bootstrapDefinition,
          ratesProvider,
          capFloorData.getStrikes(),
          volatilityData,
          startDate,
          endDate,
          metaData,
          volatilitiesFunction,
          timeList,
          strikeList,
          volList,
          capList,
          priceList);
      startIndex[i + 1] = volList.size();
    }

    InterpolatedNodalSurface surface = InterpolatedNodalSurface.of(
        metaData,
        DoubleArray.copyOf(timeList),
        DoubleArray.copyOf(strikeList),
        DoubleArray.copyOf(volList),
        bootstrapDefinition.getInterpolator());
    IborCapletFloorletVolatilities vols = volatilitiesFunction.apply(surface);
    for (int i = 1; i < nExpiries; ++i) {
      for (int j = startIndex[i]; j < startIndex[i + 1]; ++j) {
        Function<Double, double[]> func = getValueVegaFunction(capList.get(j), ratesProvider, vols, j);
        GenericImpliedVolatiltySolver solver = new GenericImpliedVolatiltySolver(func);
        double capletVol = solver.impliedVolatility(priceList.get(j), volList.get(j));
        vols = vols.withParameter(j, capletVol);
      }
    }
    return vols;
  }


  public BlackIborCapletFloorletExpiryStrikeVolatilities calibrate(
      IborCapletFloorletBootstrapDefinition definition,
      ZonedDateTime calibrationDateTime,
      List<ResolvedIborCapFloorLeg> capFloors,
      List<Double> volatilities,
      RatesProvider ratesProvider) {

    // strike common, start common
    // capFloors size == volatilities size
    // notional == 1
    // val date in ratesProvider == calib date

    int nNodes = capFloors.size();
    LocalDate calibrationDate = calibrationDateTime.toLocalDate();
    DayCount dayCount = definition.getDayCount();
    DoubleArray capVolatilities = DoubleArray.of(nNodes, i -> volatilities.get(i));
    DoubleArray expirationTimes = DoubleArray.of(nNodes, i -> dayCount.relativeYearFraction(calibrationDate,
        capFloors.get(i).getFinalFixingDate().toLocalDate()));
    DoubleArray strikes =
        DoubleArray.of(nNodes, i -> capFloors.get(i).getCapletFloorletPeriods().get(i).getStrike()); // TODO check common strike?
    SurfaceMetadata metadata = Surfaces.blackVolatilityByExpiryStrike(definition.getName().getName(), dayCount);
    InterpolatedNodalSurface surface = InterpolatedNodalSurface.of(
        metadata, expirationTimes, strikes, capVolatilities, definition.getInterpolator());
    BlackIborCapletFloorletExpiryStrikeVolatilities vols =
        BlackIborCapletFloorletExpiryStrikeVolatilities.of(definition.getIndex(), calibrationDateTime, surface);
    DoubleArray prices = DoubleArray.of(
        capFloors.size(),
        i -> pricer.presentValue(
            capFloors.get(i),
            ratesProvider,
            BlackIborCapletFloorletExpiryStrikeVolatilities.of(
                definition.getIndex(), calibrationDateTime, ConstantSurface.of(metadata, volatilities.get(i))))
            .getAmount());
    for (int i = 0; i < nNodes; ++i) {
      Function<Double, double[]> func = getValueVegaFunction(capFloors.get(i), ratesProvider, vols, i);
      GenericImpliedVolatiltySolver solver = new GenericImpliedVolatiltySolver(func);
      double capletVol = solver.impliedVolatility(prices.get(i), volatilities.get(i));
      vols = vols.withParameter(i, capletVol);
//      double prRes = legPricer.presentValue(capFloors.get(i), ratesProvider, vols).getAmount();
//      System.out.println(prices.get(i) + "\t" + prRes);
    }
//    System.out.println("");
    return vols;
  }

  private Function<Double, double[]> getValueVegaFunction(
      ResolvedIborCapFloorLeg leg,
      RatesProvider ratesProvider,
      IborCapletFloorletVolatilities vols,
      int nodeIndex) {

    Function<Double, double[]> priceAndVegaFunction = new Function<Double, double[]>() {
      @Override
      public double[] apply(Double x) {
        IborCapletFloorletVolatilities newVols = vols.withParameter(nodeIndex, x);
        double price = pricer.presentValue(leg, ratesProvider, newVols).getAmount();
        double vega = pricer.presentValueSensitivityModelParamsVolatility(leg, ratesProvider, newVols).build()
            .getSensitivities().get(0).getSensitivity();
        return new double[] {price, vega};
      }
    };
    return priceAndVegaFunction;
  }

}
