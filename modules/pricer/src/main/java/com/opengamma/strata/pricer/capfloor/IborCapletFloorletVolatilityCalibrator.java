package com.opengamma.strata.pricer.capfloor;

import java.time.LocalDate;
import java.time.ZonedDateTime;
import java.util.List;
import java.util.function.Function;

import com.opengamma.strata.basics.ReferenceData;
import com.opengamma.strata.basics.index.IborIndex;
import com.opengamma.strata.basics.schedule.Frequency;
import com.opengamma.strata.basics.schedule.PeriodicSchedule;
import com.opengamma.strata.basics.schedule.RollConventions;
import com.opengamma.strata.basics.schedule.StubConvention;
import com.opengamma.strata.basics.value.ValueSchedule;
import com.opengamma.strata.collect.ArgChecker;
import com.opengamma.strata.collect.array.DoubleArray;
import com.opengamma.strata.market.ValueType;
import com.opengamma.strata.market.surface.ConstantSurface;
import com.opengamma.strata.market.surface.Surface;
import com.opengamma.strata.market.surface.SurfaceMetadata;
import com.opengamma.strata.pricer.option.RawOptionData;
import com.opengamma.strata.pricer.rate.RatesProvider;
import com.opengamma.strata.product.capfloor.IborCapFloorLeg;
import com.opengamma.strata.product.capfloor.ResolvedIborCapFloorLeg;
import com.opengamma.strata.product.common.PayReceive;
import com.opengamma.strata.product.swap.IborRateCalculation;

abstract class IborCapletFloorletVolatilityCalibrator {

  // TODO jacobian

  // TODO shift in raw data

  protected final ReferenceData referenceData;
  protected final VolatilityIborCapFloorLegPricer pricer;

  public IborCapletFloorletVolatilityCalibrator(VolatilityIborCapFloorLegPricer pricer, ReferenceData referenceData) {
    this.pricer = ArgChecker.notNull(pricer, "pricer");
    this.referenceData = ArgChecker.notNull(referenceData, "referenceData");
  }

  public abstract IborCapletFloorletVolatilities calibrate(
      IborCapletFloorletDefinition definition,
      ZonedDateTime calibrationDateTime,
      RawOptionData capFloorData,
      RatesProvider ratesProvider);

  protected void reduceRawData(
      IborCapletFloorletDefinition definition,
      RatesProvider ratesProvider,
      DoubleArray strikes,
      DoubleArray volatilityData,
      LocalDate startDate,
      LocalDate endDate,
      SurfaceMetadata metadata,
      Function<Surface, IborCapletFloorletVolatilities> volatilityFunction,
      List<Double> timeList,
      List<Double> strikeList,
      List<Double> volList,
      List<ResolvedIborCapFloorLeg> capList,
      List<Double> priceList) {

    IborIndex index = definition.getIndex();
    int nStrikes = strikes.size();
    for (int j = 0; j < nStrikes; ++j) {
      if (Double.isFinite(volatilityData.get(j))) {
        ResolvedIborCapFloorLeg capFloor = IborCapFloorLeg.builder()
            .calculation(IborRateCalculation.of(index))
            .capSchedule(ValueSchedule.of(strikes.get(j)))
            .currency(index.getCurrency())
            .notional(ValueSchedule.ALWAYS_1)
            .paymentSchedule(
                PeriodicSchedule.of(
                    startDate,
                    endDate,
                    Frequency.of(index.getTenor().getPeriod()),
                    definition.getBusinessDayAdjustment(),
                    StubConvention.NONE,
                    RollConventions.NONE))
            .payReceive(PayReceive.RECEIVE)
            .build()
            .resolve(referenceData);
        capList.add(capFloor);
        strikeList.add(strikes.get(j));
        volList.add(volatilityData.get(j));
        ConstantSurface constVolSurface = ConstantSurface.of(metadata, volatilityData.get(j));
        IborCapletFloorletVolatilities vols = volatilityFunction.apply(constVolSurface);
        timeList.add(vols.relativeTime(capFloor.getFinalFixingDate()));
        priceList.add(pricer.presentValue(capFloor, ratesProvider, vols).getAmount());
      }
    }

  }

  protected Function<Surface, IborCapletFloorletVolatilities> volatilitiesFunction(
      IborCapletFloorletDefinition definition,
      ZonedDateTime calibrationDateTime,
      RawOptionData capFloorData) {

    IborIndex index = definition.getIndex();
    if (capFloorData.getStrikeType().equals(ValueType.STRIKE)) {
      if (capFloorData.getDataType().equals(ValueType.BLACK_VOLATILITY)) {
        return blackVolatilitiesFunction(index, calibrationDateTime);
      } else if (capFloorData.getDataType().equals(ValueType.NORMAL_VOLATILITY)) {
        return normalVolatilitiesFunction(index, calibrationDateTime);
      }
      throw new IllegalArgumentException("Data type not supported");
    }
    throw new IllegalArgumentException("strike type must be ValueType.STRIKE");
  }

  private Function<Surface, IborCapletFloorletVolatilities> blackVolatilitiesFunction(
      IborIndex index,
      ZonedDateTime calibrationDateTime) {

    Function<Surface, IborCapletFloorletVolatilities> func = new Function<Surface, IborCapletFloorletVolatilities>() {
      @Override
      public IborCapletFloorletVolatilities apply(Surface s) {
        return BlackIborCapletFloorletExpiryStrikeVolatilities.of(index, calibrationDateTime, s);
      }
    };
    return func;
  }

  private Function<Surface, IborCapletFloorletVolatilities> normalVolatilitiesFunction(
      IborIndex index,
      ZonedDateTime calibrationDateTime) {

    Function<Surface, IborCapletFloorletVolatilities> func = new Function<Surface, IborCapletFloorletVolatilities>() {
      @Override
      public IborCapletFloorletVolatilities apply(Surface s) {
        return NormalIborCapletFloorletExpiryStrikeVolatilities.of(index, calibrationDateTime, s);
      }
    };
    return func;
  }

}
