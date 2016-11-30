/**
 * Copyright (C) 2016 - present by OpenGamma Inc. and the OpenGamma group of companies
 *
 * Please see distribution for license.
 */
package com.opengamma.strata.pricer.capfloor;

import java.time.LocalDate;
import java.time.Period;
import java.time.ZonedDateTime;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

import com.opengamma.strata.basics.ReferenceData;
import com.opengamma.strata.basics.index.IborIndex;
import com.opengamma.strata.collect.ArgChecker;
import com.opengamma.strata.collect.array.DoubleArray;
import com.opengamma.strata.collect.array.DoubleMatrix;
import com.opengamma.strata.market.surface.InterpolatedNodalSurface;
import com.opengamma.strata.market.surface.Surface;
import com.opengamma.strata.market.surface.SurfaceMetadata;
import com.opengamma.strata.math.impl.linearalgebra.CholeskyDecompositionCommons;
import com.opengamma.strata.math.impl.minimization.PositiveOrZero;
import com.opengamma.strata.math.impl.statistics.leastsquare.LeastSquareResults;
import com.opengamma.strata.math.impl.statistics.leastsquare.NonLinearLeastSquareWithPenalty;
import com.opengamma.strata.pricer.option.RawOptionData;
import com.opengamma.strata.pricer.rate.RatesProvider;
import com.opengamma.strata.product.capfloor.ResolvedIborCapFloorLeg;

public class DirectIborCapletFloorletVolatilityCalibrator
    extends IborCapletFloorletVolatilityCalibrator {

  /**
   * Default implementation. 
   */
  public static final DirectIborCapletFloorletVolatilityCalibrator DEFAULT =
      of(VolatilityIborCapFloorLegPricer.DEFAULT, 1.0e-8, ReferenceData.standard());

  /**
   * The positive function.
   * <p>
   * The function returns true if the new trial position is positive or zero.
   */
  private static final Function<DoubleArray, Boolean> POSITIVE = new PositiveOrZero();
  /**
   * The non-linear square with penalty. 
   */
  private final NonLinearLeastSquareWithPenalty solver;

  //-------------------------------------------------------------------------
  /**
   * Creates an instance. 
   * <p>
   * The epsilon is the parameter used in {@link NonLinearLeastSquareWithPenalty}, where the iteration stops when certain 
   * quantities are smaller than this parameter.
   * 
   * @param pricer  the cap pricer
   * @param epsilon  the epsilon parameter
   * @param referenceData  the reference data
   * @return the instance
   */
  public static DirectIborCapletFloorletVolatilityCalibrator of(
      VolatilityIborCapFloorLegPricer pricer,
      double epsilon,
      ReferenceData referenceData) {

    return new DirectIborCapletFloorletVolatilityCalibrator(pricer, epsilon, referenceData);
  }

  // private constructor
  private DirectIborCapletFloorletVolatilityCalibrator(
      VolatilityIborCapFloorLegPricer pricer,
      double epsilon,
      ReferenceData referenceData) {

    super(pricer, referenceData);
    this.solver = new NonLinearLeastSquareWithPenalty(new CholeskyDecompositionCommons(), epsilon);
  }

  //-------------------------------------------------------------------------
  @Override
  public IborCapletFloorletVolatilityCalibrationResult calibrate(
      IborCapletFloorletDefinition definition,
      ZonedDateTime calibrationDateTime,
      RawOptionData capFloorData,
      RatesProvider ratesProvider) {

    ArgChecker.isTrue(ratesProvider.getValuationDate().equals(calibrationDateTime.toLocalDate()),
        "valuationDate of ratesProvider should be coherent to calibrationDateTime");
    ArgChecker.isTrue(definition instanceof DirectIborCapletFloorletDefinition);
    DirectIborCapletFloorletDefinition directDefinition = (DirectIborCapletFloorletDefinition) definition;
    IborIndex index = directDefinition.getIndex();
    LocalDate calibrationDate = calibrationDateTime.toLocalDate();
    LocalDate baseDate = index.getEffectiveDateOffset().adjust(calibrationDate, referenceData);
    LocalDate startDate = baseDate.plus(index.getTenor());
    Function<Surface, IborCapletFloorletVolatilities> volatilitiesFunction = volatilitiesFunction(
        directDefinition, calibrationDateTime, capFloorData);
    SurfaceMetadata metaData = directDefinition.createMetadata(capFloorData);
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
      reduceRawData(directDefinition, ratesProvider, capFloorData.getStrikes(), volatilityData, startDate, endDate, metaData,
          volatilitiesFunction, timeList, strikeList, volList, capList, priceList);
      startIndex[i + 1] = priceList.size();
    }

    List<Double> timeCapletList = new ArrayList<>();
    List<Double> strikeCapletList = new ArrayList<>();
    List<Double> volCapletList = new ArrayList<>();
    ResolvedIborCapFloorLeg cap = capList.get(capList.size() - 1);
    int nCaplets = cap.getCapletFloorletPeriods().size();
    DoubleArray capletExpiries = DoubleArray.of(nCaplets, n -> directDefinition.getDayCount().relativeYearFraction(
        calibrationDate, cap.getCapletFloorletPeriods().get(n).getFixingDateTime().toLocalDate()));
    DoubleArray strikes = capFloorData.getStrikes();
    createCapletNode(capletExpiries, strikes, timeCapletList, strikeCapletList, volCapletList);
    SurfaceMetadata metadata = directDefinition.createMetadata(capFloorData);
    InterpolatedNodalSurface baseSurface = InterpolatedNodalSurface.of(metadata, DoubleArray.copyOf(timeCapletList),
        DoubleArray.copyOf(strikeCapletList), DoubleArray.copyOf(volCapletList), directDefinition.getInterpolator());

    DoubleMatrix penaltyMatrix = directDefinition.computePenaltyMatrix(strikes, capletExpiries);
    DoubleArray errors = DoubleArray.filled(priceList.size(), directDefinition.getError()); // TODO

    LeastSquareResults res = solver.solve(
        DoubleArray.copyOf(priceList),
        errors,
        getPriceFunction(capList, ratesProvider, volatilitiesFunction, baseSurface),
        getJacobianFunction(capList, ratesProvider, volatilitiesFunction, baseSurface),
        DoubleArray.copyOf(volCapletList),
        penaltyMatrix,
        POSITIVE);
    InterpolatedNodalSurface resSurface = baseSurface.withZValues(res.getFitParameters());
    return IborCapletFloorletVolatilityCalibrationResult.ofLestSquare(volatilitiesFunction.apply(resSurface), res.getChiSq());
  }

  //-------------------------------------------------------------------------
  private void createCapletNode(
      DoubleArray capletExpiries,
      DoubleArray strikes,
      List<Double> timeCapletList,
      List<Double> strikeCapletList,
      List<Double> volCapletList) {

    int nTimes = capletExpiries.size();
    int nStrikes = strikes.size();
    for (int i = 0; i < nTimes; ++i) {
      timeCapletList.addAll(DoubleArray.filled(nStrikes, capletExpiries.get(i)).toList()); // TODO simplify
      strikeCapletList.addAll(strikes.toList());
      volCapletList.addAll(DoubleArray.filled(nStrikes, 0.3d).toList()); // TODO simplify
    }
  }

  private Function<DoubleArray, DoubleArray> getPriceFunction(
      List<ResolvedIborCapFloorLeg> capList,
      RatesProvider ratesProvider,
      Function<Surface, IborCapletFloorletVolatilities> volatilitiesFunction,
      InterpolatedNodalSurface baseSurface) {

    int nCaps = capList.size();
    Function<DoubleArray, DoubleArray> priceFunction = new Function<DoubleArray, DoubleArray>() {
      @Override
      public DoubleArray apply(DoubleArray capletVols) {
        IborCapletFloorletVolatilities newVols = volatilitiesFunction.apply(baseSurface.withZValues(capletVols));
        return DoubleArray.of(nCaps, n -> pricer.presentValue(capList.get(n), ratesProvider, newVols).getAmount());
      }
    };
    return priceFunction;
  }

  private Function<DoubleArray, DoubleMatrix> getJacobianFunction(
      List<ResolvedIborCapFloorLeg> capList,
      RatesProvider ratesProvider,
      Function<Surface, IborCapletFloorletVolatilities> volatilitiesFunction,
      InterpolatedNodalSurface baseSurface) {

    int nCaps = capList.size();
    int nNodes = baseSurface.getParameterCount();
    Function<DoubleArray, DoubleMatrix> jacobianFunction = new Function<DoubleArray, DoubleMatrix>() {
      @Override
      public DoubleMatrix apply(DoubleArray capletVols) {
        IborCapletFloorletVolatilities newVols = volatilitiesFunction.apply(baseSurface.withZValues(capletVols));
        return DoubleMatrix.ofArrayObjects(nCaps, nNodes, n -> newVols.parameterSensitivity(
            pricer.presentValueSensitivityModelParamsVolatility(capList.get(n), ratesProvider, newVols).build())
            .getSensitivities()
            .get(0)
            .getSensitivity());
      }
    };
    return jacobianFunction;
  }

}
