package com.opengamma.strata.pricer.capfloor;

import java.time.LocalDate;
import java.time.Period;
import java.time.ZonedDateTime;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.opengamma.strata.basics.ReferenceData;
import com.opengamma.strata.basics.index.IborIndex;
import com.opengamma.strata.collect.ArgChecker;
import com.opengamma.strata.collect.array.DoubleArray;
import com.opengamma.strata.collect.array.DoubleMatrix;
import com.opengamma.strata.market.ValueType;
import com.opengamma.strata.market.curve.Curve;
import com.opengamma.strata.market.curve.CurveMetadata;
import com.opengamma.strata.market.curve.InterpolatedNodalCurve;
import com.opengamma.strata.market.param.CurrencyParameterSensitivities;
import com.opengamma.strata.market.sensitivity.PointSensitivities;
import com.opengamma.strata.market.surface.Surface;
import com.opengamma.strata.market.surface.SurfaceMetadata;
import com.opengamma.strata.math.impl.linearalgebra.DecompositionFactory;
import com.opengamma.strata.math.impl.matrix.MatrixAlgebra;
import com.opengamma.strata.math.impl.matrix.OGMatrixAlgebra;
import com.opengamma.strata.math.impl.minimization.DoubleRangeLimitTransform;
import com.opengamma.strata.math.impl.minimization.NonLinearTransformFunction;
import com.opengamma.strata.math.impl.minimization.ParameterLimitsTransform;
import com.opengamma.strata.math.impl.minimization.ParameterLimitsTransform.LimitType;
import com.opengamma.strata.math.impl.minimization.SingleRangeLimitTransform;
import com.opengamma.strata.math.impl.minimization.UncoupledParameterTransforms;
import com.opengamma.strata.math.impl.statistics.leastsquare.LeastSquareResults;
import com.opengamma.strata.math.impl.statistics.leastsquare.LeastSquareResultsWithTransform;
import com.opengamma.strata.math.impl.statistics.leastsquare.NonLinearLeastSquare;
import com.opengamma.strata.pricer.model.SabrParameters;
import com.opengamma.strata.pricer.model.SabrVolatilityFormula;
import com.opengamma.strata.pricer.option.RawOptionData;
import com.opengamma.strata.pricer.rate.RatesProvider;
import com.opengamma.strata.product.capfloor.ResolvedIborCapFloorLeg;

public class SabrIborCapletFloorletVolatilityBootstrapper extends IborCapletFloorletVolatilityCalibrator {


  public static final SabrIborCapletFloorletVolatilityBootstrapper DEFAULT =
      new SabrIborCapletFloorletVolatilityBootstrapper(
          SabrVolatilityFormula.hagan(),
          VolatilityIborCapFloorLegPricer.DEFAULT,
          ReferenceData.standard());

  private final SabrIborCapFloorLegPricer sabrLegPricer = SabrIborCapFloorLegPricer.DEFAULT;
  private static final MatrixAlgebra MA = new OGMatrixAlgebra();
  private static final NonLinearLeastSquare SOLVER = new NonLinearLeastSquare(DecompositionFactory.SV_COMMONS, MA, 1e-10);

  private static final ParameterLimitsTransform[] DEFAULT_TRANSFORMS;

  private static final double RHO_LIMIT = 0.999;
  static {
    DEFAULT_TRANSFORMS = new ParameterLimitsTransform[4];
    DEFAULT_TRANSFORMS[0] = new SingleRangeLimitTransform(0, LimitType.GREATER_THAN); // alpha > 0
//    DEFAULT_TRANSFORMS[0] = new DoubleRangeLimitTransform(0.001, 3.0);
    DEFAULT_TRANSFORMS[1] = new DoubleRangeLimitTransform(0, 1.0); // 0 <= beta <= 1
    DEFAULT_TRANSFORMS[2] = new DoubleRangeLimitTransform(-RHO_LIMIT, RHO_LIMIT); // -1 <= rho <= 1
    DEFAULT_TRANSFORMS[3] = new SingleRangeLimitTransform(0, LimitType.GREATER_THAN);
//    DEFAULT_TRANSFORMS[3] = new DoubleRangeLimitTransform(0.001, 5.0);
  }

  private final SabrVolatilityFormula sabrVolatilityFormula;

  public static SabrIborCapletFloorletVolatilityBootstrapper of(
      SabrVolatilityFormula sabrVolatilityFormula,
      VolatilityIborCapFloorLegPricer pricer,
      ReferenceData referenceData) {

    return new SabrIborCapletFloorletVolatilityBootstrapper(sabrVolatilityFormula, pricer, referenceData);
  }

  private SabrIborCapletFloorletVolatilityBootstrapper(
      SabrVolatilityFormula sabrVolatilityFormula,
      VolatilityIborCapFloorLegPricer pricer,
      ReferenceData referenceData) {

    super(pricer, referenceData);
    this.sabrVolatilityFormula = ArgChecker.notNull(sabrVolatilityFormula, "sabrVolatilityFormula");
  }

  private DoubleArray computeInitialValues(
      RatesProvider ratesProvider,
      Curve betaCurve,
      List<Double> timeList,
      List<Double> strikeList,
      List<Double> volList,
      List<ResolvedIborCapFloorLeg> capList,
      int[] startIndex,
      int postion,
      boolean betaFixed,
      ValueType valueType) {

    List<Double> vols = volList.subList(startIndex[postion], startIndex[postion + 1]);
    List<Double> strikes = strikeList.subList(startIndex[postion], startIndex[postion + 1]);
    ResolvedIborCapFloorLeg cap = capList.get(startIndex[postion]);
    double factor = valueType.equals(ValueType.BLACK_VOLATILITY) ? 1d
        : 1d / ratesProvider.iborIndexRates(cap.getIndex()).rate(
            cap.getCapletFloorletPeriods().get(cap.getCapletFloorletPeriods().size() - 1).getIborRate().getObservation());
    List<Double> volsEquiv = vols.stream().map(v -> v * factor).collect(Collectors.toList());
    Map<Double, Double> map = IntStream.range(0, strikes.size()).boxed().collect(Collectors.toMap(strikes::get, volsEquiv::get));
    Double minKey =
        map.entrySet().stream().min((entry1, entry2) -> entry1.getValue() > entry2.getValue() ? 1 : -1).get().getKey();
    double nuFirst = 0.1d;
    double betaInitial;
    double alphaInitial = DoubleArray.copyOf(volsEquiv).min();
    if (minKey == strikes.get(0)) {
      betaInitial = betaFixed ? betaCurve.yValue(timeList.get(startIndex[postion])) : 0.95d;
      alphaInitial *= 0.9d;
    } else if (minKey == strikes.get(strikes.size() - 1)) {
      betaInitial = betaFixed ? betaCurve.yValue(timeList.get(startIndex[postion])) : 0.05;
      alphaInitial *= 0.9d;
    } else {
      double grad2 = (volsEquiv.get(strikes.size() - 1) - map.get(minKey)) / (strikes.get(strikes.size() - 1) - minKey);
      double grad1 = (map.get(minKey) - volsEquiv.get(0)) / (minKey - strikes.get(0));
      double cvt = (grad2 - grad1) / (minKey - strikes.get(0));
      nuFirst = cvt * minKey * minKey;
      betaInitial = 0.9d;
    }
    return DoubleArray.of(alphaInitial, betaInitial, -0.5 * betaInitial + 0.5 * (1d - betaInitial), nuFirst);
  }

  @Override
  public SabrIborCapletFloorletVolatilities calibrate(
      IborCapletFloorletDefinition definition,
      ZonedDateTime calibrationDateTime,
      RawOptionData capFloorData,
      RatesProvider ratesProvider) {

    ArgChecker.isTrue(definition instanceof SabrIborCapletFloorletBootstrapDefinition);
    SabrIborCapletFloorletBootstrapDefinition bootstrapDefinition = (SabrIborCapletFloorletBootstrapDefinition) definition;
    IborIndex index = bootstrapDefinition.getIndex();
    LocalDate calibrationDate = calibrationDateTime.toLocalDate();
    LocalDate baseDate = index.getEffectiveDateOffset().adjust(calibrationDate, referenceData);
    LocalDate startDate = baseDate.plus(index.getTenor());
    Function<Surface, IborCapletFloorletVolatilities> volatilitiesFunction = volatilitiesFunction(
        bootstrapDefinition, calibrationDateTime, capFloorData);
    SurfaceMetadata metaData = bootstrapDefinition.createMetadata(capFloorData);;
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

    List<CurveMetadata> metadataList = bootstrapDefinition.createSabrParameterMetadata(capFloorData);
    DoubleArray timeToExpiries = DoubleArray.of(nExpiries, i -> timeList.get(startIndex[i]));
    final BitSet fixed = new BitSet();
    Curve betaCurve;
    if (bootstrapDefinition.getBetaCurve().isPresent()) {
      fixed.set(1); // Beta fixed
      betaCurve = bootstrapDefinition.getBetaCurve().get();
    } else {
      betaCurve = InterpolatedNodalCurve.of(
          metadataList.get(1), timeToExpiries, DoubleArray.filled(nExpiries), bootstrapDefinition.getTimeInterpolator());
    }
    InterpolatedNodalCurve alphaCurve = InterpolatedNodalCurve.of(
        metadataList.get(0), timeToExpiries, DoubleArray.filled(nExpiries), bootstrapDefinition.getTimeInterpolator());
    InterpolatedNodalCurve rhoCurve = InterpolatedNodalCurve.of(
        metadataList.get(2), timeToExpiries, DoubleArray.filled(nExpiries), bootstrapDefinition.getTimeInterpolator());
    InterpolatedNodalCurve nuCurve = InterpolatedNodalCurve.of(
        metadataList.get(3), timeToExpiries, DoubleArray.filled(nExpiries), bootstrapDefinition.getTimeInterpolator());
    SabrParameters sabrParams = bootstrapDefinition.getShiftCurve().isPresent()
        ? SabrParameters.of(
            alphaCurve,
            betaCurve,
            rhoCurve,
            nuCurve,
            bootstrapDefinition.getShiftCurve().get(),
            bootstrapDefinition.getSabrVolatilityFormula())
        : SabrParameters.of(
            alphaCurve,
            betaCurve,
            rhoCurve,
            nuCurve,
            bootstrapDefinition.getSabrVolatilityFormula());
    SabrParametersIborCapletFloorletVolatilities vols =
        SabrParametersIborCapletFloorletVolatilities.of(bootstrapDefinition.getName(), index, calibrationDateTime, sabrParams);
    for (int i = 0; i < nExpiries; ++i) {
      DoubleArray start = computeInitialValues(ratesProvider, betaCurve, timeList, strikeList, volList, capList, startIndex, i,
          fixed.get(1), capFloorData.getDataType());
      UncoupledParameterTransforms transform = new UncoupledParameterTransforms(start, DEFAULT_TRANSFORMS, fixed);
      int size = startIndex[i + 1] - startIndex[i];
      final int timIndex = i;
      final SabrParametersIborCapletFloorletVolatilities finalVols1 = vols;
      final SabrParametersIborCapletFloorletVolatilities finalVols2 = vols;
      Function<DoubleArray, DoubleArray> valueFunction = new Function<DoubleArray, DoubleArray>() {
        @Override
        public DoubleArray apply(DoubleArray x) {
          SabrParametersIborCapletFloorletVolatilities volsNew = fixed.get(1)
              ? finalVols1
                  .withParameter(timIndex, x.get(0))
                  .withParameter(timIndex + nExpiries + betaCurve.getParameterCount(), x.get(2))
                  .withParameter(timIndex + 2 * nExpiries + betaCurve.getParameterCount(), x.get(3))
              : finalVols1
                  .withParameter(timIndex, x.get(0))
                  .withParameter(timIndex + nExpiries, x.get(1))
                  .withParameter(timIndex + 2 * nExpiries, x.get(2))
                  .withParameter(timIndex + 3 * nExpiries, x.get(3));
          return DoubleArray.of(size,
              n -> sabrLegPricer.presentValue(capList.get(startIndex[timIndex] + n), ratesProvider, volsNew).getAmount());
        }
      };
      Function<DoubleArray, DoubleMatrix> jacobianFunction = new Function<DoubleArray, DoubleMatrix>() {
        @Override
        public DoubleMatrix apply(DoubleArray x) {
          SabrParametersIborCapletFloorletVolatilities volsNew = fixed.get(1)
              ? finalVols2
                  .withParameter(timIndex, x.get(0))
                  .withParameter(timIndex + nExpiries + betaCurve.getParameterCount(), x.get(2))
                  .withParameter(timIndex + 2 * nExpiries + betaCurve.getParameterCount(), x.get(3))
              : finalVols2
                  .withParameter(timIndex, x.get(0))
                  .withParameter(timIndex + nExpiries, x.get(1))
                  .withParameter(timIndex + 2 * nExpiries, x.get(2))
                  .withParameter(timIndex + 3 * nExpiries, x.get(3));
          double[][] jacobian = new double[size][4];
          for (int n = 0; n < size; ++n) {
            PointSensitivities point = sabrLegPricer.presentValueSensitivityModelParamsSabr(
                capList.get(startIndex[timIndex] + n), ratesProvider, volsNew).build();
            CurrencyParameterSensitivities sensi = volsNew.parameterSensitivity(point);
            jacobian[n][0] = sensi.getSensitivity(alphaCurve.getName(), index.getCurrency()).getSensitivity().get(timIndex);
            jacobian[n][1] = fixed.get(1) ? 0d
                : sensi.getSensitivity(betaCurve.getName(), index.getCurrency()).getSensitivity().get(timIndex);
            jacobian[n][2] = sensi.getSensitivity(rhoCurve.getName(), index.getCurrency()).getSensitivity().get(timIndex);
            jacobian[n][3] = sensi.getSensitivity(nuCurve.getName(), index.getCurrency()).getSensitivity().get(timIndex);
          }
          return DoubleMatrix.ofUnsafe(jacobian);
        }
      };

      NonLinearTransformFunction transFunc = new NonLinearTransformFunction(valueFunction, jacobianFunction, transform);
      DoubleArray capPrices = DoubleArray.of(size, n -> priceList.get(startIndex[timIndex] + n));
      DoubleArray errors = DoubleArray.filled(size, 1.0);
      LeastSquareResults res = SOLVER.solve(
          capPrices, errors, transFunc.getFittingFunction(), transFunc.getFittingJacobian(), transform.transform(start));
      LeastSquareResultsWithTransform resTransform = new LeastSquareResultsWithTransform(res, transform);
      vols = fixed.get(1)
          ? vols.withParameter(timIndex, resTransform.getModelParameters().get(0))
              .withParameter(timIndex + nExpiries + betaCurve.getParameterCount(), resTransform.getModelParameters().get(2))
              .withParameter(timIndex + 2 * nExpiries + betaCurve.getParameterCount(), resTransform.getModelParameters().get(3))
          : vols.withParameter(timIndex, resTransform.getModelParameters().get(0))
              .withParameter(timIndex + nExpiries, resTransform.getModelParameters().get(1))
              .withParameter(timIndex + 2 * nExpiries, resTransform.getModelParameters().get(2))
              .withParameter(timIndex + 3 * nExpiries, resTransform.getModelParameters().get(3));
//      System.out.println("final chi sq: " + resTransform.getChiSq());
    }
    
//    System.out.println(vols.getParameters().getAlphaCurve());
//    System.out.println(vols.getParameters().getBetaCurve());
//    System.out.println(vols.getParameters().getRhoCurve());
//    System.out.println(vols.getParameters().getNuCurve());
    return vols;
  }

//
//  Function<double[], Double> createFunction() {
//
//  }


}
