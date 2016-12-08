package com.opengamma.strata.pricer.capfloor;

import static com.opengamma.strata.math.impl.linearalgebra.DecompositionFactory.SV_COMMONS;
import static com.opengamma.strata.math.impl.matrix.MatrixAlgebraFactory.OG_ALGEBRA;

import java.time.LocalDate;
import java.time.Period;
import java.time.ZonedDateTime;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.function.Function;

import com.opengamma.strata.basics.ReferenceData;
import com.opengamma.strata.basics.currency.Currency;
import com.opengamma.strata.basics.index.IborIndex;
import com.opengamma.strata.collect.ArgChecker;
import com.opengamma.strata.collect.array.DoubleArray;
import com.opengamma.strata.collect.array.DoubleMatrix;
import com.opengamma.strata.market.curve.CurveMetadata;
import com.opengamma.strata.market.curve.NodalCurve;
import com.opengamma.strata.market.param.CurrencyParameterSensitivities;
import com.opengamma.strata.market.sensitivity.PointSensitivities;
import com.opengamma.strata.market.surface.Surface;
import com.opengamma.strata.market.surface.SurfaceMetadata;
import com.opengamma.strata.math.impl.minimization.DoubleRangeLimitTransform;
import com.opengamma.strata.math.impl.minimization.NonLinearTransformFunction;
import com.opengamma.strata.math.impl.minimization.ParameterLimitsTransform;
import com.opengamma.strata.math.impl.minimization.UncoupledParameterTransforms;
import com.opengamma.strata.math.impl.statistics.leastsquare.LeastSquareResults;
import com.opengamma.strata.math.impl.statistics.leastsquare.LeastSquareResultsWithTransform;
import com.opengamma.strata.math.impl.statistics.leastsquare.NonLinearLeastSquare;
import com.opengamma.strata.pricer.model.SabrParameters;
import com.opengamma.strata.pricer.option.RawOptionData;
import com.opengamma.strata.pricer.rate.RatesProvider;
import com.opengamma.strata.product.capfloor.ResolvedIborCapFloorLeg;

public class SabrIborCapletFloorletVolatilityCalibrator
    extends IborCapletFloorletVolatilityCalibrator {

  /**
   * Default implementation.
   */
  public static final SabrIborCapletFloorletVolatilityCalibrator DEFAULT = of(
      VolatilityIborCapFloorLegPricer.DEFAULT, SabrIborCapFloorLegPricer.DEFAULT, 1.0e-10, ReferenceData.standard());

  /**
   * Transformation for SABR parameters.
   */
  private static final ParameterLimitsTransform[] TRANSFORMS;
  /**
   * SABR parameter range. 
   */
  private static final double RHO_LIMIT = 0.999;
  static {
    TRANSFORMS = new ParameterLimitsTransform[4];
//    TRANSFORMS[0] = new SingleRangeLimitTransform(0.01, LimitType.GREATER_THAN); // alpha > 0
    TRANSFORMS[0] = new DoubleRangeLimitTransform(0.001, 2d); // alpha > 0
    TRANSFORMS[1] = new DoubleRangeLimitTransform(0.001, 0.999); // 0 <= beta <= 1
    TRANSFORMS[2] = new DoubleRangeLimitTransform(-RHO_LIMIT, RHO_LIMIT); // -1 <= rho <= 1
//    TRANSFORMS[3] = new SingleRangeLimitTransform(0.01, LimitType.GREATER_THAN);
    TRANSFORMS[3] = new DoubleRangeLimitTransform(0.001, 2.5d);
  }

  /**
   * The nonlinear least square solver.
   */
  private final NonLinearLeastSquare solver;
  /**
   * SABR pricer for cap/floor leg.
   */
  private final SabrIborCapFloorLegPricer sabrLegPricer;

  public static SabrIborCapletFloorletVolatilityCalibrator of(
      VolatilityIborCapFloorLegPricer pricer,
      SabrIborCapFloorLegPricer sabrLegPricer,
      double epsilon,
      ReferenceData referenceData) {

    NonLinearLeastSquare solver = new NonLinearLeastSquare(SV_COMMONS, OG_ALGEBRA, epsilon);
    return new SabrIborCapletFloorletVolatilityCalibrator(pricer, sabrLegPricer, solver, referenceData);
  }

  // private constructor
  private SabrIborCapletFloorletVolatilityCalibrator(
      VolatilityIborCapFloorLegPricer pricer,
      SabrIborCapFloorLegPricer sabrLegPricer,
      NonLinearLeastSquare solver,
      ReferenceData referenceData) {

    super(pricer, referenceData);
    this.sabrLegPricer = ArgChecker.notNull(sabrLegPricer, "sabrLegPricer");
    this.solver = ArgChecker.notNull(solver, "solver");
  }

  @Override
  public IborCapletFloorletVolatilityCalibrationResult calibrate(
      IborCapletFloorletDefinition definition,
      ZonedDateTime calibrationDateTime,
      RawOptionData capFloorData,
      RatesProvider ratesProvider) {

    ArgChecker.isTrue(ratesProvider.getValuationDate().equals(calibrationDateTime.toLocalDate()),
        "valuationDate of ratesProvider should be coherent to calibrationDateTime");
    ArgChecker.isTrue(definition instanceof SabrIborCapletFloorletCalibrationDefinition);
    SabrIborCapletFloorletCalibrationDefinition sabrDefinition = (SabrIborCapletFloorletCalibrationDefinition) definition;
    IborIndex index = sabrDefinition.getIndex();
    LocalDate calibrationDate = calibrationDateTime.toLocalDate();
    LocalDate baseDate = index.getEffectiveDateOffset().adjust(calibrationDate, referenceData);
    LocalDate startDate = baseDate.plus(index.getTenor());
    Function<Surface, IborCapletFloorletVolatilities> volatilitiesFunction = volatilitiesFunction(
        sabrDefinition, calibrationDateTime, capFloorData);
    SurfaceMetadata metadata = sabrDefinition.createMetadata(capFloorData);
    List<Period> expiries = capFloorData.getExpiries();
    DoubleArray strikes = capFloorData.getStrikes();
    int nExpiries = expiries.size();
    List<Double> timeList = new ArrayList<>();
    List<Double> strikeList = new ArrayList<>();
    List<Double> volList = new ArrayList<>();
    List<ResolvedIborCapFloorLeg> capList = new ArrayList<>();
    List<Double> priceList = new ArrayList<>();
    List<Double> errorList = new ArrayList<>();
    DoubleMatrix errorMatrix = capFloorData.getError().orElse(DoubleMatrix.filled(nExpiries, strikes.size(), 1d));
    for (int i = 0; i < nExpiries; ++i) {
      LocalDate endDate = baseDate.plus(expiries.get(i));
      DoubleArray volatilityForTime = capFloorData.getData().row(i);
      DoubleArray errorForTime = errorMatrix.row(i);
      reduceRawData(sabrDefinition, ratesProvider, capFloorData.getStrikes(), volatilityForTime, errorForTime, startDate,
          endDate, metadata, volatilitiesFunction, timeList, strikeList, volList, capList, priceList, errorList);
    }

    List<CurveMetadata> metadataList = sabrDefinition.createSabrParameterMetadata();
    DoubleArray initialValues = sabrDefinition.createFullInitialValues();
    List<NodalCurve> curveList = sabrDefinition.createSabrParameterCurve(metadataList, initialValues);
    SabrParameters sabrParams = SabrParameters.of(
        curveList.get(0), curveList.get(1), curveList.get(2), curveList.get(3),
        sabrDefinition.getShiftCurve(), sabrDefinition.getSabrVolatilityFormula());
    SabrParametersIborCapletFloorletVolatilities vols =
        SabrParametersIborCapletFloorletVolatilities.of(sabrDefinition.getName(), index, calibrationDateTime, sabrParams);

    BitSet fixed = new BitSet();
    UncoupledParameterTransforms transform = new UncoupledParameterTransforms(
        initialValues, sabrDefinition.createFullTransform(TRANSFORMS), fixed);
    Function<DoubleArray, DoubleArray> valueFunction = createPriceFunction(
        sabrDefinition, ratesProvider, metadataList, vols, capList);
    Function<DoubleArray, DoubleMatrix> jacobianFunction = createJacobianFunction(
        sabrDefinition, ratesProvider, metadataList, vols, capList, index.getCurrency());
    NonLinearTransformFunction transFunc = new NonLinearTransformFunction(valueFunction, jacobianFunction, transform);
    LeastSquareResults res = solver.solve(
        DoubleArray.copyOf(priceList),
        DoubleArray.copyOf(errorList), transFunc.getFittingFunction(), transFunc.getFittingJacobian(),
        transform.transform(initialValues));
    LeastSquareResultsWithTransform resTransform = new LeastSquareResultsWithTransform(res, transform);
    vols = updateParameters(sabrDefinition, metadataList, vols, resTransform.getModelParameters());
    return IborCapletFloorletVolatilityCalibrationResult.ofLestSquare(vols, res.getChiSq());
  }

  // price function
  private Function<DoubleArray, DoubleArray> createPriceFunction(  // TODO duplicated information
      SabrIborCapletFloorletCalibrationDefinition sabrDefinition,
      RatesProvider ratesProvider,
      List<CurveMetadata> metadataList,
      SabrParametersIborCapletFloorletVolatilities volatilities,
      List<ResolvedIborCapFloorLeg> capList) {

    Function<DoubleArray, DoubleArray> priceFunction = new Function<DoubleArray, DoubleArray>() {
      @Override
      public DoubleArray apply(DoubleArray x) {
        SabrParametersIborCapletFloorletVolatilities volsNew =
            updateParameters(sabrDefinition, metadataList, volatilities, x);
        return DoubleArray.of(capList.size(),
            n -> sabrLegPricer.presentValue(capList.get(n), ratesProvider, volsNew).getAmount());
      }
    };
    return priceFunction;
  }

  // node sensitivity function
  private Function<DoubleArray, DoubleMatrix> createJacobianFunction( // TODO duplicated information
      SabrIborCapletFloorletCalibrationDefinition sabrDefinition,
      RatesProvider ratesProvider,
      List<CurveMetadata> metadataList,
      SabrParametersIborCapletFloorletVolatilities volatilities,
      List<ResolvedIborCapFloorLeg> capList,
      Currency currency) {

    int nCaps = capList.size();
    Function<DoubleArray, DoubleMatrix> jacobianFunction = new Function<DoubleArray, DoubleMatrix>() {
      @Override
      public DoubleMatrix apply(DoubleArray x) {
        SabrParametersIborCapletFloorletVolatilities volsNew =
            updateParameters(sabrDefinition, metadataList, volatilities, x);
        double[][] jacobian = new double[nCaps][];
        for (int i = 0; i < nCaps; ++i) {
          PointSensitivities point = sabrLegPricer.presentValueSensitivityModelParamsSabr(
              capList.get(i), ratesProvider, volsNew).build();
          CurrencyParameterSensitivities sensi = volsNew.parameterSensitivity(point);
          jacobian[i] = sensi.getSensitivity(metadataList.get(0).getCurveName(), currency).getSensitivity()
              .concat(sensi.getSensitivity(metadataList.get(1).getCurveName(), currency).getSensitivity())
              .concat(sensi.getSensitivity(metadataList.get(2).getCurveName(), currency).getSensitivity())
              .concat(sensi.getSensitivity(metadataList.get(3).getCurveName(), currency).getSensitivity())
              .toArray();
        }
        return DoubleMatrix.ofUnsafe(jacobian);
      }
    };
    return jacobianFunction;
  }

  // update vols
  private SabrParametersIborCapletFloorletVolatilities updateParameters( // TODO duplicated information
      SabrIborCapletFloorletCalibrationDefinition sabrDefinition,
      List<CurveMetadata> metadataList,
      SabrParametersIborCapletFloorletVolatilities volatilities,
      DoubleArray newValues) {

    List<NodalCurve> curveList = sabrDefinition.createSabrParameterCurve(metadataList, newValues);
    SabrParameters sabrParams = SabrParameters.of(
        curveList.get(0), curveList.get(1), curveList.get(2), curveList.get(3),
        sabrDefinition.getShiftCurve(), sabrDefinition.getSabrVolatilityFormula());
    SabrParametersIborCapletFloorletVolatilities newVols = SabrParametersIborCapletFloorletVolatilities.of(
        volatilities.getName(), volatilities.getIndex(), volatilities.getValuationDateTime(), sabrParams);
    return newVols;
  }

}
