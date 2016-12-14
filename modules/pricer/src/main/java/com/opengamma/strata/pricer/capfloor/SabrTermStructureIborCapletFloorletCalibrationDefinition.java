/**
 * Copyright (C) 2016 - present by OpenGamma Inc. and the OpenGamma group of companies
 *
 * Please see distribution for license.
 */
package com.opengamma.strata.pricer.capfloor;

import static com.opengamma.strata.market.ValueType.BLACK_VOLATILITY;
import static com.opengamma.strata.market.ValueType.NORMAL_VOLATILITY;
import static com.opengamma.strata.market.ValueType.SABR_ALPHA;
import static com.opengamma.strata.market.ValueType.SABR_BETA;
import static com.opengamma.strata.market.ValueType.SABR_NU;
import static com.opengamma.strata.market.ValueType.SABR_RHO;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Optional;
import java.util.Set;

import org.joda.beans.Bean;
import org.joda.beans.BeanDefinition;
import org.joda.beans.ImmutableBean;
import org.joda.beans.ImmutableValidator;
import org.joda.beans.JodaBeanUtils;
import org.joda.beans.MetaProperty;
import org.joda.beans.Property;
import org.joda.beans.PropertyDefinition;
import org.joda.beans.impl.direct.DirectFieldsBeanBuilder;
import org.joda.beans.impl.direct.DirectMetaBean;
import org.joda.beans.impl.direct.DirectMetaProperty;
import org.joda.beans.impl.direct.DirectMetaPropertyMap;

import com.google.common.collect.ImmutableList;
import com.opengamma.strata.basics.date.DayCount;
import com.opengamma.strata.basics.index.IborIndex;
import com.opengamma.strata.collect.ArgChecker;
import com.opengamma.strata.collect.array.DoubleArray;
import com.opengamma.strata.market.curve.ConstantCurve;
import com.opengamma.strata.market.curve.ConstantNodalCurve;
import com.opengamma.strata.market.curve.Curve;
import com.opengamma.strata.market.curve.CurveMetadata;
import com.opengamma.strata.market.curve.Curves;
import com.opengamma.strata.market.curve.InterpolatedNodalCurve;
import com.opengamma.strata.market.curve.interpolator.CurveExtrapolator;
import com.opengamma.strata.market.curve.interpolator.CurveInterpolator;
import com.opengamma.strata.market.surface.SurfaceMetadata;
import com.opengamma.strata.market.surface.Surfaces;
import com.opengamma.strata.math.impl.minimization.ParameterLimitsTransform;
import com.opengamma.strata.pricer.model.SabrVolatilityFormula;
import com.opengamma.strata.pricer.option.RawOptionData;

/**
 * Definition of caplet volatilities calibration.
 * <p>
 * This definition is used with {@link SabrTermStructureIborCapletFloorletVolatilityCalibrator}. 
 * The term structure of SABR model parameters is calibrated to cap volatilities. 
 * The SABR parameters are represented by {@code NodalCurve} and the node positions on the curves are flexible.
 * The resulting volatilities object will be {@link SabrParametersIborCapletFloorletVolatilities}.
 */
@BeanDefinition
public final class SabrTermStructureIborCapletFloorletCalibrationDefinition
    implements IborCapletFloorletDefinition, ImmutableBean, Serializable {

  /**
   * The name of the volatilities.
   */
  @PropertyDefinition(validate = "notNull", overrideGet = true)
  private final IborCapletFloorletVolatilitiesName name;
  /**
   * The Ibor index for which the data is valid.
   */
  @PropertyDefinition(validate = "notNull", overrideGet = true)
  private final IborIndex index;
  /**
   * The day count to use.
   */
  @PropertyDefinition(validate = "notNull", overrideGet = true)
  private final DayCount dayCount;
  /**
   * The beta (elasticity) curve.
   * <p>
   * This represents the beta parameter of SABR model.
   * The x value of the curve is the expiry.
   * <p>
   * The beta will be treated as one of the calibration parameters if this field is not specified.
   */
  @PropertyDefinition(get = "optional")
  private final Curve betaCurve;
  /**
   * The shift curve.
   * <p>
   * This represents the shift parameter of shifted SABR model.
   * The x value of the curve is the expiry.
   * <p>
   * The shift is set to be zero if this field is not specified.
   */
  @PropertyDefinition(validate = "notNull")
  private final Curve shiftCurve;
  /**
   * The nodes of SABR parameter curves.
   * <p>
   * The size of the list must be 4, ordered as alpha, beta, rho and nu. 
   * The second element corresponding to beta will be ignored if beta is fixed, i.e., {@code betaCurve} is present.
   * <p>
   * If the number of nodes is greater than 1, the curve will be created with {@code CurveInterpolator} and 
   * {@code CurveExtrapolator} specified below. Otherwise, {@code ConstantNodalCurve} will be created.
   */
  @PropertyDefinition(validate = "notNull")
  private final ImmutableList<DoubleArray> parameterCurveNodes;
  /**
   * The interpolator for the SABR parameters.
   */
  @PropertyDefinition(validate = "notNull")
  private final CurveInterpolator interpolator;
  /**
   * The left extrapolator for the SABR parameters.
   * <p>
   * The flat extrapolation is used if not specified.
   */
  @PropertyDefinition(validate = "notNull")
  private final CurveExtrapolator extrapolatorLeft;
  /**
   * The right extrapolator for the SABR parameters.
   * <p>
   * The flat extrapolation is used if not specified.
   */
  @PropertyDefinition(validate = "notNull")
  private final CurveExtrapolator extrapolatorRight;
  /**
   * The SABR formula.
   */
  @PropertyDefinition(validate = "notNull")
  private final SabrVolatilityFormula sabrVolatilityFormula;
  /**
   * The initial parameter values used in calibration. 
   * <p>
   * Default values will be used if not specified. 
   * The size of this field must be 4, ordered as alpha, beta, rho and nu. 
   * The second element corresponding to beta will be ignored if beta is fixed, i.e., {@code betaCurve} is present.
   */
  @PropertyDefinition(get = "optional")
  private final DoubleArray initialParameters;

  //-------------------------------------------------------------------------
  /**
   * Obtains an instance with fixed beta and nonzero shift.
   * <p>
   * The beta and shift are constant in time.
   * The default initial values will be used in the calibration.
   * 
   * @param name  the name of volatilities
   * @param index  the Ibor index
   * @param dayCount  the day count
   * @param beta  the beta
   * @param shift  the shift
   * @param alphaCurveNodes  the alpha curve nodes
   * @param rhoCurveNodes  the rho curve nodes
   * @param nuCurveNodes  the nu curve nodes
   * @param interpolator  the interpolator
   * @param extrapolatorLeft  the left extrapolator
   * @param extrapolatorRight  the right extrapolator
   * @param sabrVolatilityFormula  the SABR formula
   * @return the instance
   */
  public static SabrTermStructureIborCapletFloorletCalibrationDefinition of(
      IborCapletFloorletVolatilitiesName name,
      IborIndex index,
      DayCount dayCount,
      double beta,
      double shift,
      DoubleArray alphaCurveNodes,
      DoubleArray rhoCurveNodes,
      DoubleArray nuCurveNodes,
      CurveInterpolator interpolator,
      CurveExtrapolator extrapolatorLeft,
      CurveExtrapolator extrapolatorRight,
      SabrVolatilityFormula sabrVolatilityFormula) {

    ConstantCurve betaCurve = ConstantCurve.of(
        Curves.sabrParameterByExpiry(name.getName() + "-Beta", dayCount, SABR_BETA), beta);
    ConstantCurve shiftCurve = ConstantCurve.of("Shift curve", shift);
    return new SabrTermStructureIborCapletFloorletCalibrationDefinition(
        name,
        index,
        dayCount,
        betaCurve,
        shiftCurve,
        ImmutableList.of(alphaCurveNodes, DoubleArray.of(), rhoCurveNodes, nuCurveNodes),
        interpolator,
        extrapolatorLeft,
        extrapolatorRight,
        sabrVolatilityFormula,
        null);
  }

  /**
   * Obtains an instance with fixed beta and zero shift.
   * <p>
   * The default initial values will be used in the calibration.
   * 
   * @param name  the name of volatilities
   * @param index  the Ibor index
   * @param dayCount  the day count
   * @param beta  the beta
   * @param alphaCurveNodes  the alpha curve nodes
   * @param rhoCurveNodes  the rho curve nodes
   * @param nuCurveNodes  the nu curve nodes
   * @param interpolator  the interpolator
   * @param extrapolatorLeft  the left extrapolator
   * @param extrapolatorRight  the right extrapolator
   * @param sabrVolatilityFormula  the SABR formula
   * @return the instance
   */
  public static SabrTermStructureIborCapletFloorletCalibrationDefinition of(
      IborCapletFloorletVolatilitiesName name,
      IborIndex index,
      DayCount dayCount,
      double beta,
      DoubleArray alphaCurveNodes,
      DoubleArray rhoCurveNodes,
      DoubleArray nuCurveNodes,
      CurveInterpolator interpolator,
      CurveExtrapolator extrapolatorLeft,
      CurveExtrapolator extrapolatorRight,
      SabrVolatilityFormula sabrVolatilityFormula) {

    ConstantCurve betaCurve = ConstantCurve.of(
        Curves.sabrParameterByExpiry(name.getName() + "-Beta", dayCount, SABR_BETA), beta);
    Curve shiftCurve = ConstantCurve.of("Zero shift", 0d);
    return new SabrTermStructureIborCapletFloorletCalibrationDefinition(
        name,
        index,
        dayCount,
        betaCurve,
        shiftCurve,
        ImmutableList.of(alphaCurveNodes, DoubleArray.of(), rhoCurveNodes, nuCurveNodes),
        interpolator,
        extrapolatorLeft,
        extrapolatorRight,
        sabrVolatilityFormula,
        null);
  }

  /**
   * Obtains an instance with zero shift and initial values.
   * 
   * @param name  the name of volatilities
   * @param index  the Ibor index
   * @param dayCount  the day count
   * @param alphaCurveNodes  the alpha curve nodes
   * @param betaCurveNodes  the beta curve nodes
   * @param rhoCurveNodes  the rho curve nodes
   * @param nuCurveNodes  the nu curve nodes
   * @param interpolator  the interpolator
   * @param extrapolatorLeft  the left extrapolator
   * @param extrapolatorRight  the right extrapolator
   * @param sabrVolatilityFormula  the SABR formula
   * @param initialParameters  the initial parameters
   * @return the instance
   */
  public static SabrTermStructureIborCapletFloorletCalibrationDefinition of(
      IborCapletFloorletVolatilitiesName name,
      IborIndex index,
      DayCount dayCount,
      DoubleArray alphaCurveNodes,
      DoubleArray betaCurveNodes,
      DoubleArray rhoCurveNodes,
      DoubleArray nuCurveNodes,
      CurveInterpolator interpolator,
      CurveExtrapolator extrapolatorLeft,
      CurveExtrapolator extrapolatorRight,
      SabrVolatilityFormula sabrVolatilityFormula,
      DoubleArray initialParameters) {

    Curve shiftCurve = ConstantCurve.of("Zero shift", 0d);
    return new SabrTermStructureIborCapletFloorletCalibrationDefinition(
        name,
        index,
        dayCount,
        null,
        shiftCurve,
        ImmutableList.of(alphaCurveNodes, betaCurveNodes, rhoCurveNodes, nuCurveNodes),
        interpolator,
        extrapolatorLeft,
        extrapolatorRight,
        sabrVolatilityFormula,
        initialParameters);
  }

  /**
   * Obtains an instance with fixed beta, nonzero shift and initial values.
   * <p>
   * The beta and shift are constant in time.
   * The default initial values will be used in the calibration.
   * 
   * @param name  the name of volatilities
   * @param index  the Ibor index
   * @param dayCount  the day count
   * @param beta  the beta
   * @param shift  the shift
   * @param alphaCurveNodes  the alpha curve nodes
   * @param rhoCurveNodes  the rho curve nodes
   * @param nuCurveNodes  the nu curve nodes
   * @param interpolator  the interpolator
   * @param extrapolatorLeft  the left extrapolator
   * @param extrapolatorRight  the right extrapolator
   * @param sabrVolatilityFormula  the SABR formula
   * @param initialParameters  the initial parameters
   * @return the instance
   */
  public static SabrTermStructureIborCapletFloorletCalibrationDefinition of(
      IborCapletFloorletVolatilitiesName name,
      IborIndex index,
      DayCount dayCount,
      double beta,
      double shift,
      DoubleArray alphaCurveNodes,
      DoubleArray rhoCurveNodes,
      DoubleArray nuCurveNodes,
      CurveInterpolator interpolator,
      CurveExtrapolator extrapolatorLeft,
      CurveExtrapolator extrapolatorRight,
      SabrVolatilityFormula sabrVolatilityFormula,
      DoubleArray initialParameters) {

    ConstantCurve betaCurve = ConstantCurve.of(
        Curves.sabrParameterByExpiry(name.getName() + "-Beta", dayCount, SABR_BETA), beta);
    ConstantCurve shiftCurve = ConstantCurve.of("Shift curve", shift);
    return new SabrTermStructureIborCapletFloorletCalibrationDefinition(
        name, index, dayCount, betaCurve, shiftCurve,
        ImmutableList.of(alphaCurveNodes, DoubleArray.of(), rhoCurveNodes, nuCurveNodes),
        interpolator, extrapolatorLeft,
        extrapolatorRight, sabrVolatilityFormula,
        initialParameters);
  }

  /**
   * Obtains an instance with fixed beta, zero shift and initial values.
   * <p>
   * The beta and shift are constant in time.
   * The default initial values will be used in the calibration.
   * 
   * @param name  the name of volatilities
   * @param index  the Ibor index
   * @param dayCount  the day count
   * @param beta  the beta
   * @param alphaCurveNodes  the alpha curve nodes
   * @param rhoCurveNodes  the rho curve nodes
   * @param nuCurveNodes  the nu curve nodes
   * @param interpolator  the interpolator
   * @param extrapolatorLeft  the left extrapolator
   * @param extrapolatorRight  the right extrapolator
   * @param sabrVolatilityFormula  the SABR formula
   * @param initialParameters  the initial parameters
   * @return the instance
   */
  public static SabrTermStructureIborCapletFloorletCalibrationDefinition of(
      IborCapletFloorletVolatilitiesName name,
      IborIndex index,
      DayCount dayCount,
      double beta,
      DoubleArray alphaCurveNodes,
      DoubleArray rhoCurveNodes,
      DoubleArray nuCurveNodes,
      CurveInterpolator interpolator,
      CurveExtrapolator extrapolatorLeft,
      CurveExtrapolator extrapolatorRight,
      SabrVolatilityFormula sabrVolatilityFormula,
      DoubleArray initialParameters) {

    ConstantCurve betaCurve = ConstantCurve.of(
        Curves.sabrParameterByExpiry(name.getName() + "-Beta", dayCount, SABR_BETA), beta);
    Curve shiftCurve = ConstantCurve.of("Zero shift", 0d);
    return new SabrTermStructureIborCapletFloorletCalibrationDefinition(
        name, index, dayCount, betaCurve, shiftCurve,
        ImmutableList.of(alphaCurveNodes, DoubleArray.of(), rhoCurveNodes, nuCurveNodes),
        interpolator, extrapolatorLeft,
        extrapolatorRight, sabrVolatilityFormula,
        initialParameters);
  }

  @ImmutableValidator
  private void validate() {
    if (parameterCurveNodes != null) {
      ArgChecker.isTrue(parameterCurveNodes.size() == 4, "The size of parameterCurveNodes must be 4");
      ArgChecker.isFalse(parameterCurveNodes.get(0).isEmpty(), "The alpha curve nodes must not be empty");
      ArgChecker.isFalse(parameterCurveNodes.get(2).isEmpty(), "The rho curve nodes must not be empty");
      ArgChecker.isFalse(parameterCurveNodes.get(3).isEmpty(), "The nu curve nodes must not be empty");
      if (betaCurve != null) {
        ArgChecker.isFalse(parameterCurveNodes.get(1).isEmpty(), "The beta curve nodes must not be empty");
      }
    }
  }

  //-------------------------------------------------------------------------
  @Override
  public SurfaceMetadata createMetadata(RawOptionData capFloorData) {
    SurfaceMetadata metadata;
    if (capFloorData.getDataType().equals(BLACK_VOLATILITY)) {
      metadata = Surfaces.blackVolatilityByExpiryStrike(name.getName(), dayCount);
    } else if (capFloorData.getDataType().equals(NORMAL_VOLATILITY)) {
      metadata = Surfaces.normalVolatilityByExpiryStrike(name.getName(), dayCount);
    } else {
      throw new IllegalArgumentException("Data type not supported");
    }
    return metadata;
  }

  /**
   * Creates curve metadata for SABR parameters.
   * <p>
   * The metadata in the list are order as alpha, beta, rho, then nu.  
   * 
   * @return the curve metadata
   */
  public ImmutableList<CurveMetadata> createSabrParameterMetadata() {
    CurveMetadata alphaMetadata = Curves.sabrParameterByExpiry(name.getName() + "-Alpha", dayCount, SABR_ALPHA);
    CurveMetadata betaMetadata = Curves.sabrParameterByExpiry(name.getName() + "-Beta", dayCount, SABR_BETA);
    CurveMetadata rhoMetadata = Curves.sabrParameterByExpiry(name.getName() + "-Nu", dayCount, SABR_RHO);
    CurveMetadata nuMetadata = Curves.sabrParameterByExpiry(name.getName() + "-Rho", dayCount, SABR_NU);
    return ImmutableList.of(alphaMetadata, betaMetadata, rhoMetadata, nuMetadata);
  }

  /**
   * Creates the parameter curves with parameter node values. 
   * <p>
   * The node values must be combined nodes ordered as alpha, beta (if beta is not fixed), rho, then nu. 
   * The returned curves are ordered in the same manner. 
   * If the beta is fixed, {@code betaCurve} is returned as the second element.
   * 
   * @param metadata  the metadata
   * @param nodeValues  the parameter node values 
   * @return the curves
   */
  public List<Curve> createSabrParameterCurve(List<CurveMetadata> metadata, DoubleArray nodeValues) {
    List<Curve> res = new ArrayList<>();
    int offset = 0;
    for (int i = 0; i < 4; ++i) {
      if (i == 1 && getBetaCurve().isPresent()) {
        res.add(getBetaCurve().get());  // beta fixed
      } else {
        int nNodes = parameterCurveNodes.get(i).size();
        int currentOffset = offset;
        if (nNodes > 1) {
          res.add(InterpolatedNodalCurve.of(
              metadata.get(i),
              parameterCurveNodes.get(i),
              DoubleArray.of(nNodes, n -> nodeValues.get(n + currentOffset)),
              interpolator,
              extrapolatorLeft,
              extrapolatorRight));
        } else {
          res.add(ConstantNodalCurve.of(
              metadata.get(i),
              parameterCurveNodes.get(i).get(0),
              nodeValues.get(currentOffset)));
        }
        offset += nNodes;
      }
    }
    return res;
  }

  /**
   * Creates the transformation definition for all the curve parameters. 
   * <p>
   * The elements in {@code transform} must be ordered as alpha, beta, rho, then nu.
   * 
   * @param transform  the transform
   * @return the full transform
   */
  public ParameterLimitsTransform[] createFullTransform(ParameterLimitsTransform[] transform) {
    ArgChecker.isTrue(transform.length == 4,
        "transform must contain transformation defintion for alpha, beta, rho and nu");
    List<ParameterLimitsTransform> fullTransformList = new ArrayList<>();
    int length = 0;
    for (int i = 0; i < 4; ++i) {
      if (i == 1 && getBetaCurve().isPresent()) {
        // beta fixed
      } else {
        int nNodes = parameterCurveNodes.get(i).size();
        fullTransformList.addAll(Collections.nCopies(nNodes, transform[i]));
        length += nNodes;
      }
    }
    return fullTransformList.toArray(new ParameterLimitsTransform[length]);
  }

  /**
   * Create initial values for all the curve parameters. 
   * <p>
   * Default values are used if {@code initialParameters} is not specified.
   * 
   * @return the initial values
   */
  public DoubleArray createFullInitialValues() {
    double[] initialValues =
        getInitialParameters().isPresent() ? getInitialParameters().get().toArray() : new double[] {0.1, 0.7, -0.2, 0.5};
    List<Double> fullInitialValues = new ArrayList<>();
    for (int i = 0; i < 4; ++i) {
      if (i == 1 && getBetaCurve().isPresent()) {
        // beta fixed
      } else {
        int nNodes = parameterCurveNodes.get(i).size();
        fullInitialValues.addAll(Collections.nCopies(nNodes, initialValues[i]));
      }
    }
    return DoubleArray.copyOf(fullInitialValues);
  }

  //------------------------- AUTOGENERATED START -------------------------
  ///CLOVER:OFF
  /**
   * The meta-bean for {@code SabrTermStructureIborCapletFloorletCalibrationDefinition}.
   * @return the meta-bean, not null
   */
  public static SabrTermStructureIborCapletFloorletCalibrationDefinition.Meta meta() {
    return SabrTermStructureIborCapletFloorletCalibrationDefinition.Meta.INSTANCE;
  }

  static {
    JodaBeanUtils.registerMetaBean(SabrTermStructureIborCapletFloorletCalibrationDefinition.Meta.INSTANCE);
  }

  /**
   * The serialization version id.
   */
  private static final long serialVersionUID = 1L;

  /**
   * Returns a builder used to create an instance of the bean.
   * @return the builder, not null
   */
  public static SabrTermStructureIborCapletFloorletCalibrationDefinition.Builder builder() {
    return new SabrTermStructureIborCapletFloorletCalibrationDefinition.Builder();
  }

  private SabrTermStructureIborCapletFloorletCalibrationDefinition(
      IborCapletFloorletVolatilitiesName name,
      IborIndex index,
      DayCount dayCount,
      Curve betaCurve,
      Curve shiftCurve,
      List<DoubleArray> parameterCurveNodes,
      CurveInterpolator interpolator,
      CurveExtrapolator extrapolatorLeft,
      CurveExtrapolator extrapolatorRight,
      SabrVolatilityFormula sabrVolatilityFormula,
      DoubleArray initialParameters) {
    JodaBeanUtils.notNull(name, "name");
    JodaBeanUtils.notNull(index, "index");
    JodaBeanUtils.notNull(dayCount, "dayCount");
    JodaBeanUtils.notNull(shiftCurve, "shiftCurve");
    JodaBeanUtils.notNull(parameterCurveNodes, "parameterCurveNodes");
    JodaBeanUtils.notNull(interpolator, "interpolator");
    JodaBeanUtils.notNull(extrapolatorLeft, "extrapolatorLeft");
    JodaBeanUtils.notNull(extrapolatorRight, "extrapolatorRight");
    JodaBeanUtils.notNull(sabrVolatilityFormula, "sabrVolatilityFormula");
    this.name = name;
    this.index = index;
    this.dayCount = dayCount;
    this.betaCurve = betaCurve;
    this.shiftCurve = shiftCurve;
    this.parameterCurveNodes = ImmutableList.copyOf(parameterCurveNodes);
    this.interpolator = interpolator;
    this.extrapolatorLeft = extrapolatorLeft;
    this.extrapolatorRight = extrapolatorRight;
    this.sabrVolatilityFormula = sabrVolatilityFormula;
    this.initialParameters = initialParameters;
    validate();
  }

  @Override
  public SabrTermStructureIborCapletFloorletCalibrationDefinition.Meta metaBean() {
    return SabrTermStructureIborCapletFloorletCalibrationDefinition.Meta.INSTANCE;
  }

  @Override
  public <R> Property<R> property(String propertyName) {
    return metaBean().<R>metaProperty(propertyName).createProperty(this);
  }

  @Override
  public Set<String> propertyNames() {
    return metaBean().metaPropertyMap().keySet();
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the name of the volatilities.
   * @return the value of the property, not null
   */
  @Override
  public IborCapletFloorletVolatilitiesName getName() {
    return name;
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the Ibor index for which the data is valid.
   * @return the value of the property, not null
   */
  @Override
  public IborIndex getIndex() {
    return index;
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the day count to use.
   * @return the value of the property, not null
   */
  @Override
  public DayCount getDayCount() {
    return dayCount;
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the beta (elasticity) curve.
   * <p>
   * This represents the beta parameter of SABR model.
   * The x value of the curve is the expiry.
   * <p>
   * The beta will be treated as one of the calibration parameters if this field is not specified.
   * @return the optional value of the property, not null
   */
  public Optional<Curve> getBetaCurve() {
    return Optional.ofNullable(betaCurve);
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the shift curve.
   * <p>
   * This represents the shift parameter of shifted SABR model.
   * The x value of the curve is the expiry.
   * <p>
   * The shift is set to be zero if this field is not specified.
   * @return the value of the property, not null
   */
  public Curve getShiftCurve() {
    return shiftCurve;
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the nodes of SABR parameter curves.
   * <p>
   * The size of the list must be 4, ordered as alpha, beta, rho and nu.
   * The second element corresponding to beta will be ignored if beta is fixed, i.e., {@code betaCurve} is present.
   * <p>
   * If the number of nodes is greater than 1, the curve will be created with {@code CurveInterpolator} and
   * {@code CurveExtrapolator} specified below. Otherwise, {@code ConstantNodalCurve} will be created.
   * @return the value of the property, not null
   */
  public ImmutableList<DoubleArray> getParameterCurveNodes() {
    return parameterCurveNodes;
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the interpolator for the SABR parameters.
   * @return the value of the property, not null
   */
  public CurveInterpolator getInterpolator() {
    return interpolator;
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the left extrapolator for the SABR parameters.
   * <p>
   * The flat extrapolation is used if not specified.
   * @return the value of the property, not null
   */
  public CurveExtrapolator getExtrapolatorLeft() {
    return extrapolatorLeft;
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the right extrapolator for the SABR parameters.
   * <p>
   * The flat extrapolation is used if not specified.
   * @return the value of the property, not null
   */
  public CurveExtrapolator getExtrapolatorRight() {
    return extrapolatorRight;
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the SABR formula.
   * @return the value of the property, not null
   */
  public SabrVolatilityFormula getSabrVolatilityFormula() {
    return sabrVolatilityFormula;
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the initial parameter values used in calibration.
   * <p>
   * Default values will be used if not specified.
   * The size of this field must be 4, ordered as alpha, beta, rho and nu.
   * The second element corresponding to beta will be ignored if beta is fixed, i.e., {@code betaCurve} is present.
   * @return the optional value of the property, not null
   */
  public Optional<DoubleArray> getInitialParameters() {
    return Optional.ofNullable(initialParameters);
  }

  //-----------------------------------------------------------------------
  /**
   * Returns a builder that allows this bean to be mutated.
   * @return the mutable builder, not null
   */
  public Builder toBuilder() {
    return new Builder(this);
  }

  @Override
  public boolean equals(Object obj) {
    if (obj == this) {
      return true;
    }
    if (obj != null && obj.getClass() == this.getClass()) {
      SabrTermStructureIborCapletFloorletCalibrationDefinition other = (SabrTermStructureIborCapletFloorletCalibrationDefinition) obj;
      return JodaBeanUtils.equal(name, other.name) &&
          JodaBeanUtils.equal(index, other.index) &&
          JodaBeanUtils.equal(dayCount, other.dayCount) &&
          JodaBeanUtils.equal(betaCurve, other.betaCurve) &&
          JodaBeanUtils.equal(shiftCurve, other.shiftCurve) &&
          JodaBeanUtils.equal(parameterCurveNodes, other.parameterCurveNodes) &&
          JodaBeanUtils.equal(interpolator, other.interpolator) &&
          JodaBeanUtils.equal(extrapolatorLeft, other.extrapolatorLeft) &&
          JodaBeanUtils.equal(extrapolatorRight, other.extrapolatorRight) &&
          JodaBeanUtils.equal(sabrVolatilityFormula, other.sabrVolatilityFormula) &&
          JodaBeanUtils.equal(initialParameters, other.initialParameters);
    }
    return false;
  }

  @Override
  public int hashCode() {
    int hash = getClass().hashCode();
    hash = hash * 31 + JodaBeanUtils.hashCode(name);
    hash = hash * 31 + JodaBeanUtils.hashCode(index);
    hash = hash * 31 + JodaBeanUtils.hashCode(dayCount);
    hash = hash * 31 + JodaBeanUtils.hashCode(betaCurve);
    hash = hash * 31 + JodaBeanUtils.hashCode(shiftCurve);
    hash = hash * 31 + JodaBeanUtils.hashCode(parameterCurveNodes);
    hash = hash * 31 + JodaBeanUtils.hashCode(interpolator);
    hash = hash * 31 + JodaBeanUtils.hashCode(extrapolatorLeft);
    hash = hash * 31 + JodaBeanUtils.hashCode(extrapolatorRight);
    hash = hash * 31 + JodaBeanUtils.hashCode(sabrVolatilityFormula);
    hash = hash * 31 + JodaBeanUtils.hashCode(initialParameters);
    return hash;
  }

  @Override
  public String toString() {
    StringBuilder buf = new StringBuilder(384);
    buf.append("SabrTermStructureIborCapletFloorletCalibrationDefinition{");
    buf.append("name").append('=').append(name).append(',').append(' ');
    buf.append("index").append('=').append(index).append(',').append(' ');
    buf.append("dayCount").append('=').append(dayCount).append(',').append(' ');
    buf.append("betaCurve").append('=').append(betaCurve).append(',').append(' ');
    buf.append("shiftCurve").append('=').append(shiftCurve).append(',').append(' ');
    buf.append("parameterCurveNodes").append('=').append(parameterCurveNodes).append(',').append(' ');
    buf.append("interpolator").append('=').append(interpolator).append(',').append(' ');
    buf.append("extrapolatorLeft").append('=').append(extrapolatorLeft).append(',').append(' ');
    buf.append("extrapolatorRight").append('=').append(extrapolatorRight).append(',').append(' ');
    buf.append("sabrVolatilityFormula").append('=').append(sabrVolatilityFormula).append(',').append(' ');
    buf.append("initialParameters").append('=').append(JodaBeanUtils.toString(initialParameters));
    buf.append('}');
    return buf.toString();
  }

  //-----------------------------------------------------------------------
  /**
   * The meta-bean for {@code SabrTermStructureIborCapletFloorletCalibrationDefinition}.
   */
  public static final class Meta extends DirectMetaBean {
    /**
     * The singleton instance of the meta-bean.
     */
    static final Meta INSTANCE = new Meta();

    /**
     * The meta-property for the {@code name} property.
     */
    private final MetaProperty<IborCapletFloorletVolatilitiesName> name = DirectMetaProperty.ofImmutable(
        this, "name", SabrTermStructureIborCapletFloorletCalibrationDefinition.class, IborCapletFloorletVolatilitiesName.class);
    /**
     * The meta-property for the {@code index} property.
     */
    private final MetaProperty<IborIndex> index = DirectMetaProperty.ofImmutable(
        this, "index", SabrTermStructureIborCapletFloorletCalibrationDefinition.class, IborIndex.class);
    /**
     * The meta-property for the {@code dayCount} property.
     */
    private final MetaProperty<DayCount> dayCount = DirectMetaProperty.ofImmutable(
        this, "dayCount", SabrTermStructureIborCapletFloorletCalibrationDefinition.class, DayCount.class);
    /**
     * The meta-property for the {@code betaCurve} property.
     */
    private final MetaProperty<Curve> betaCurve = DirectMetaProperty.ofImmutable(
        this, "betaCurve", SabrTermStructureIborCapletFloorletCalibrationDefinition.class, Curve.class);
    /**
     * The meta-property for the {@code shiftCurve} property.
     */
    private final MetaProperty<Curve> shiftCurve = DirectMetaProperty.ofImmutable(
        this, "shiftCurve", SabrTermStructureIborCapletFloorletCalibrationDefinition.class, Curve.class);
    /**
     * The meta-property for the {@code parameterCurveNodes} property.
     */
    @SuppressWarnings({"unchecked", "rawtypes" })
    private final MetaProperty<ImmutableList<DoubleArray>> parameterCurveNodes = DirectMetaProperty.ofImmutable(
        this, "parameterCurveNodes", SabrTermStructureIborCapletFloorletCalibrationDefinition.class, (Class) ImmutableList.class);
    /**
     * The meta-property for the {@code interpolator} property.
     */
    private final MetaProperty<CurveInterpolator> interpolator = DirectMetaProperty.ofImmutable(
        this, "interpolator", SabrTermStructureIborCapletFloorletCalibrationDefinition.class, CurveInterpolator.class);
    /**
     * The meta-property for the {@code extrapolatorLeft} property.
     */
    private final MetaProperty<CurveExtrapolator> extrapolatorLeft = DirectMetaProperty.ofImmutable(
        this, "extrapolatorLeft", SabrTermStructureIborCapletFloorletCalibrationDefinition.class, CurveExtrapolator.class);
    /**
     * The meta-property for the {@code extrapolatorRight} property.
     */
    private final MetaProperty<CurveExtrapolator> extrapolatorRight = DirectMetaProperty.ofImmutable(
        this, "extrapolatorRight", SabrTermStructureIborCapletFloorletCalibrationDefinition.class, CurveExtrapolator.class);
    /**
     * The meta-property for the {@code sabrVolatilityFormula} property.
     */
    private final MetaProperty<SabrVolatilityFormula> sabrVolatilityFormula = DirectMetaProperty.ofImmutable(
        this, "sabrVolatilityFormula", SabrTermStructureIborCapletFloorletCalibrationDefinition.class, SabrVolatilityFormula.class);
    /**
     * The meta-property for the {@code initialParameters} property.
     */
    private final MetaProperty<DoubleArray> initialParameters = DirectMetaProperty.ofImmutable(
        this, "initialParameters", SabrTermStructureIborCapletFloorletCalibrationDefinition.class, DoubleArray.class);
    /**
     * The meta-properties.
     */
    private final Map<String, MetaProperty<?>> metaPropertyMap$ = new DirectMetaPropertyMap(
        this, null,
        "name",
        "index",
        "dayCount",
        "betaCurve",
        "shiftCurve",
        "parameterCurveNodes",
        "interpolator",
        "extrapolatorLeft",
        "extrapolatorRight",
        "sabrVolatilityFormula",
        "initialParameters");

    /**
     * Restricted constructor.
     */
    private Meta() {
    }

    @Override
    protected MetaProperty<?> metaPropertyGet(String propertyName) {
      switch (propertyName.hashCode()) {
        case 3373707:  // name
          return name;
        case 100346066:  // index
          return index;
        case 1905311443:  // dayCount
          return dayCount;
        case 1607020767:  // betaCurve
          return betaCurve;
        case 1908090253:  // shiftCurve
          return shiftCurve;
        case -1431162997:  // parameterCurveNodes
          return parameterCurveNodes;
        case 2096253127:  // interpolator
          return interpolator;
        case 1271703994:  // extrapolatorLeft
          return extrapolatorLeft;
        case 773779145:  // extrapolatorRight
          return extrapolatorRight;
        case -683564541:  // sabrVolatilityFormula
          return sabrVolatilityFormula;
        case 1451864142:  // initialParameters
          return initialParameters;
      }
      return super.metaPropertyGet(propertyName);
    }

    @Override
    public SabrTermStructureIborCapletFloorletCalibrationDefinition.Builder builder() {
      return new SabrTermStructureIborCapletFloorletCalibrationDefinition.Builder();
    }

    @Override
    public Class<? extends SabrTermStructureIborCapletFloorletCalibrationDefinition> beanType() {
      return SabrTermStructureIborCapletFloorletCalibrationDefinition.class;
    }

    @Override
    public Map<String, MetaProperty<?>> metaPropertyMap() {
      return metaPropertyMap$;
    }

    //-----------------------------------------------------------------------
    /**
     * The meta-property for the {@code name} property.
     * @return the meta-property, not null
     */
    public MetaProperty<IborCapletFloorletVolatilitiesName> name() {
      return name;
    }

    /**
     * The meta-property for the {@code index} property.
     * @return the meta-property, not null
     */
    public MetaProperty<IborIndex> index() {
      return index;
    }

    /**
     * The meta-property for the {@code dayCount} property.
     * @return the meta-property, not null
     */
    public MetaProperty<DayCount> dayCount() {
      return dayCount;
    }

    /**
     * The meta-property for the {@code betaCurve} property.
     * @return the meta-property, not null
     */
    public MetaProperty<Curve> betaCurve() {
      return betaCurve;
    }

    /**
     * The meta-property for the {@code shiftCurve} property.
     * @return the meta-property, not null
     */
    public MetaProperty<Curve> shiftCurve() {
      return shiftCurve;
    }

    /**
     * The meta-property for the {@code parameterCurveNodes} property.
     * @return the meta-property, not null
     */
    public MetaProperty<ImmutableList<DoubleArray>> parameterCurveNodes() {
      return parameterCurveNodes;
    }

    /**
     * The meta-property for the {@code interpolator} property.
     * @return the meta-property, not null
     */
    public MetaProperty<CurveInterpolator> interpolator() {
      return interpolator;
    }

    /**
     * The meta-property for the {@code extrapolatorLeft} property.
     * @return the meta-property, not null
     */
    public MetaProperty<CurveExtrapolator> extrapolatorLeft() {
      return extrapolatorLeft;
    }

    /**
     * The meta-property for the {@code extrapolatorRight} property.
     * @return the meta-property, not null
     */
    public MetaProperty<CurveExtrapolator> extrapolatorRight() {
      return extrapolatorRight;
    }

    /**
     * The meta-property for the {@code sabrVolatilityFormula} property.
     * @return the meta-property, not null
     */
    public MetaProperty<SabrVolatilityFormula> sabrVolatilityFormula() {
      return sabrVolatilityFormula;
    }

    /**
     * The meta-property for the {@code initialParameters} property.
     * @return the meta-property, not null
     */
    public MetaProperty<DoubleArray> initialParameters() {
      return initialParameters;
    }

    //-----------------------------------------------------------------------
    @Override
    protected Object propertyGet(Bean bean, String propertyName, boolean quiet) {
      switch (propertyName.hashCode()) {
        case 3373707:  // name
          return ((SabrTermStructureIborCapletFloorletCalibrationDefinition) bean).getName();
        case 100346066:  // index
          return ((SabrTermStructureIborCapletFloorletCalibrationDefinition) bean).getIndex();
        case 1905311443:  // dayCount
          return ((SabrTermStructureIborCapletFloorletCalibrationDefinition) bean).getDayCount();
        case 1607020767:  // betaCurve
          return ((SabrTermStructureIborCapletFloorletCalibrationDefinition) bean).betaCurve;
        case 1908090253:  // shiftCurve
          return ((SabrTermStructureIborCapletFloorletCalibrationDefinition) bean).getShiftCurve();
        case -1431162997:  // parameterCurveNodes
          return ((SabrTermStructureIborCapletFloorletCalibrationDefinition) bean).getParameterCurveNodes();
        case 2096253127:  // interpolator
          return ((SabrTermStructureIborCapletFloorletCalibrationDefinition) bean).getInterpolator();
        case 1271703994:  // extrapolatorLeft
          return ((SabrTermStructureIborCapletFloorletCalibrationDefinition) bean).getExtrapolatorLeft();
        case 773779145:  // extrapolatorRight
          return ((SabrTermStructureIborCapletFloorletCalibrationDefinition) bean).getExtrapolatorRight();
        case -683564541:  // sabrVolatilityFormula
          return ((SabrTermStructureIborCapletFloorletCalibrationDefinition) bean).getSabrVolatilityFormula();
        case 1451864142:  // initialParameters
          return ((SabrTermStructureIborCapletFloorletCalibrationDefinition) bean).initialParameters;
      }
      return super.propertyGet(bean, propertyName, quiet);
    }

    @Override
    protected void propertySet(Bean bean, String propertyName, Object newValue, boolean quiet) {
      metaProperty(propertyName);
      if (quiet) {
        return;
      }
      throw new UnsupportedOperationException("Property cannot be written: " + propertyName);
    }

  }

  //-----------------------------------------------------------------------
  /**
   * The bean-builder for {@code SabrTermStructureIborCapletFloorletCalibrationDefinition}.
   */
  public static final class Builder extends DirectFieldsBeanBuilder<SabrTermStructureIborCapletFloorletCalibrationDefinition> {

    private IborCapletFloorletVolatilitiesName name;
    private IborIndex index;
    private DayCount dayCount;
    private Curve betaCurve;
    private Curve shiftCurve;
    private List<DoubleArray> parameterCurveNodes = ImmutableList.of();
    private CurveInterpolator interpolator;
    private CurveExtrapolator extrapolatorLeft;
    private CurveExtrapolator extrapolatorRight;
    private SabrVolatilityFormula sabrVolatilityFormula;
    private DoubleArray initialParameters;

    /**
     * Restricted constructor.
     */
    private Builder() {
    }

    /**
     * Restricted copy constructor.
     * @param beanToCopy  the bean to copy from, not null
     */
    private Builder(SabrTermStructureIborCapletFloorletCalibrationDefinition beanToCopy) {
      this.name = beanToCopy.getName();
      this.index = beanToCopy.getIndex();
      this.dayCount = beanToCopy.getDayCount();
      this.betaCurve = beanToCopy.betaCurve;
      this.shiftCurve = beanToCopy.getShiftCurve();
      this.parameterCurveNodes = beanToCopy.getParameterCurveNodes();
      this.interpolator = beanToCopy.getInterpolator();
      this.extrapolatorLeft = beanToCopy.getExtrapolatorLeft();
      this.extrapolatorRight = beanToCopy.getExtrapolatorRight();
      this.sabrVolatilityFormula = beanToCopy.getSabrVolatilityFormula();
      this.initialParameters = beanToCopy.initialParameters;
    }

    //-----------------------------------------------------------------------
    @Override
    public Object get(String propertyName) {
      switch (propertyName.hashCode()) {
        case 3373707:  // name
          return name;
        case 100346066:  // index
          return index;
        case 1905311443:  // dayCount
          return dayCount;
        case 1607020767:  // betaCurve
          return betaCurve;
        case 1908090253:  // shiftCurve
          return shiftCurve;
        case -1431162997:  // parameterCurveNodes
          return parameterCurveNodes;
        case 2096253127:  // interpolator
          return interpolator;
        case 1271703994:  // extrapolatorLeft
          return extrapolatorLeft;
        case 773779145:  // extrapolatorRight
          return extrapolatorRight;
        case -683564541:  // sabrVolatilityFormula
          return sabrVolatilityFormula;
        case 1451864142:  // initialParameters
          return initialParameters;
        default:
          throw new NoSuchElementException("Unknown property: " + propertyName);
      }
    }

    @SuppressWarnings("unchecked")
    @Override
    public Builder set(String propertyName, Object newValue) {
      switch (propertyName.hashCode()) {
        case 3373707:  // name
          this.name = (IborCapletFloorletVolatilitiesName) newValue;
          break;
        case 100346066:  // index
          this.index = (IborIndex) newValue;
          break;
        case 1905311443:  // dayCount
          this.dayCount = (DayCount) newValue;
          break;
        case 1607020767:  // betaCurve
          this.betaCurve = (Curve) newValue;
          break;
        case 1908090253:  // shiftCurve
          this.shiftCurve = (Curve) newValue;
          break;
        case -1431162997:  // parameterCurveNodes
          this.parameterCurveNodes = (List<DoubleArray>) newValue;
          break;
        case 2096253127:  // interpolator
          this.interpolator = (CurveInterpolator) newValue;
          break;
        case 1271703994:  // extrapolatorLeft
          this.extrapolatorLeft = (CurveExtrapolator) newValue;
          break;
        case 773779145:  // extrapolatorRight
          this.extrapolatorRight = (CurveExtrapolator) newValue;
          break;
        case -683564541:  // sabrVolatilityFormula
          this.sabrVolatilityFormula = (SabrVolatilityFormula) newValue;
          break;
        case 1451864142:  // initialParameters
          this.initialParameters = (DoubleArray) newValue;
          break;
        default:
          throw new NoSuchElementException("Unknown property: " + propertyName);
      }
      return this;
    }

    @Override
    public Builder set(MetaProperty<?> property, Object value) {
      super.set(property, value);
      return this;
    }

    @Override
    public Builder setString(String propertyName, String value) {
      setString(meta().metaProperty(propertyName), value);
      return this;
    }

    @Override
    public Builder setString(MetaProperty<?> property, String value) {
      super.setString(property, value);
      return this;
    }

    @Override
    public Builder setAll(Map<String, ? extends Object> propertyValueMap) {
      super.setAll(propertyValueMap);
      return this;
    }

    @Override
    public SabrTermStructureIborCapletFloorletCalibrationDefinition build() {
      return new SabrTermStructureIborCapletFloorletCalibrationDefinition(
          name,
          index,
          dayCount,
          betaCurve,
          shiftCurve,
          parameterCurveNodes,
          interpolator,
          extrapolatorLeft,
          extrapolatorRight,
          sabrVolatilityFormula,
          initialParameters);
    }

    //-----------------------------------------------------------------------
    /**
     * Sets the name of the volatilities.
     * @param name  the new value, not null
     * @return this, for chaining, not null
     */
    public Builder name(IborCapletFloorletVolatilitiesName name) {
      JodaBeanUtils.notNull(name, "name");
      this.name = name;
      return this;
    }

    /**
     * Sets the Ibor index for which the data is valid.
     * @param index  the new value, not null
     * @return this, for chaining, not null
     */
    public Builder index(IborIndex index) {
      JodaBeanUtils.notNull(index, "index");
      this.index = index;
      return this;
    }

    /**
     * Sets the day count to use.
     * @param dayCount  the new value, not null
     * @return this, for chaining, not null
     */
    public Builder dayCount(DayCount dayCount) {
      JodaBeanUtils.notNull(dayCount, "dayCount");
      this.dayCount = dayCount;
      return this;
    }

    /**
     * Sets the beta (elasticity) curve.
     * <p>
     * This represents the beta parameter of SABR model.
     * The x value of the curve is the expiry.
     * <p>
     * The beta will be treated as one of the calibration parameters if this field is not specified.
     * @param betaCurve  the new value
     * @return this, for chaining, not null
     */
    public Builder betaCurve(Curve betaCurve) {
      this.betaCurve = betaCurve;
      return this;
    }

    /**
     * Sets the shift curve.
     * <p>
     * This represents the shift parameter of shifted SABR model.
     * The x value of the curve is the expiry.
     * <p>
     * The shift is set to be zero if this field is not specified.
     * @param shiftCurve  the new value, not null
     * @return this, for chaining, not null
     */
    public Builder shiftCurve(Curve shiftCurve) {
      JodaBeanUtils.notNull(shiftCurve, "shiftCurve");
      this.shiftCurve = shiftCurve;
      return this;
    }

    /**
     * Sets the nodes of SABR parameter curves.
     * <p>
     * The size of the list must be 4, ordered as alpha, beta, rho and nu.
     * The second element corresponding to beta will be ignored if beta is fixed, i.e., {@code betaCurve} is present.
     * <p>
     * If the number of nodes is greater than 1, the curve will be created with {@code CurveInterpolator} and
     * {@code CurveExtrapolator} specified below. Otherwise, {@code ConstantNodalCurve} will be created.
     * @param parameterCurveNodes  the new value, not null
     * @return this, for chaining, not null
     */
    public Builder parameterCurveNodes(List<DoubleArray> parameterCurveNodes) {
      JodaBeanUtils.notNull(parameterCurveNodes, "parameterCurveNodes");
      this.parameterCurveNodes = parameterCurveNodes;
      return this;
    }

    /**
     * Sets the {@code parameterCurveNodes} property in the builder
     * from an array of objects.
     * @param parameterCurveNodes  the new value, not null
     * @return this, for chaining, not null
     */
    public Builder parameterCurveNodes(DoubleArray... parameterCurveNodes) {
      return parameterCurveNodes(ImmutableList.copyOf(parameterCurveNodes));
    }

    /**
     * Sets the interpolator for the SABR parameters.
     * @param interpolator  the new value, not null
     * @return this, for chaining, not null
     */
    public Builder interpolator(CurveInterpolator interpolator) {
      JodaBeanUtils.notNull(interpolator, "interpolator");
      this.interpolator = interpolator;
      return this;
    }

    /**
     * Sets the left extrapolator for the SABR parameters.
     * <p>
     * The flat extrapolation is used if not specified.
     * @param extrapolatorLeft  the new value, not null
     * @return this, for chaining, not null
     */
    public Builder extrapolatorLeft(CurveExtrapolator extrapolatorLeft) {
      JodaBeanUtils.notNull(extrapolatorLeft, "extrapolatorLeft");
      this.extrapolatorLeft = extrapolatorLeft;
      return this;
    }

    /**
     * Sets the right extrapolator for the SABR parameters.
     * <p>
     * The flat extrapolation is used if not specified.
     * @param extrapolatorRight  the new value, not null
     * @return this, for chaining, not null
     */
    public Builder extrapolatorRight(CurveExtrapolator extrapolatorRight) {
      JodaBeanUtils.notNull(extrapolatorRight, "extrapolatorRight");
      this.extrapolatorRight = extrapolatorRight;
      return this;
    }

    /**
     * Sets the SABR formula.
     * @param sabrVolatilityFormula  the new value, not null
     * @return this, for chaining, not null
     */
    public Builder sabrVolatilityFormula(SabrVolatilityFormula sabrVolatilityFormula) {
      JodaBeanUtils.notNull(sabrVolatilityFormula, "sabrVolatilityFormula");
      this.sabrVolatilityFormula = sabrVolatilityFormula;
      return this;
    }

    /**
     * Sets the initial parameter values used in calibration.
     * <p>
     * Default values will be used if not specified.
     * The size of this field must be 4, ordered as alpha, beta, rho and nu.
     * The second element corresponding to beta will be ignored if beta is fixed, i.e., {@code betaCurve} is present.
     * @param initialParameters  the new value
     * @return this, for chaining, not null
     */
    public Builder initialParameters(DoubleArray initialParameters) {
      this.initialParameters = initialParameters;
      return this;
    }

    //-----------------------------------------------------------------------
    @Override
    public String toString() {
      StringBuilder buf = new StringBuilder(384);
      buf.append("SabrTermStructureIborCapletFloorletCalibrationDefinition.Builder{");
      buf.append("name").append('=').append(JodaBeanUtils.toString(name)).append(',').append(' ');
      buf.append("index").append('=').append(JodaBeanUtils.toString(index)).append(',').append(' ');
      buf.append("dayCount").append('=').append(JodaBeanUtils.toString(dayCount)).append(',').append(' ');
      buf.append("betaCurve").append('=').append(JodaBeanUtils.toString(betaCurve)).append(',').append(' ');
      buf.append("shiftCurve").append('=').append(JodaBeanUtils.toString(shiftCurve)).append(',').append(' ');
      buf.append("parameterCurveNodes").append('=').append(JodaBeanUtils.toString(parameterCurveNodes)).append(',').append(' ');
      buf.append("interpolator").append('=').append(JodaBeanUtils.toString(interpolator)).append(',').append(' ');
      buf.append("extrapolatorLeft").append('=').append(JodaBeanUtils.toString(extrapolatorLeft)).append(',').append(' ');
      buf.append("extrapolatorRight").append('=').append(JodaBeanUtils.toString(extrapolatorRight)).append(',').append(' ');
      buf.append("sabrVolatilityFormula").append('=').append(JodaBeanUtils.toString(sabrVolatilityFormula)).append(',').append(' ');
      buf.append("initialParameters").append('=').append(JodaBeanUtils.toString(initialParameters));
      buf.append('}');
      return buf.toString();
    }

  }

  ///CLOVER:ON
  //-------------------------- AUTOGENERATED END --------------------------
}
