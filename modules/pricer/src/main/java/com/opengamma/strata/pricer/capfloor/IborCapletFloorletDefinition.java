package com.opengamma.strata.pricer.capfloor;

import com.opengamma.strata.basics.date.BusinessDayAdjustment;
import com.opengamma.strata.basics.date.BusinessDayConventions;
import com.opengamma.strata.basics.date.DayCount;
import com.opengamma.strata.basics.index.IborIndex;

public interface IborCapletFloorletDefinition {

  /**
   * Gets the name of these volatilities.
   * 
   * @return the name
   */
  public abstract IborCapletFloorletVolatilitiesName getName();

  /**
   * Gets the Ibor index for which the data is valid.
   * 
   * @return the Ibor index
   */
  public abstract IborIndex getIndex();

  /**
   * Gets the day count to use.
   * 
   * @return the day count
   */
  public abstract DayCount getDayCount();

  public default BusinessDayAdjustment getBusinessDayAdjustment() {
    return BusinessDayAdjustment.of(BusinessDayConventions.MODIFIED_FOLLOWING, getIndex().getFixingCalendar());
  }

}
