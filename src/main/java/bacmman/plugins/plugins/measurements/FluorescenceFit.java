/* 
 * Copyright (C) 2018 Jean Ollion
 *
 * This File is part of BACMMAN
 *
 * BACMMAN is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BACMMAN is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with BACMMAN.  If not, see <http://www.gnu.org/licenses/>.
 */
package bacmman.plugins.plugins.measurements;

import bacmman.configuration.parameters.BoundedNumberParameter;
import bacmman.configuration.parameters.ObjectClassParameter;
import bacmman.configuration.parameters.Parameter;
import bacmman.configuration.parameters.TextParameter;
import bacmman.data_structure.SegmentedObject;
import bacmman.image.Image;
import bacmman.measurement.MeasurementKey;
import bacmman.measurement.MeasurementKeyObject;
import bacmman.plugins.Measurement;
import bacmman.plugins.Tip;
import bacmman.utils.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Function;
import java.util.function.IntToDoubleFunction;
import java.util.stream.IntStream;

/**
 *
 * @author Jean Ollion
 */
public class FluorescenceFit implements Measurement, Tip {
    private ObjectClassParameter bacteriaClass = new ObjectClassParameter("Bacteria");
    private BoundedNumberParameter fitXRange = new BoundedNumberParameter("Fit X-range", 0, 100, 0, null).setTip("Length of segment (in pixels) on which fit will be performed");
    private ObjectClassParameter fluo = new ObjectClassParameter("Fluorescence Channel").setTip("Fluorescence channel to measure");
    private TextParameter suffix = new TextParameter("Suffix", "", false).setTip("Suffix to add to the measurement column name");
    private Parameter[] parameters = new Parameter[]{bacteriaClass, fluo, fitXRange, suffix};
    boolean verbose;
    
    // implementation of the ToolTip interface (help displayed on hover)
    @Override
    public String getToolTipText() {
        return "Estimation of total fluorescence and background of a cell as described in: Kaiser, Matthias, Florian Jug, Thomas Julou, Siddharth Deshpande, Thomas Pfohl, Olin K. Silander, Gene Myers, and Erik van Nimwegen. 2018. “Monitoring Single-Cell Gene Regulation Under Dynamically Controllable Conditions with Integrated Microfluidics and Software.” Nature Communications 9 (1):212. https://doi.org/10.1038/s41467-017-02505-0.";
    }
    
    // implementation of the Measurement interface
    @Override
    public int getCallObjectClassIdx() {
        return bacteriaClass.getSelectedClassIdx();
    }

    // this method asks if the measurement should be applied on each object of the call structure track or only on the first element (ie trackhead)
    // in this case we apply the measurement to each object   
    @Override
    public boolean callOnlyOnTrackHeads() {
        return false;
    }
    
    // name and object class of output values. 
    // in our case all fitted values are returned and associated to bacteria object (4 measurement per bacteria)
    @Override
    public List<MeasurementKey> getMeasurementKeys() {
        String suf = suffix.getValue();
        int bIdx = bacteriaClass.getSelectedClassIdx();
        return new ArrayList<MeasurementKey>(2) {{
            add(new MeasurementKeyObject("FluoAmplitude"+suf, bIdx));
            add(new MeasurementKeyObject("FluoBackground"+suf, bIdx));
        }};
    }

    @Override
    public void performMeasurement(SegmentedObject bacteria) {
        // get observed fluo distribution
        Image fluoImage = bacteria.getRoot().getRawImage(fluo.getSelectedClassIdx()); // the "root" object corresponds to the whole viewfield. the fit extends beyond the range of the microchannel object in wich the bacteria is located
        int yMin = bacteria.getBounds().yMin();
        int yMax = bacteria.getBounds().yMax();
        int xMin = Math.max(fluoImage.xMin(), (int)(bacteria.getBounds().xMean() - fitXRange.getValue().intValue()/2d)); // enshures xMin is within bounds of image
        int xMax = Math.min(fluoImage.xMax(), (int)(bacteria.getBounds().xMean() + fitXRange.getValue().intValue()/2d)); // enshures xMin is within bounds of image
        int z = bacteria.getBounds().zMin(); // here we are assuming 2D images. if not we would need to also include loop in z direction
        double[] observedFluo = IntStream.rangeClosed(xMin, xMax).mapToDouble(x -> {
            return IntStream.rangeClosed(yMin, yMax)
                    .mapToDouble(y -> fluoImage.getPixelWithOffset(x, y, z)).sum(); // sum of fluorescence for each column
        }).toArray();
        
        // initial parameter values & bounds, as in Kaiser 2018. Might be included as parameters in order to better adapt to other systems (e.g with wider bacteria) ? 
        double muStart = bacteria.getBounds().xMean() - xMin; // middle of peak starts at middle of bacteria object. -xMin because in the fitting function the index is 0-based
        double wStart = 5.5; // in a more general case, should it be a value depending on bacteria width ? 
        double[] muBounds = new double[]{40, 60};
        double[] wBounds = new double[] {3, 12};
        double precision = 1e-3;
        int maxIterations = 500;
        
        double[] fittedParams = fitFluo(observedFluo, muStart, muBounds, wStart, wBounds, maxIterations, precision);
        if (verbose) plot(observedFluo, fittedParams);
        
        // store fitted parameters in measurements
        bacteria.getMeasurements().setValue("FluoAmplitude"+suffix.getValue(), fittedParams[0]);
        bacteria.getMeasurements().setValue("FluoBackground"+suffix.getValue(), fittedParams[1]);
    }

    @Override
    public Parameter[] getParameters() {
        return parameters;
    }
    // funtion to set parameters 
    public FluorescenceFit setObjectClasses(int bacteria, int fluo) {
        this.bacteriaClass.setSelectedClassIdx(bacteria);
        this.fluo.setSelectedClassIdx(fluo);
        return this;
    }
    public FluorescenceFit setXFitSegmentLength(int length) {
        this.fitXRange.setValue(length);
        return this;
    }
    public FluorescenceFit setVerbose(boolean verbose) {
        this.verbose= verbose;
        return this;
    }
    
    // processing functions
    
    /**
     * Fits fluorescence profile {@param c_i} with a Lorenzian distribution using expectation maximization procedure as described in: Kaiser, Matthias, Florian Jug, Thomas Julou, Siddharth Deshpande, Thomas Pfohl, Olin K. Silander, Gene Myers, and Erik van Nimwegen. 2018. “Monitoring Single-Cell Gene Regulation Under Dynamically Controllable Conditions with Integrated Microfluidics and Software.” Nature Communications 9 (1):212. https://doi.org/10.1038/s41467-017-02505-0.
     * ci = noise + B + A / ( (1 + i - mu) / w )^2
     * @param c_i observed fluorescence profile
     * @param wStart starting value for width of distribution
     * @param muBounds range [muMin, muMax] for mu optimization
     * @param muStart starting value for center of peak
     * @param wBounds range [wMin, wMax] for w optimization
     * @param maxIterationNumber maximum iteration number
     * @param precision precision of parameter search
     * @return array with fitted values of  A, B, i_mid, w
     */
    public static double[] fitFluo(double[] c_i, double muStart, double[] muBounds, double wStart, double[] wBounds, int maxIterationNumber, double precision) {
        // step 1
        double cMax = Arrays.stream(c_i).max().getAsDouble();
        double cMin = Arrays.stream(c_i).min().getAsDouble();
        // step 2
        double[] params = new double[] {cMax - cMin, cMin, muStart, wStart};
        double[] ro_i = new double[c_i.length];
        List<Double> diffList = new ArrayList<>(maxIterationNumber);
        int iteration = 0;
        double diff = Double.POSITIVE_INFINITY;
        while (iteration++ < maxIterationNumber && diff > precision) diffList.add(diff = iterate(c_i, ro_i, params, muBounds, wBounds, precision));
        
        //Utils.plotProfile("Convergence", ArrayUtil.toPrimitive(diffList), "iteration number", "sum of parameter relative error");
        return params;
    }
    /**
     * Updates {@param parameters} array after one round of parameter fitting
     * @param c_i observed values
     * @param ro_i estimated values of ro_i. ro_i = 1 / (1 + ( (i - mu) / w )^2 )
     * @param parameters parameter array (0=A, 1=B, 2=iMid (mu), 3=w)
     * @param muBounds range for mu optimisation 
     * @param wBounds range for w optimisation
     * @param precision precision for mi and w search
     */
    private static double iterate(double[] c_i, double[] ro_i, double[] parameters, double[] muBounds, double[] wBounds, double precision) {
        // update ro_i array using previous parameters
        IntStream.range(0, c_i.length).forEach(i -> ro_i[i] = 1/(1+Math.pow((i-parameters[2])/parameters[3], 2)));
        
        // step 3: compute ro
        Function<IntToDoubleFunction, Double> sum = fun -> IntStream.range(0, c_i.length).mapToDouble(fun).sum();
        double ro = Arrays.stream(ro_i).sum();
        
        // step 4: compute new value of B
        double B =  1d/c_i.length * sum.apply(i -> parameters[1] * c_i[i] / (parameters[1] + parameters[0] * ro_i[i]) );
        
        // step 5: compute new value of A
        double A =  1d/ro * sum.apply(i -> parameters[0] * c_i[i] * ro_i[i] / (parameters[1] + parameters[0] * ro_i[i]) );
        
        // step 6
        // find mu using root of derivative function with updated values of A and B
        DoubleBinaryOperator ro_f_i_mu = (i, mu) -> 1/(1+Math.pow((i-mu)/parameters[3], 2));
        DoubleUnaryOperator d_mu = mu -> sum.apply(i-> {
                    double ro_i_ = ro_f_i_mu.applyAsDouble(i, mu);
                    return (i-mu)*ro_i_*ro_i_ * (c_i[i] - (B + A * ro_i_) );
                } );
        double mu = getRootBisection(d_mu, 1000, precision, muBounds);
        
        // step 7 find new w using root of derivative function with updated values of A, B, and mu
        DoubleBinaryOperator ro_f_i_w = (i, w) -> 1/(1+Math.pow((i-mu)/w, 2));
        DoubleUnaryOperator d_w = w -> sum.apply(i-> {
            double ro_i_ = ro_f_i_w.applyAsDouble(i, w);
            return ro_i_ * (1 - ro_i_) * (-1 + c_i[i] / (B + A * ro_i_) );
        });
        double w = getRootBisection(d_w, 1000, precision, wBounds);
        
        // compute parameter relative error to estimate fit precision
        double diff = Math.abs((A-parameters[0])/(A+parameters[0])) + Math.abs((B-parameters[1])/(B+parameters[1])) + Math.abs((mu-parameters[2]) / (mu+parameters[2])) + Math.abs((w-parameters[3]) / (w+parameters[3]));
        
        // update parameters
        parameters[0] = A;
        parameters[1] = B;
        parameters[2] = mu;
        parameters[3] = w;
        //logger.debug("fit: Bck{} Fluo: {}, mu: {}, width: {}, ro: {}", B, A, mu, w, ro);
        return diff;
    }
    private static void plot(String title, DoubleUnaryOperator function, double[] x) {
        double[] values = Arrays.stream(x).map(function).toArray();
        Utils.plotProfile(title, values, x, null, null);
    }
    private static void plot(double[] observed, double[] parameters) {
        double[] estimated_c_i = IntStream.range(0, observed.length).mapToDouble(i -> parameters[1] + parameters[0] / (1 + Math.pow((i - parameters[2])/parameters[3], 2) )).toArray();
        double[] x = IntStream.range(0, observed.length).mapToDouble(i -> i).toArray(); // centered on iMid
        Utils.plotProfile("observed and estimated fluo", estimated_c_i, x, observed, x);
    }
    
    
    /**
     * Finds a root of {@param function} within range {@param bounds} using bisection method
     * @param function
     * @param maxIteration maximum number of iteration
     * @param precision precision of the search 
     * @param bounds range of the search [left bound, right bound]
     * @return a root of {@param function}
     */
    public static double getRootBisection(DoubleUnaryOperator function, int maxIteration, double precision, double[] bounds) {
        double x=bounds[0], xLeft = bounds[0], xRight = bounds[1];
        double error = xRight-xLeft;
        int iter = 0;
        
        double fLeft = function.applyAsDouble(xLeft);
        double fx = fLeft;
        while (Math.abs(error) > precision && iter<maxIteration && fx!=0 ) {
            x = ((xLeft+xRight)/2);
            fx = function.applyAsDouble(x);
            if ((fLeft*fx) < 0) {
                xRight  = x;
                error = xRight-xLeft;
            } else {
                xLeft = x;
                error = xRight-xLeft;
                fLeft = fx;
            }
            ++iter;
        }
        return x;
    }

    

}
