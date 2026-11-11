#ifndef EMI_BDDC_PRECISION_HH
#define EMI_BDDC_PRECISION_HH

#include <stdexcept>

namespace EmiBddc
{
  // EMI uses uint32_t quantized symbols. -1 inherits the common bit setting;
  // explicit directional settings must fit that symbol type.
  inline bool validDirectionalCompressionBits(int bits)
  {
    return bits == -1 || (bits >= 1 && bits <= 32);
  }

  /**
   * Validate the explicit CG (one step) / FCG (1..100 steps) constant-split experiment.
   *
   * The existing transfer represents W*v as W*(v-c*1) + W*(c*1), and
   * quantizes only the first term. Nodal coding reconstructs Q(W*(v-c*1));
   * floating transform-first reconstructs T^-1 Q(T(W*(v-c*1))). Both add
   * W*c afterwards with the wire constant. Reducing transfer error does not establish a
   * linear or positive-definite preconditioner. CG/FCG recurrences and
   * breakdown checks remain unchanged; convergence must be measured.
   *
   * Disabled by default. Existing FGMRES/probe validation is independent.
   */
  template <class Options>
  void validateCgConstantSplitExperiment(Options const& options)
  {
    if (options.bddcCgConstantSplitExperiment < 0 || options.bddcCgConstantSplitExperiment > 1)
      throw std::runtime_error("bddcCgConstantSplitExperiment must be 0 or 1");
    if (!options.bddcCgConstantSplitExperiment)
      return;
    int const rBits = options.bddcRestrictCompressionBits > 0
                   ? options.bddcRestrictCompressionBits : options.bddcCompressionBits;
    int const pBits = options.bddcProlongateCompressionBits > 0
                   ? options.bddcProlongateCompressionBits : options.bddcCompressionBits;
    if (!options.bddc || options.sdc || options.direct || !options.bddcUseCg ||
        options.bddcFgmres || options.bddcPreconditionerProbe ||
        !options.bddcCgOriginalResidual || options.bddcOperatorAudit ||
        options.bddcConditionEstimate || options.bddcConstantCorrection || options.masterWorker ||
        !options.bddcProlongateConstantSplit || !options.bddcCompression ||
        (options.maximumNumberOfTimeSteps < 1 ||
         options.maximumNumberOfTimeSteps > (options.bddcFlexibleCg ? 100 : 1)) ||
        options.bddcRestrictionPasses != 1 ||
        options.bddcRelativeTolerance != 1e-7 || rBits != pBits ||
        (rBits != 4 && rBits != 6 && rBits != 8 && rBits != 12 && rBits != 16 && rBits != 20 && rBits != 24 && rBits != 32) || options.bddcRestrictZeroPreserving ||
        options.bddcProlongateZeroPreserving || options.bddcRestrictGroupSize ||
        options.bddcProlongateGroupSize || options.bddcHistogram ||
        options.bddcCgDiagnostics || options.bddcProlongationDiagnostics)
      throw std::runtime_error("CG/FCG constant-split experiment requires compressed non-SDC CG, "
          "original-residual rho, CG one step or FCG 1..100 steps, one restriction pass, tolerance 1e-7, "
          "uniform 4, 6, 8, 12, 16, 20, 24 or 32-bit affine grids, split enabled and diagnostics/corrections disabled");
  }

  /**
   * Configure optional restriction/prolongation quantization overrides.
   * The default retains the existing common precision in both directions.
   * Mixed precision is a separate numerical configuration used to isolate
   * transfer sensitivity; it must not be labelled as uniform low precision.
   */
  template <class Transfer, class Options>
  void configureTransferPrecision(Transfer& transfer, Options const& options)
  {
    if constexpr (requires { transfer.setQuantizationBits(options.bddcCompressionBits); })
    {
      transfer.setQuantizationBits(options.bddcCompressionBits);
      if constexpr (requires { options.bddcProlongateConstantSplit; })
      {
        if constexpr (requires { transfer.setProlongateConstantSplit(false); })
          transfer.setProlongateConstantSplit(options.bddcProlongateConstantSplit!=0);
        else if (options.bddcProlongateConstantSplit)
          throw std::runtime_error("Kaskade headers lack constant-split prolongation");
      }
      if constexpr (requires { options.bddcProlongateGroupSize; })
      {
        if constexpr (requires { transfer.setProlongateGroupSize(0); })
          transfer.setProlongateGroupSize(options.bddcProlongateGroupSize);
        else if (options.bddcProlongateGroupSize)
          throw std::runtime_error("Kaskade headers lack grouped prolongation");
      }
      if constexpr (requires { options.bddcProlongateZeroPreserving; })
      {
        if constexpr (requires { transfer.setProlongateZeroPreserving(false); })
          transfer.setProlongateZeroPreserving(options.bddcProlongateZeroPreserving!=0);
        else if (options.bddcProlongateZeroPreserving)
          throw std::runtime_error("Kaskade headers lack zero-preserving prolongation");
      }
      if constexpr (requires { options.bddcRestrictGroupSize; })
      {
        if constexpr (requires { transfer.setRestrictGroupSize(0); })
          transfer.setRestrictGroupSize(options.bddcRestrictGroupSize);
        else if (options.bddcRestrictGroupSize)
          throw std::runtime_error("Kaskade headers lack grouped restriction quantization");
      }
      if constexpr (requires { options.bddcRestrictZeroPreserving; })
      {
        if constexpr (requires { transfer.setRestrictZeroPreserving(false); })
          transfer.setRestrictZeroPreserving(options.bddcRestrictZeroPreserving!=0);
        else if (options.bddcRestrictZeroPreserving)
          throw std::runtime_error("Kaskade headers lack the zero-preserving restriction experiment");
      }
      if constexpr (requires {
                      options.bddcRestrictCompressionBits;
                      transfer.setRestrictQuantizationBits(1);
                    })
        if (options.bddcRestrictCompressionBits > 0)
          transfer.setRestrictQuantizationBits(options.bddcRestrictCompressionBits);
      if constexpr (requires {
                      options.bddcProlongateCompressionBits;
                      transfer.setProlongateQuantizationBits(1);
                    })
        if (options.bddcProlongateCompressionBits > 0)
          transfer.setProlongateQuantizationBits(options.bddcProlongateCompressionBits);
    }
  }
}

#endif
