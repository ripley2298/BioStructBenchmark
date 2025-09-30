"""
biostructbenchmark/__main__.py
Main entry point for BioStructBenchmark
"""

import sys
import traceback


def main():
    """Main entry point with comprehensive error handling"""
    try:
        # Import CLI module
        from biostructbenchmark.cli import (
            arg_parser, 
            get_analysis_flags, 
            get_structure_pairs,
            get_version
        )
        from biostructbenchmark.core.application import run_analyses
        from biostructbenchmark.core.metrics import generate_batch_report
        
        # Parse arguments
        args = arg_parser()
        analysis_flags = get_analysis_flags(args)
        
        # Setup logging level
        if args.verbose:
            import logging
            logging.basicConfig(level=logging.DEBUG)
        elif args.quiet:
            import logging
            logging.basicConfig(level=logging.ERROR)
        
        # Print header
        if not args.quiet:
            print(f"BioStructBenchmark v{get_version()}")
            print("=" * 70)
            print(f"Experimental: {args.experimental}")
            print(f"Predicted: {args.predicted}")
            print(f"Output: {args.output}")
            print("=" * 70)
        
        # Get structure pairs to process
        structure_pairs = get_structure_pairs(args)
        
        if not structure_pairs:
            print("Error: No valid structure pairs found", file=sys.stderr)
            return 1
        
        # Track overall results
        all_results = []
        failed_pairs = []
        
        # Process each structure pair
        for i, (exp_path, pred_path) in enumerate(structure_pairs, 1):
            if not args.quiet:
                print(f"\n[{i}/{len(structure_pairs)}] Processing: {exp_path.name} vs {pred_path.name}")
            
            # Create output directory for this pair
            pair_name = f"{exp_path.stem}_vs_{pred_path.stem}"
            pair_output = args.output / pair_name
            pair_output.mkdir(parents=True, exist_ok=True)
            
            try:
                # Run analyses based on flags
                pair_results = run_analyses(
                    exp_path, pred_path, pair_output, 
                    analysis_flags, args
                )
                
                all_results.append({
                    'pair': pair_name,
                    'experimental': str(exp_path),
                    'predicted': str(pred_path),
                    'results': pair_results
                })
                
            except Exception as e:
                print(f"Error processing {pair_name}: {e}", file=sys.stderr)
                if args.verbose:
                    traceback.print_exc()
                failed_pairs.append(pair_name)
                continue
        
        # Generate summary report if multiple pairs
        if len(structure_pairs) > 1:
            generate_batch_report(all_results, failed_pairs, args.output)
        
        # Final summary
        if not args.quiet:
            print("\n" + "=" * 70)
            print("ANALYSIS COMPLETE")
            print(f"Processed: {len(all_results)}/{len(structure_pairs)} pairs")
            if failed_pairs:
                print(f"Failed: {', '.join(failed_pairs)}")
            print(f"Results saved to: {args.output}")
            print("=" * 70)
        
        return 0 if not failed_pairs else 1
        
    except KeyboardInterrupt:
        print("\nAnalysis interrupted by user", file=sys.stderr)
        return 130
    except Exception as e:
        print(f"Fatal error: {e}", file=sys.stderr)
        if '--verbose' in sys.argv or '-v' in sys.argv:
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())