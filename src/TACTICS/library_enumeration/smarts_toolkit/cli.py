"""
Command-line interface for SMARTS toolkit
"""

import argparse
import json
from pathlib import Path
from .pattern_validator import SMARTSValidator

def main():
    parser = argparse.ArgumentParser(
        description='SMARTS Pattern Troubleshooting Toolkit'
    )
    
    parser.add_argument(
        'reaction_smarts',
        help='Reaction SMARTS pattern'
    )
    
    parser.add_argument(
        'reagent_files',
        nargs='+',
        help='Reagent files (.smi or .csv)'
    )
    
    parser.add_argument(
        '--validate',
        action='store_true',
        help='Run validation'
    )

    parser.add_argument(
        '--export-compatible',
        help='Export compatible reagents to directory'
    )
    
    parser.add_argument(
        '--test-reactions',
        action='store_true',
        help='Test actual reaction combinations'
    )
    
    parser.add_argument(
        '--sample-size',
        type=int,
        default=100,
        help='Number of combinations to test'
    )
    
    parser.add_argument(
        '--output',
        default='smarts_analysis.json',
        help='Output file for results'
    )
    
    args = parser.parse_args()
    
    # Initialize validator
    validator = SMARTSValidator(args.reaction_smarts, args.reagent_files)
    
    results = {}
    
    if args.validate:
        print("Running validation...")
        validation_result = validator.validate(
            test_reactions=args.test_reactions,
            sample_size=args.sample_size if args.test_reactions else None
        )
        
        print(f"\n=== Validation Results ===")
        for position in range(len(args.reagent_files)):
            coverage = validation_result.coverage_stats.get(position, 0)
            print(f"Position {position}: {coverage:.1f}% coverage")
        
        if validation_result.reaction_success_rate > 0:
            print(f"Reaction success rate: {validation_result.reaction_success_rate:.1f}%")
        
        results['validation'] = {
            'coverage': validation_result.coverage_stats,
            'success_rate': validation_result.reaction_success_rate,
            'errors': validation_result.error_messages,
            'warnings': validation_result.warnings
        }

    if args.export_compatible:
        print(f"\nExporting compatible reagents to {args.export_compatible}...")
        validator.export_compatible_reagents(args.export_compatible)
    
    # Save results to JSON
    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {args.output}")

if __name__ == '__main__':
    main()