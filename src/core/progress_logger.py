#!/usr/bin/env python3
"""
Progress logging utility for panDecay.

This module provides dynamic console progress updates using carriage return
to overwrite lines, reducing verbose output while maintaining user feedback.
"""

import sys
import logging
from typing import Optional, List, Dict, Any
from pathlib import Path

logger = logging.getLogger(__name__)


class ProgressLogger:
    """
    Handles dynamic progress logging with overwriting capabilities.
    
    Uses carriage return (\\r) to overwrite progress lines and provide
    clean, concise progress updates instead of verbose per-file logging.
    """
    
    def __init__(self, show_progress: bool = True, verbose: bool = False):
        """
        Initialize progress logger.
        
        Args:
            show_progress: Whether to show dynamic progress updates
            verbose: Whether to show detailed logging (disables progress overwriting)
        """
        self.show_progress = show_progress and not verbose
        self.verbose = verbose
        self.current_line = ""
        self.completed_tasks = []
        
    def progress(self, message: str, current: Optional[int] = None, 
                total: Optional[int] = None, overwrite: bool = True):
        """
        Show progress message with optional counter.
        
        Args:
            message: Progress message to display
            current: Current item number (1-based)
            total: Total number of items
            overwrite: Whether to overwrite the previous line
        """
        if not self.show_progress:
            return
            
        # Build progress message
        if current is not None and total is not None:
            progress_msg = f"{message} [{current}/{total}]"
        else:
            progress_msg = message
            
        # Clear previous line if overwriting
        if overwrite and self.current_line:
            sys.stdout.write('\r' + ' ' * len(self.current_line) + '\r')
            
        # Write new progress message
        sys.stdout.write(progress_msg)
        sys.stdout.flush()
        
        self.current_line = progress_msg if overwrite else ""
        
    def complete(self, final_message: str, count: Optional[int] = None, 
                file_type: Optional[str] = None):
        """
        Complete current progress with final message.
        
        Args:
            final_message: Final completion message
            count: Number of items processed
            file_type: Type of files processed
        """
        if not self.show_progress:
            if self.verbose:
                logger.info(final_message)
            return
            
        # Clear current progress line
        if self.current_line:
            sys.stdout.write('\r' + ' ' * len(self.current_line) + '\r')
            
        # Build completion message
        if count is not None:
            if file_type:
                completion_msg = f"✓ {final_message} ({count} {file_type})"
            else:
                completion_msg = f"✓ {final_message} ({count} items)"
        else:
            completion_msg = f"✓ {final_message}"
            
        # Write completion message
        print(completion_msg)
        self.completed_tasks.append(completion_msg)
        self.current_line = ""
        
    def info(self, message: str):
        """
        Show informational message without progress formatting.
        
        Args:
            message: Information message to display
        """
        # Clear current progress line if needed
        if self.current_line:
            sys.stdout.write('\r' + ' ' * len(self.current_line) + '\r')
            
        print(message)
        self.current_line = ""
        
    def warning(self, message: str):
        """
        Show warning message.
        
        Args:
            message: Warning message to display
        """
        # Clear current progress line if needed
        if self.current_line:
            sys.stdout.write('\r' + ' ' * len(self.current_line) + '\r')
            
        print(f"⚠️  {message}")
        self.current_line = ""
        
    def error(self, message: str):
        """
        Show error message.
        
        Args:
            message: Error message to display
        """
        # Clear current progress line if needed
        if self.current_line:
            sys.stdout.write('\r' + ' ' * len(self.current_line) + '\r')
            
        print(f"❌ {message}")
        self.current_line = ""
        
    def section_header(self, title: str):
        """
        Display a section header.
        
        Args:
            title: Section title to display
        """
        # Clear current progress line if needed
        if self.current_line:
            sys.stdout.write('\r' + ' ' * len(self.current_line) + '\r')
            
        print(f"\n{title}")
        print("-" * len(title))
        self.current_line = ""
        
    def file_summary(self, output_dir: Path, file_structure: Dict[str, List[str]]):
        """
        Display organized summary of created files.
        
        Args:
            output_dir: Base output directory
            file_structure: Dictionary mapping subdirs to file lists
        """
        # Clear current progress line if needed
        if self.current_line:
            sys.stdout.write('\r' + ' ' * len(self.current_line) + '\r')
            
        print(f"\n📁 Output files organized in: {output_dir}")
        
        for subdir, files in file_structure.items():
            if files:
                subdir_path = output_dir / subdir if subdir != "." else output_dir
                print(f"  {subdir}/")
                for file in files[:3]:  # Show first 3 files
                    print(f"    - {file}")
                if len(files) > 3:
                    print(f"    ... and {len(files) - 3} more")
                    
        self.current_line = ""


class FileTracker:
    """
    Tracks files created during analysis for organized output summary.
    """
    
    def __init__(self, output_path: Path, base_name: str):
        """
        Initialize file tracker with organized directory structure.
        
        Args:
            output_path: Parent directory where organized folder will be created
            base_name: Base name for the analysis (e.g., alignment name)
        """
        import datetime
        
        # Create timestamp-based directory name
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        self.output_path = output_path
        self.base_name = base_name
        self.organized_dir_name = f"{timestamp}_panDecay_{base_name}"
        self.base_path = output_path / self.organized_dir_name
        
        # Create the organized directory structure
        self.base_path.mkdir(parents=True, exist_ok=True)
        self.files_by_category = {
            "results": [],
            "trees": [],
            "reports": [],
            "visualizations": [],
            "site_analysis": [],
        }
        
    def get_organized_path(self, category: str, filename: str) -> Path:
        """
        Get organized path for a file in the appropriate subdirectory.
        
        Args:
            category: File category (results, trees, reports, etc.)
            filename: Original filename
            
        Returns:
            Organized path in appropriate subdirectory
        """
        category_dirs = {
            "results": "results",
            "trees": "trees",
            "reports": "reports",
            "visualizations": "visualizations", 
            "site_analysis": "site_analysis"
        }
        
        subdir = category_dirs.get(category, "misc")
        organized_dir = self.base_path / subdir
        organized_dir.mkdir(parents=True, exist_ok=True)
        
        return organized_dir / filename
        
    def track_file(self, category: str, filepath: Path, description: str = ""):
        """
        Track a created file in the appropriate category.
        
        Args:
            category: File category
            filepath: Path to the created file
            description: Optional description of the file
        """
        self.files_by_category[category].append(filepath.name)
        
    def get_summary(self) -> Dict[str, List[str]]:
        """
        Get summary of all tracked files by category.
        
        Returns:
            Dictionary mapping categories to file lists
        """
        return {k: v for k, v in self.files_by_category.items() if v}