#!/usr/bin/env python3
"""
Async constraint processing for panDecay.

This module provides asynchronous parallel processing for constraint tree
generation and scoring, significantly improving performance for multi-constraint
analyses.
"""

import asyncio
import logging
import time
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Callable, Awaitable
from dataclasses import dataclass
import threading

logger = logging.getLogger(__name__)


@dataclass
class ConstraintTask:
    """Represents a single constraint analysis task."""
    
    branch_num: int
    total_branches: int
    clade_log_idx: int
    clade_taxa_names: List[str]
    tree_idx: int
    
    def __repr__(self) -> str:
        return f"ConstraintTask(branch={self.branch_num}/{self.total_branches}, clade={self.clade_log_idx})"


@dataclass
class ConstraintResult:
    """Represents the result of a constraint analysis."""
    
    task: ConstraintTask
    success: bool
    tree_filename: Optional[str] = None
    likelihood: Optional[float] = None
    error_message: Optional[str] = None
    processing_time: float = 0.0
    
    def __repr__(self) -> str:
        status = "SUCCESS" if self.success else "FAILED"
        return f"ConstraintResult({self.task.clade_log_idx}, {status}, {self.processing_time:.2f}s)"


class ProgressTracker:
    """Thread-safe progress tracking for async constraint processing."""
    
    def __init__(self, total_tasks: int, update_callback: Optional[Callable] = None):
        self.total_tasks = total_tasks
        self.completed_tasks = 0
        self.failed_tasks = 0
        self.start_time = time.time()
        self.update_callback = update_callback
        self._lock = threading.Lock()
        
    def update(self, result: ConstraintResult) -> None:
        """Update progress with a completed task result."""
        with self._lock:
            self.completed_tasks += 1
            if not result.success:
                self.failed_tasks += 1
            
            # Calculate progress statistics
            progress_pct = (self.completed_tasks / self.total_tasks) * 100
            elapsed_time = time.time() - self.start_time
            
            if self.completed_tasks > 0:
                avg_time_per_task = elapsed_time / self.completed_tasks
                eta = avg_time_per_task * (self.total_tasks - self.completed_tasks)
            else:
                eta = 0
            
            # Log progress
            logger.info(
                f"Progress: {self.completed_tasks}/{self.total_tasks} "
                f"({progress_pct:.1f}%) • "
                f"Failed: {self.failed_tasks} • "
                f"ETA: {eta:.0f}s • "
                f"Task: {result.task.clade_log_idx}"
            )
            
            # Call update callback if provided
            if self.update_callback:
                try:
                    self.update_callback(self.completed_tasks, self.total_tasks, result)
                except Exception as e:
                    logger.warning(f"Progress callback failed: {e}")
    
    def get_summary(self) -> Dict[str, Any]:
        """Get progress summary statistics."""
        with self._lock:
            elapsed_time = time.time() - self.start_time
            return {
                'total_tasks': self.total_tasks,
                'completed_tasks': self.completed_tasks,
                'failed_tasks': self.failed_tasks,
                'success_rate': (self.completed_tasks - self.failed_tasks) / max(1, self.completed_tasks),
                'elapsed_time': elapsed_time,
                'avg_time_per_task': elapsed_time / max(1, self.completed_tasks)
            }


class AsyncConstraintProcessor:
    """Asynchronous constraint tree processor with parallel execution."""
    
    def __init__(self, 
                 max_workers: Optional[int] = None,
                 timeout_per_task: int = 600,
                 progress_callback: Optional[Callable] = None):
        """
        Initialize async constraint processor.
        
        Args:
            max_workers: Maximum number of parallel workers (None = auto)
            timeout_per_task: Timeout per constraint task in seconds
            progress_callback: Optional callback for progress updates
        """
        if max_workers is None:
            # Use 75% of available cores, minimum 1, maximum 8
            max_workers = max(1, min(8, int(multiprocessing.cpu_count() * 0.75)))
        
        self.max_workers = max_workers
        self.timeout_per_task = timeout_per_task
        self.progress_callback = progress_callback
        self._semaphore = None
        self._progress_tracker = None
        
        logger.info(f"AsyncConstraintProcessor initialized with {max_workers} workers")
    
    async def process_constraints(self,
                                  testable_branches: List[Tuple[int, Any, int, List[str]]],
                                  constraint_generator: Callable,
                                  use_process_pool: bool = False) -> Tuple[List[ConstraintResult], Dict[str, Any]]:
        """
        Process multiple constraint analyses in parallel.
        
        Args:
            testable_branches: List of testable branch information
            constraint_generator: Function to generate constraint trees
            use_process_pool: Whether to use process pool instead of thread pool
            
        Returns:
            Tuple of (results_list, summary_stats)
        """
        if not testable_branches:
            logger.warning("No testable branches provided for constraint processing")
            return [], {}
        
        # Create constraint tasks
        tasks = self._create_constraint_tasks(testable_branches)
        
        # Initialize progress tracking
        self._progress_tracker = ProgressTracker(
            len(tasks), 
            update_callback=self.progress_callback
        )
        
        # Create semaphore for controlling concurrency
        self._semaphore = asyncio.Semaphore(self.max_workers)
        
        logger.info(
            f"Starting parallel constraint processing: "
            f"{len(tasks)} tasks, {self.max_workers} workers"
        )
        
        start_time = time.time()
        
        try:
            # Execute tasks in parallel
            if use_process_pool:
                results = await self._process_with_process_pool(tasks, constraint_generator)
            else:
                results = await self._process_with_thread_pool(tasks, constraint_generator)
            
            # Calculate final statistics
            processing_time = time.time() - start_time
            summary = self._progress_tracker.get_summary()
            summary.update({
                'total_processing_time': processing_time,
                'tasks_per_second': len(tasks) / processing_time if processing_time > 0 else 0,
                'max_workers': self.max_workers,
                'execution_mode': 'process_pool' if use_process_pool else 'thread_pool'
            })
            
            logger.info(
                f"Constraint processing completed: "
                f"{len(results)} results in {processing_time:.2f}s "
                f"({summary['tasks_per_second']:.2f} tasks/sec)"
            )
            
            return results, summary
            
        except Exception as e:
            logger.error(f"Constraint processing failed: {e}")
            raise
    
    def _create_constraint_tasks(self, testable_branches: List[Tuple[int, Any, int, List[str]]]) -> List[ConstraintTask]:
        """Create ConstraintTask objects from testable branches."""
        tasks = []
        
        for branch_num, (i, clade_obj, clade_log_idx, clade_taxa_names) in enumerate(testable_branches, 1):
            task = ConstraintTask(
                branch_num=branch_num,
                total_branches=len(testable_branches),
                clade_log_idx=clade_log_idx,
                clade_taxa_names=clade_taxa_names,
                tree_idx=clade_log_idx  # Use clade_log_idx as tree_idx
            )
            tasks.append(task)
        
        return tasks
    
    async def _process_with_thread_pool(self,
                                       tasks: List[ConstraintTask],
                                       constraint_generator: Callable) -> List[ConstraintResult]:
        """Process tasks using ThreadPoolExecutor for I/O-bound operations."""
        
        results = []
        
        # Create async tasks for parallel execution
        async_tasks = [
            self._process_single_constraint_async(task, constraint_generator)
            for task in tasks
        ]
        
        # Execute all tasks and collect results as they complete
        for coro in asyncio.as_completed(async_tasks):
            try:
                result = await coro
                results.append(result)
                self._progress_tracker.update(result)
            except Exception as e:
                logger.error(f"Constraint task failed: {e}")
                # Create a failure result
                failed_result = ConstraintResult(
                    task=ConstraintTask(0, 0, -1, [], -1),  # Dummy task
                    success=False,
                    error_message=str(e)
                )
                results.append(failed_result)
                self._progress_tracker.update(failed_result)
        
        return results
    
    async def _process_with_process_pool(self,
                                        tasks: List[ConstraintTask],
                                        constraint_generator: Callable) -> List[ConstraintResult]:
        """Process tasks using ProcessPoolExecutor for CPU-bound operations."""
        
        # Note: ProcessPoolExecutor is more complex due to serialization requirements
        # For now, we'll use ThreadPoolExecutor as PAUP* calls are I/O bound
        logger.info("Process pool requested, but using thread pool due to PAUP* I/O bound nature")
        return await self._process_with_thread_pool(tasks, constraint_generator)
    
    async def _process_single_constraint_async(self,
                                             task: ConstraintTask,
                                             constraint_generator: Callable) -> ConstraintResult:
        """Process a single constraint task asynchronously."""
        
        async with self._semaphore:  # Limit concurrent tasks
            start_time = time.time()
            
            try:
                # Run the constraint generation in a thread pool since it's I/O bound
                loop = asyncio.get_event_loop()
                
                with ThreadPoolExecutor(max_workers=1) as executor:
                    future = loop.run_in_executor(
                        executor,
                        self._execute_constraint_task,
                        task,
                        constraint_generator
                    )
                    
                    # Wait for completion with timeout
                    tree_filename, likelihood = await asyncio.wait_for(
                        future,
                        timeout=self.timeout_per_task
                    )
                
                processing_time = time.time() - start_time
                
                # Create success result
                result = ConstraintResult(
                    task=task,
                    success=True,
                    tree_filename=tree_filename,
                    likelihood=likelihood,
                    processing_time=processing_time
                )
                
                logger.debug(f"Constraint task completed: {result}")
                return result
                
            except asyncio.TimeoutError:
                processing_time = time.time() - start_time
                error_msg = f"Constraint task timed out after {self.timeout_per_task}s"
                logger.warning(f"{error_msg}: {task}")
                
                return ConstraintResult(
                    task=task,
                    success=False,
                    error_message=error_msg,
                    processing_time=processing_time
                )
                
            except Exception as e:
                processing_time = time.time() - start_time
                error_msg = f"Constraint task failed: {e}"
                logger.warning(f"{error_msg}: {task}")
                
                return ConstraintResult(
                    task=task,
                    success=False,
                    error_message=error_msg,
                    processing_time=processing_time
                )
    
    def _execute_constraint_task(self,
                                task: ConstraintTask,
                                constraint_generator: Callable) -> Tuple[Optional[str], Optional[float]]:
        """Execute a single constraint task synchronously."""
        try:
            # Call the constraint generator function
            # This function should match the signature of _generate_and_score_constraint_tree
            tree_filename, likelihood = constraint_generator(
                task.clade_taxa_names,
                task.tree_idx
            )
            
            return tree_filename, likelihood
            
        except Exception as e:
            logger.error(f"Constraint generation failed for task {task}: {e}")
            return None, None
    
    def get_performance_summary(self) -> Dict[str, Any]:
        """Get performance summary for the last processing run."""
        if self._progress_tracker is None:
            return {}
        
        return self._progress_tracker.get_summary()


def create_constraint_processor(config: Optional[Dict[str, Any]] = None) -> AsyncConstraintProcessor:
    """
    Factory function to create AsyncConstraintProcessor with configuration.
    
    Args:
        config: Optional configuration dictionary
        
    Returns:
        Configured AsyncConstraintProcessor instance
    """
    if config is None:
        config = {}
    
    max_workers = config.get('max_workers', None)
    timeout_per_task = config.get('timeout_per_task', 600)
    progress_callback = config.get('progress_callback', None)
    
    return AsyncConstraintProcessor(
        max_workers=max_workers,
        timeout_per_task=timeout_per_task,
        progress_callback=progress_callback
    )


# Backward compatibility function for existing code
async def process_constraints_async(testable_branches: List[Tuple[int, Any, int, List[str]]],
                                   constraint_generator: Callable,
                                   max_workers: Optional[int] = None,
                                   timeout_per_task: int = 600) -> Tuple[List[ConstraintResult], Dict[str, Any]]:
    """
    Convenience function for async constraint processing.
    
    This provides a simple interface for existing code to use async processing
    without needing to manage the AsyncConstraintProcessor directly.
    """
    processor = AsyncConstraintProcessor(
        max_workers=max_workers,
        timeout_per_task=timeout_per_task
    )
    
    return await processor.process_constraints(
        testable_branches=testable_branches,
        constraint_generator=constraint_generator
    )