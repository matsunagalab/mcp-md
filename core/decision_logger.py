"""
Decision Logger - Logs all workflow decisions for reproducibility.

Records:
- User queries
- Tool selections
- Parameter choices
- Rationales
"""

import json
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, List


logger = logging.getLogger(__name__)


class DecisionLogger:
    """Logger for workflow decisions"""
    
    def __init__(self, run_dir: Path):
        """Initialize decision logger
        
        Args:
            run_dir: Directory to store decision logs
        """
        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)
        
        self.log_file = self.run_dir / "decisions.jsonl"
        self.summary_file = self.run_dir / "decisions_summary.json"
        
        # Initialize log file
        if not self.log_file.exists():
            self.log_file.touch()
        
        logger.info(f"Decision logger initialized: {self.log_file}")
    
    def log_decision(
        self,
        step: str,
        tool: str,
        params: Dict[str, Any],
        reason: str,
        metadata: Dict[str, Any] = None
    ):
        """Log a single decision
        
        Args:
            step: Workflow step name
            tool: Tool/function name
            params: Parameters passed to tool
            reason: Rationale for this decision
            metadata: Additional metadata (optional)
        """
        entry = {
            "timestamp": datetime.utcnow().isoformat(),
            "step": step,
            "tool": tool,
            "params": params,
            "reason": reason,
            "metadata": metadata or {}
        }
        
        # Append to JSONL file
        with open(self.log_file, 'a') as f:
            f.write(json.dumps(entry) + '\n')
        
        logger.debug(f"Logged decision: {step} -> {tool}")
    
    def log_error(
        self,
        step: str,
        tool: str,
        error: str,
        params: Dict[str, Any] = None
    ):
        """Log an error
        
        Args:
            step: Workflow step name
            tool: Tool/function name
            error: Error message
            params: Parameters that caused error (optional)
        """
        entry = {
            "timestamp": datetime.utcnow().isoformat(),
            "step": step,
            "tool": tool,
            "params": params or {},
            "error": error,
            "type": "error"
        }
        
        with open(self.log_file, 'a') as f:
            f.write(json.dumps(entry) + '\n')
        
        logger.error(f"Logged error: {step} -> {tool}: {error}")
    
    def get_all_decisions(self) -> List[Dict[str, Any]]:
        """Get all logged decisions
        
        Returns:
            List of decision entries
        """
        decisions = []
        
        if not self.log_file.exists():
            return decisions
        
        with open(self.log_file, 'r') as f:
            for line in f:
                if line.strip():
                    decisions.append(json.loads(line))
        
        return decisions
    
    def generate_summary(self) -> Dict[str, Any]:
        """Generate summary of all decisions
        
        Returns:
            Summary dict with statistics
        """
        decisions = self.get_all_decisions()
        
        summary = {
            "total_decisions": len(decisions),
            "steps": {},
            "tools_used": {},
            "errors": [],
            "timeline": []
        }
        
        for decision in decisions:
            # Count by step
            step = decision.get("step", "unknown")
            if step not in summary["steps"]:
                summary["steps"][step] = 0
            summary["steps"][step] += 1
            
            # Count by tool
            tool = decision.get("tool", "unknown")
            if tool not in summary["tools_used"]:
                summary["tools_used"][tool] = 0
            summary["tools_used"][tool] += 1
            
            # Collect errors
            if decision.get("type") == "error":
                summary["errors"].append({
                    "timestamp": decision["timestamp"],
                    "step": step,
                    "tool": tool,
                    "error": decision.get("error")
                })
            
            # Timeline
            summary["timeline"].append({
                "timestamp": decision["timestamp"],
                "step": step,
                "tool": tool
            })
        
        # Save summary
        with open(self.summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"Generated summary: {self.summary_file}")
        return summary
    
    def export_markdown_report(self, output_file: Path = None) -> str:
        """Export decisions as Markdown report
        
        Args:
            output_file: Output file path (default: run_dir/decisions_report.md)
        
        Returns:
            Path to generated report
        """
        if output_file is None:
            output_file = self.run_dir / "decisions_report.md"
        
        decisions = self.get_all_decisions()
        summary = self.generate_summary()
        
        # Generate Markdown
        md = f"""# Workflow Decisions Report

**Run Directory**: {self.run_dir}
**Total Decisions**: {summary['total_decisions']}
**Errors**: {len(summary['errors'])}

## Summary by Step

"""
        
        for step, count in summary["steps"].items():
            md += f"- **{step}**: {count} decisions\n"
        
        md += "\n## Summary by Tool\n\n"
        
        for tool, count in summary["tools_used"].items():
            md += f"- `{tool}`: {count} uses\n"
        
        if summary["errors"]:
            md += "\n## Errors\n\n"
            for error in summary["errors"]:
                md += f"### {error['timestamp']}\n"
                md += f"- **Step**: {error['step']}\n"
                md += f"- **Tool**: {error['tool']}\n"
                md += f"- **Error**: {error['error']}\n\n"
        
        md += "\n## Timeline\n\n"
        
        for event in summary["timeline"]:
            md += f"- `{event['timestamp']}` - {event['step']} -> {event['tool']}\n"
        
        # Write report
        with open(output_file, 'w') as f:
            f.write(md)
        
        logger.info(f"Exported Markdown report: {output_file}")
        return str(output_file)

