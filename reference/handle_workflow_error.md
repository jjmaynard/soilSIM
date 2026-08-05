# Handle Workflow Error

Comprehensive error handling and recovery for workflow components.

## Usage

``` r
handle_workflow_error(
  error,
  context = "Unknown",
  recovery_action = "stop",
  max_retries = 3
)
```

## Arguments

- error:

  Error object

- context:

  Context information

- recovery_action:

  Recovery action ("stop", "warn", "continue", "retry")

- max_retries:

  Maximum number of retries

## Value

Recovery action result
