//! Optional tracing setup (feature `logging`).
//!
//! Installs a `tracing-subscriber` formatter filtered by the `KIMIYA_LOG`
//! environment variable (same syntax as `RUST_LOG`; default `warn`).

/// Initialize the global tracing subscriber for kimiya logging.
///
/// The filter is read from `KIMIYA_LOG` (e.g. `KIMIYA_LOG=kimiya=debug`),
/// defaulting to `warn`. Safe to call more than once or from a host
/// application that already installed its own global subscriber: in that
/// case this function does nothing (the existing subscriber is kept).
pub fn init() {
    use tracing_subscriber::EnvFilter;
    let filter = EnvFilter::try_from_env("KIMIYA_LOG").unwrap_or_else(|_| EnvFilter::new("warn"));
    // `try_init` fails only if a global subscriber is already set; keep that one.
    let _ = tracing_subscriber::fmt().with_env_filter(filter).try_init();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn init_twice_does_not_panic() {
        // `.init()` panicked on the second call ("a global default trace
        // dispatcher has already been set").
        init();
        init();
    }
}
