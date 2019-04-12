//! Contains core routines and types. The types defined in this module are
//! subject to change and the API should not be relied on.
//!
//! This concerns especially the trait [`PositionStore`](crate::core::PositionStore),
//! which is the trait implemented by the different objects storing sequence
//! record positions in the buffer.
//! All parsers have a generic `PositionStore` parameter, which is useful for
//! implementing a FASTX readers with dynamic dispatch
//! ([`seq_io::fastx::dynamic`](crate::fastx::dynamic)).
//! However, the trait methods of `PositionStore` may not be used directly.
#[macro_use]
mod util;
mod position;
#[macro_use]
mod record;
mod bufreader;
mod inner;
#[macro_use]
mod reader;

pub(crate) use self::inner::*;
pub(crate) use self::util::*;

pub use self::bufreader::*;
pub use self::position::*;
pub(crate) use self::record::*;
