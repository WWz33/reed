//! Per-record evaluation on raw bytes (§4 eval layer): head → mask → tail
//! → edit plan. Pure functions over the transport's byte views; no I/O.

pub mod edit_plan;
pub mod head;
pub mod tail;
