# MiraDACSymEngineExt — SymEngine.jl interop (plan A.8, T5.5). SymEngine.jl links its own
# libsymengine (another version), so expressions cross as text: the printer of one side, the
# parser of the other.
module MiraDACSymEngineExt

using MiraDAC: MiraDAC, SymExpr, SDA
using SymEngine: SymEngine, Basic

MiraDAC.SymExpr(b::Basic) = SymExpr(string(b))
SymEngine.Basic(x::SymExpr) = Basic(string(x))
MiraDAC.SDA(b::Basic) = SDA(SymExpr(b))

end
