using ODEInterface, SciMLTesting

# Allow-listed non-public accesses, kept so that tightening QA does not change
# the released API:
#
# * `Markdown.parse` is the only way to build a `Markdown.MD` from a string
#   assembled at runtime (`@md_str` needs a literal), and the stdlib neither
#   exports it nor declares it public.
# * `Base.WrappedException` is the supertype of `WrappedODEException`. It is
#   undocumented, but it is part of this package's released type hierarchy, so
#   it can only be dropped in a breaking release.
# * `Base.dump`'s 4-argument form is Base's internal recursion hook, not public
#   API. `dump(::IO, ::BandedMatrix, ::Integer, ::Any)` is ambiguous with
#   `Base.dump(::IOContext, ::Any, ::Int, ::Any)`, so `ambiguities` excludes
#   `dump` rather than the method being removed.
run_qa(
    ODEInterface;
    aqua_kwargs = (; ambiguities = (; exclude = [Base.dump])),
    ei_kwargs = (;
        all_qualified_accesses_are_public = (; ignore = (:parse, :WrappedException)),
    )
)
