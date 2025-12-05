module AnythingLLMDocs

import Documenter
using Documenter: RawHTMLHeadContent
using Documenter.HTMLWriter
using Documenter.DocSystem
using HTTP
using JSON
using Markdown
using TOML

const DEFAULT_API_BASE = get(
    ENV,
    "ANYTHINGLLM_API_BASE",
    "https://anythingllm.krastanov.org/api/v1",
)
const DEFAULT_API_KEY = get(
    ENV,
    "ANYTHINGLLM_API_KEY",
    "",
)
const DEFAULT_ALLOWLIST = get(
    ENV,
    "ANYTHINGLLM_ALLOWLIST",
    "",
)

"""
Configuration for connecting to an AnythingLLM instance.
`api_base` should point to the `/api/v1` root, `api_key` is the bearer token,
and `embed_host` is the host used to serve the embed script.
"""
struct AnythingLLMConfig
    api_base::String
    api_key::String
    embed_host::String
end

"""Remove a trailing slash from `s` if present."""
strip_trailing_slash(s::AbstractString) = endswith(s, "/") ? s[1:end-1] : s

"""
Derive the host root from an API base, stripping `/api/v1` and trailing slashes.
"""
function host_from_api(api_base::String)
    stripped = replace(strip_trailing_slash(api_base), r"/api/v1$" => "")
    return strip_trailing_slash(stripped)
end

"""
Load configuration from environment defaults and compute the embed host.
"""
function load_config()
    api_base = strip_trailing_slash(DEFAULT_API_BASE)
    AnythingLLMConfig(api_base, DEFAULT_API_KEY, host_from_api(api_base))
end

"""Build a fully-qualified API URL for `path` relative to the configured base."""
api_url(cfg::AnythingLLMConfig, path::AbstractString) = string(cfg.api_base, path)

"""
Perform an HTTP request against AnythingLLM and parse JSON.

Throws if no API key is configured or if the status is non-2xx. Returns the
parsed JSON object or a `Dict` with `"raw"` when parsing fails.
"""
function request_json(cfg::AnythingLLMConfig, method::String, path::String; body=nothing)
    isempty(cfg.api_key) && error("ANYTHINGLLM_API_KEY must be provided")

    headers = ["Authorization" => "Bearer $(cfg.api_key)"]
    payload = body === nothing ? nothing : JSON.json(body)
    if payload !== nothing
        push!(headers, "Content-Type" => "application/json")
    end
    resp = HTTP.request(method, api_url(cfg, path), headers, payload === nothing ? "" : payload)
    if resp.status ÷ 100 == 2
        text = String(resp.body)
        if isempty(text)
            return Dict{String,Any}()
        end
        try
            return JSON.parse(text)
        catch
            return Dict{String,Any}("raw" => text)
        end
    end
    error("AnythingLLM $(method) $(path) failed with status $(resp.status): $(String(resp.body))")
end

"""Create a URL-friendly slug from an arbitrary package name."""
slugify(name::AbstractString) = replace(lowercase(strip(name)), r"[^a-z0-9]+" => "-")

"""List all workspaces from AnythingLLM."""
function current_workspaces(cfg)
    data = request_json(cfg, "GET", "/workspaces")
    get(data, "workspaces", Vector{Any}())
end

"""Delete a workspace by slug."""
function delete_workspace!(cfg, slug::AbstractString)
    request_json(cfg, "DELETE", "/workspace/$slug")
    return slug
end

"""Return all document entries (flattened) from the AnythingLLM file browser."""
function list_documents(cfg)
    data = request_json(cfg, "GET", "/documents")
    items = get(get(data, "localFiles", Dict{String,Any}()), "items", Vector{Any}())
    docs = Any[]
    stack = collect(items)
    while !isempty(stack)
        entry = pop!(stack)
        if get(entry, "type", "") == "file"
            push!(docs, entry)
        elseif get(entry, "type", "") == "folder"
            append!(stack, get(entry, "items", Any[]))
        end
    end
    return docs
end

"""Delete the provided document `names` from AnythingLLM."""
function delete_documents!(cfg, names::Vector{String})
    isempty(names) && return
    request_json(cfg, "DELETE", "/system/remove-documents"; body=Dict("names" => names))
end

"""
Delete uploaded documentation artifacts associated with a workspace.

Matches on known doc and docstring `docSource` markers containing the workspace
name or slug.
"""
function delete_workspace_documents!(cfg, name::String, slug::String)
    docs = list_documents(cfg)
    targets = String[]
    for doc in docs
        docsource = get(doc, "docSource", "")
        docloc = lowercase(String(get(doc, "location", get(doc, "name", ""))))
        if isempty(docsource)
            continue
        end
        if startswith(docsource, "docs/") || startswith(docsource, "docstring:")
            if occursin(lowercase(slug), lowercase(docsource)) || occursin(lowercase(name), lowercase(docsource)) || occursin(lowercase(slug), docloc) || occursin(lowercase(name), docloc)
                push!(targets, String(get(doc, "location", get(doc, "name", ""))))
            end
        end
    end
    delete_documents!(cfg, targets)
end

"""Ensure a folder exists (best effort)."""
function ensure_folder!(cfg, folder::String)
    try
        request_json(cfg, "POST", "/document/create-folder"; body=Dict("name" => folder))
    catch err
        @debug "create-folder failed (possibly exists)" exception = (err, catch_backtrace())
    end
end

"""
Move uploaded documents into a target folder.

`locations` should be the `location` paths returned by AnythingLLM uploads.
"""
function move_files_to_folder!(cfg, folder::String, locations::Vector{String})
    isempty(locations) && return
    ensure_folder!(cfg, folder)
    files = [
        Dict("from" => loc, "to" => string(folder, "/", split(loc, "/")[end]))
        for loc in locations
    ]
    request_json(cfg, "POST", "/document/move-files"; body=Dict("files" => files))
end

function delete_embeds_for_workspace!(cfg, workspace_slugs::Set{String}, workspace_by_id::AbstractDict)
    embeds = get(request_json(cfg, "GET", "/embed"), "embeds", Vector{Any}())
    for embed in embeds
        ws_id = get(get(embed, "workspace", Dict{Any,Any}()), "id", nothing)
        ws_slug = ws_id === nothing ? nothing : get(workspace_by_id, ws_id, nothing)
        if ws_slug !== nothing && ws_slug in workspace_slugs
            uuid = get(embed, "uuid", nothing)
            uuid === nothing && continue
            request_json(cfg, "DELETE", "/embed/$uuid")
        end
    end
end

"""
Recreate the workspace named `name`, ensuring embeds and uploaded docs are removed first.

Returns the newly created workspace object.
"""
function recreate_workspace!(cfg, name::String)
    workspaces = current_workspaces(cfg)
    target_slug = slugify(name)
    workspace_by_id = Dict(ws["id"] => String(ws["slug"]) for ws in workspaces if haskey(ws, "id") && haskey(ws, "slug"))
    matches = Set{String}(
        String(ws["slug"]) for ws in workspaces if (
            get(ws, "name", "") == name || get(ws, "slug", "") == target_slug
        ) && haskey(ws, "slug")
    )
    if !isempty(matches)
        delete_embeds_for_workspace!(cfg, matches, workspace_by_id)
        for slug in matches
            delete_workspace_documents!(cfg, name, slug)
            delete_workspace!(cfg, slug)
        end
    end
    workspace = request_json(cfg, "POST", "/workspace/new"; body=Dict("name" => name))
    return get(workspace, "workspace", workspace)
end

"""
Upload a text blob as a document to AnythingLLM, attaching metadata.

Returns the locations of the created document artifacts.
"""
function upload_raw!(cfg, workspace_slug::String; title::String, text::String, source::String)
    body = Dict(
        "textContent" => text,
        "addToWorkspaces" => workspace_slug,
        "metadata" => Dict(
            "title" => title,
            "docSource" => source,
        ),
    )
    resp = request_json(cfg, "POST", "/document/raw-text"; body=body)
    docs = get(resp, "documents", Any[])
    return String[
        String(get(doc, "location", get(doc, "name", "")))
        for doc in docs if get(doc, "location", nothing) !== nothing || get(doc, "name", nothing) !== nothing
    ]
end

"""Extract the first Markdown heading text, if present."""
function first_heading(text::AbstractString)
    for line in split(text, '\n')
        m = match(r"^\\s*#\\s+(.*)$", line)
        m !== nothing && return strip(m.captures[1])
    end
    return ""
end

"""
Upload all Markdown sources under `docs_root` into the target workspace.

Returns the locations of uploaded artifacts.
"""
function upload_markdown_sources!(cfg, workspace_slug::String, docs_root::String)
    locations = String[]
    for (root, _, files) in walkdir(docs_root)
        for file in files
            endswith(file, ".md") || continue
            path = joinpath(root, file)
            text = read(path, String)
            title = first_heading(text)
            title = isempty(title) ? file : title
            rel = relpath(path, docs_root)
            append!(locations, upload_raw!(cfg, workspace_slug; title=title, text=text, source="docs/$rel"))
        end
    end
    return locations
end

"""Convert Markdown to a plain text representation for upload."""
markdown_string(md::Markdown.MD) = Markdown.plain(md)

"""
Collect docstrings for the given module into a Dict keyed by fully-qualified name.
"""
function docstrings_for_module(mod::Module)
    entries = Dict{String,String}()
    for name in names(mod; all=false, imported=false)
        binding = DocSystem.binding(mod, name)
        docs = DocSystem.getdocs(binding)
        isempty(docs) && continue
        key = string(mod, ".", name)
        texts = [markdown_string(DocSystem.parsedoc(d)) for d in docs]
        entries[key] = join(texts, "\n\n---\n\n")
    end
    return entries
end

"""
Upload collected docstrings for all provided modules to AnythingLLM.

Returns the locations of uploaded artifacts.
"""
function upload_docstrings!(cfg, workspace_slug::String, modules::Vector{Module})
    seen = Set{String}()
    locations = String[]
    for mod in modules
        for (key, text) in docstrings_for_module(mod)
            if key in seen || isempty(strip(text))
                continue
            end
            push!(seen, key)
            append!(
                locations,
                upload_raw!(cfg, workspace_slug; title=key, text=text, source="docstring:$key"),
            )
        end
    end
    return locations
end

"""Create an embed configuration for the workspace and return its UUID."""
function create_embed!(cfg, workspace_slug::String; allowlist::String=DEFAULT_ALLOWLIST)
    body = Dict(
        "workspace_slug" => workspace_slug,
        "chat_mode" => "query",
        "allowlist_domains" => allowlist,
    )
    data = request_json(cfg, "POST", "/embed/new"; body=body)
    embed = get(data, "embed", data)
    uuid = get(embed, "uuid", nothing)
    uuid === nothing && error("Embed creation did not return a UUID")
    return String(uuid)
end

"""Generate the script tag that loads the AnythingLLM embed widget."""
function embed_script(cfg::AnythingLLMConfig, embed_uuid::String)
    api_base = string(cfg.embed_host, "/api/embed")
    script_src = string(cfg.embed_host, "/embed/anythingllm-chat-widget.min.js")
    return """
<script
  data-embed-id="$embed_uuid"
  data-base-api-url="$api_base"
  src="$script_src">
</script>
"""
end

"""
Full AnythingLLM integration hook: reset workspace, upload docs and docstrings,
create an embed, and return the head asset to inject into Documenter.
"""
function integrate_anythingllm(name::String, modules::Vector{Module}, docs_root::String)
    cfg = load_config()
    try
        root = realpath(joinpath(docs_root, ".."))
        devbranch = Documenter.git_remote_head_branch("deploydocs(devbranch = ...)", root)
        decision = Documenter.deploy_folder(
            Documenter.auto_detect_deploy_system();
            branch = "gh-pages",
            branch_previews = "gh-pages",
            devbranch = devbranch,
            devurl = "dev",
            push_preview = false,
            repo = "github.com/QuantumSavory/QuantumClifford.jl.git",
            repo_previews = nothing,
            deploy_repo = nothing,
            tag_prefix = "",
        )
        versions = Any["stable" => "v^", "v#.#", "dev" => "dev"]
        deploy_subfolder = Documenter.determine_deploy_subfolder(decision, versions)

        if !decision.all_ok
            @info "Skipping AnythingLLM deployment; deploy conditions not met."
            return RawHTMLHeadContent[]
        end

        pkg_version = TOML.parsefile(joinpath(root, "Project.toml"))["version"]
        version_tag = isempty(deploy_subfolder) ? string(pkg_version) : string(pkg_version, "-", deploy_subfolder)
        workspace_title = string(name, " ", version_tag)

        workspace = recreate_workspace!(cfg, workspace_title)
        slug = String(get(workspace, "slug", slugify(workspace_title)))
        locations = String[]
        append!(locations, upload_markdown_sources!(cfg, slug, docs_root))
        append!(locations, upload_docstrings!(cfg, slug, modules))
        move_files_to_folder!(cfg, slug, locations)
        uuid = create_embed!(cfg, slug)
        return [RawHTMLHeadContent(embed_script(cfg, uuid))]
    catch err
        @warn "AnythingLLM integration failed; docs will build without chat" exception=(err, catch_backtrace())
        return RawHTMLHeadContent[]
    end
end

end # module
