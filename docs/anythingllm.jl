using HTTP
using JSON
using Documenter

"""
    configure_anythingllm(package_name::String;
                          api_url::String="https://anythingllm.krastanov.org",
                          api_key::String=get(ENV, "ANYTHINGLLM_API_KEY", ""),
                          modules::Vector=[])

Configure AnythingLLM chat for the documentation.
This function will:
1. Delete any existing workspace with the same name
2. Create a new workspace
3. Upload all markdown documentation files
4. Upload all public docstrings
5. Create or get an embed for the workspace
6. Return the embed script HTML

# Arguments
- `package_name`: Name of the package (used for workspace name)
- `api_url`: Base URL of the AnythingLLM instance
- `api_key`: API key for authentication (defaults to ANYTHINGLLM_API_KEY env var)
- `modules`: Vector of modules to extract docstrings from (default: empty)

# Returns
- HTML string containing the embed script tag
"""
function configure_anythingllm(package_name::String;
                               api_url::String="https://anythingllm.krastanov.org",
                               api_key::String=get(ENV, "ANYTHINGLLM_API_KEY", ""),
                               modules::Vector=[])

    if isempty(api_key)
        @warn "AnythingLLM API key not provided. Skipping LLM chat configuration."
        return ""
    end

    headers = Dict(
        "Authorization" => "Bearer $api_key",
        "Content-Type" => "application/json"
    )

    workspace_slug = lowercase(package_name)

    @info "Configuring AnythingLLM for $package_name..."

    # Step 1: Delete existing workspace and embed (if any)
    @info "Cleaning up existing workspace and embed..."
    try
        response = HTTP.get("$api_url/api/v1/workspaces", headers=headers)
        workspaces = JSON.parse(String(response.body))

        for ws in get(workspaces, "workspaces", [])
            if lowercase(get(ws, "slug", "")) == workspace_slug
                @info "Deleting existing workspace: $(ws["slug"])"

                # First, get all documents in the workspace and delete them
                try
                    ws_response = HTTP.get("$api_url/api/v1/workspace/$(ws["slug"])", headers=headers)
                    ws_data = JSON.parse(String(ws_response.body))
                    ws_info = ws_data["workspace"]

                    # If ws_info is an array, take the first element
                    if isa(ws_info, Vector)
                        ws_info = ws_info[1]
                    end

                    documents = get(ws_info, "documents", [])
                    if !isempty(documents)
                        @info "  Unlinking $(length(documents)) documents from workspace before deletion"
                        # Unlink all documents from workspace
                        doc_paths = [doc["docpath"] for doc in documents]
                        remove_body = JSON.json(Dict(
                            "adds" => [],
                            "deletes" => doc_paths
                        ))
                        HTTP.post("$api_url/api/v1/workspace/$(ws["slug"])/update-embeddings",
                                 headers=headers,
                                 body=remove_body)
                        sleep(1)  # Give server time to process
                    end
                catch doc_e
                    @warn "  Error removing documents from workspace" exception=doc_e
                end

                # Now delete the workspace
                HTTP.delete("$api_url/api/v1/workspace/$(ws["slug"])", headers=headers)
                sleep(1)  # Give server time to clean up
                break
            end
        end
    catch e
        @warn "Error checking/deleting workspace" exception=e
    end

    # Delete any existing embeds for this workspace
    try
        response = HTTP.get("$api_url/api/v1/embed", headers=headers)
        embeds = JSON.parse(String(response.body))

        for embed in get(embeds, "embeds", [])
            if haskey(embed, "workspace")
                workspace_name = lowercase(get(embed["workspace"], "name", ""))
                if contains(workspace_name, lowercase(package_name))
                    @info "Deleting existing embed: $(embed["uuid"])"
                    HTTP.delete("$api_url/api/v1/embed/$(embed["uuid"])", headers=headers)
                    break
                end
            end
        end
    catch e
        @warn "Error checking/deleting embed" exception=e
    end

    # Step 2: Create new workspace
    @info "Creating new workspace: $workspace_slug"
    workspace_id = nothing
    try
        body = JSON.json(Dict("name" => package_name))
        response = HTTP.post("$api_url/api/v1/workspace/new",
                           headers=headers,
                           body=body)
        result = JSON.parse(String(response.body))
        workspace_id = result["workspace"]["id"]
        @info "Created workspace: $(result["workspace"]["slug"]) (ID: $workspace_id)"
    catch e
        @error "Failed to create workspace" exception=e
        return ""
    end

    # Step 3: Upload markdown documentation files
    @info "Uploading markdown documentation files..."
    docs_src_dir = joinpath(@__DIR__, "src")
    markdown_files = filter(f -> endswith(f, ".md"), readdir(docs_src_dir))

    for md_file in markdown_files
        file_path = joinpath(docs_src_dir, md_file)
        try
            @info "  Uploading: $md_file"
            upload_document_from_file(api_url, api_key, workspace_slug, file_path, md_file)
        catch e
            @warn "Failed to upload $md_file" exception=e
        end
    end

    # Step 4: Upload public docstrings
    if !isempty(modules)
        @info "Collecting and uploading public docstrings..."
        try
            upload_docstrings(api_url, api_key, workspace_slug, modules)
        catch e
            @warn "Failed to upload docstrings" exception=e
        end
    else
        @info "No modules specified, skipping docstring upload"
    end

    # Step 5: Update embeddings to include all uploaded documents
    @info "Updating embeddings..."
    sleep(2)  # Give server time to process uploads
    try
        update_workspace_embeddings(api_url, api_key, workspace_slug)
    catch e
        @warn "Failed to update embeddings" exception=e
    end

    # Step 6: Get or create embed
    @info "Looking for embed..."
    embed_uuid = ""

    try
        response = HTTP.get("$api_url/api/v1/embed", headers=headers)
        embeds = JSON.parse(String(response.body))

        for embed in get(embeds, "embeds", [])
            if haskey(embed, "workspace")
                workspace_name = lowercase(get(embed["workspace"], "name", ""))
                if contains(workspace_name, lowercase(package_name))
                    embed_uuid = embed["uuid"]
                    @info "Found embed: $embed_uuid"
                    break
                end
            end
        end

        if isempty(embed_uuid)
            @info "No embed found, creating new embed..."
            # Create new embed for the workspace
            embed_body = JSON.json(Dict(
                "workspace_slug" => workspace_slug,
                "max_chats_per_day" => nothing,
                "max_chats_per_session" => nothing
            ))

            create_response = HTTP.post("$api_url/api/v1/embed/new",
                                      headers=headers,
                                      body=embed_body)

            embed_result = JSON.parse(String(create_response.body))

            if haskey(embed_result, "embed") && haskey(embed_result["embed"], "uuid")
                embed_uuid = embed_result["embed"]["uuid"]
                @info "Created new embed: $embed_uuid"
            else
                @error "Failed to create embed: $embed_result"
                error("Embed creation failed")
            end
        end
    catch e
        @warn "Failed to get embed" exception=e
    end

    # Return the embed script
    embed_script = """
    <script
      data-embed-id="$embed_uuid"
      data-base-api-url="$api_url/api/embed"
      src="$api_url/embed/anythingllm-chat-widget.min.js">
    </script>
    """

    @info "AnythingLLM configuration complete!"
    return embed_script
end

"""
    upload_document_from_file(api_url, api_key, workspace_slug, file_path, filename)

Upload a document file to a workspace using the upload-link approach with a data URI.
"""
function upload_document_from_file(api_url, api_key, workspace_slug, file_path, filename)
    headers = Dict(
        "Authorization" => "Bearer $api_key",
        "Content-Type" => "application/json"
    )

    # Read file content
    content = read(file_path, String)

    # Create a temporary file path that can be used as a document location
    # We'll upload as raw text using the workspace upload endpoint
    upload_raw_text(api_url, api_key, workspace_slug, content, filename)
end

"""
    upload_raw_text(api_url, api_key, workspace_slug, content, title)

Upload raw text content as a document to the workspace.
"""
function upload_raw_text(api_url, api_key, workspace_slug, content, title)
    # Create a temporary file with the content
    temp_dir = mktempdir()
    temp_file = joinpath(temp_dir, title)
    write(temp_file, content)

    try
        # Upload using multipart form data
        headers = Dict(
            "Authorization" => "Bearer $api_key"
        )

        # Read file and create multipart form data
        file_io = open(temp_file, "r")
        file_data = read(file_io)
        close(file_io)

        # Create multipart form data with file data as IOBuffer
        form = HTTP.Form(Dict(
            "file" => HTTP.Multipart(title, IOBuffer(file_data), "text/markdown")
        ))

        response = HTTP.post("$api_url/api/v1/document/upload",
                           headers=headers,
                           body=form)

        result = JSON.parse(String(response.body))
        return result
    finally
        # Clean up temporary file
        rm(temp_dir, recursive=true, force=true)
    end
end

"""
    upload_docstrings(api_url, api_key, workspace_slug, modules)

Extract and upload all public docstrings from the specified modules.
"""
function upload_docstrings(api_url, api_key, workspace_slug, modules)
    # Use provided modules
    modules_to_document = modules

    docstrings = Dict{String, String}()

    for mod in modules_to_document
        # Get all names exported by the module
        exported_names = names(mod, all=false, imported=false)

        for name in exported_names
            try
                # Get the binding
                obj = getfield(mod, name)

                # Get documentation
                docs = Docs.doc(obj)

                # Convert to string
                doc_str = string(docs)

                # Skip if no real documentation
                if !isempty(doc_str) && !contains(doc_str, "No documentation found")
                    full_name = "$(mod).$(name)"
                    docstrings[full_name] = doc_str
                end
            catch e
                @debug "Could not get docstring for $name" exception=e
            end
        end
    end

    @info "Found $(length(docstrings)) documented symbols"

    # Upload each docstring as a separate document
    for (name, doc) in docstrings
        filename = "docstring_$(replace(name, "." => "_")).md"
        try
            upload_raw_text(api_url, api_key, workspace_slug, doc, filename)
        catch e
            @debug "Failed to upload docstring for $name" exception=e
        end
    end
end

"""
    update_workspace_embeddings(api_url, api_key, workspace_slug)

Update embeddings for all documents in the workspace.
"""
function update_workspace_embeddings(api_url, api_key, workspace_slug)
    headers = Dict(
        "Authorization" => "Bearer $api_key",
        "Content-Type" => "application/json"
    )

    # Get list of documents
    response = HTTP.get("$api_url/api/v1/workspace/$workspace_slug", headers=headers)
    workspace_data = JSON.parse(String(response.body))

    # Get the documents directory listing
    response = HTTP.get("$api_url/api/v1/documents", headers=headers)
    docs_data = JSON.parse(String(response.body))

    # Collect all document paths
    doc_paths = collect_document_paths(get(docs_data, "localFiles", Dict()))

    # Strip "documents/" prefix if present (AnythingLLM expects paths relative to documents folder)
    doc_paths = [startswith(path, "documents/") ? path[11:end] : path for path in doc_paths]

    if !isempty(doc_paths)
        # Update embeddings to include all documents
        body = JSON.json(Dict(
            "adds" => doc_paths,
            "deletes" => []
        ))

        HTTP.post("$api_url/api/v1/workspace/$workspace_slug/update-embeddings",
                 headers=headers,
                 body=body)

        @info "Updated embeddings with $(length(doc_paths)) documents"
    end
end

"""
    collect_document_paths(file_tree, current_path="")

Recursively collect all document paths from the file tree.
"""
function collect_document_paths(file_tree, current_path="")
    paths = String[]

    if haskey(file_tree, "type") && file_tree["type"] == "folder"
        folder_name = get(file_tree, "name", "")
        new_path = isempty(current_path) ? folder_name : "$current_path/$folder_name"

        for item in get(file_tree, "items", [])
            append!(paths, collect_document_paths(item, new_path))
        end
    elseif haskey(file_tree, "type") && file_tree["type"] == "file"
        file_name = get(file_tree, "name", "")
        if !isempty(file_name)
            full_path = isempty(current_path) ? file_name : "$current_path/$file_name"
            push!(paths, full_path)
        end
    end

    return paths
end
