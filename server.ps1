$listener = New-Object System.Net.HttpListener
$listener.Prefixes.Add('http://localhost:8080/')
$listener.Start()

Write-Host 'Server running at http://localhost:8000/'
Write-Host 'Press Ctrl+C to stop'
Write-Host ''

while ($listener.IsListening) {
    $context = $listener.GetContext()
    $request = $context.Request
    $response = $context.Response
    
    $filepath = Join-Path (Get-Location) ($request.Url.LocalPath.TrimStart('/'))
    
    if (Test-Path $filepath -PathType Container) {
        $filepath = Join-Path $filepath 'index.html'
    }
    
    Write-Host "$($request.HttpMethod) $($request.Url.LocalPath)"
    
    if (Test-Path $filepath -PathType Leaf) {
        $content = [System.IO.File]::ReadAllBytes($filepath)
        $response.ContentLength64 = $content.Length
        
        # Set content type
        $ext = [System.IO.Path]::GetExtension($filepath)
        switch ($ext) {
            '.html' { $response.ContentType = 'text/html' }
            '.css' { $response.ContentType = 'text/css' }
            '.js' { $response.ContentType = 'application/javascript' }
            '.csv' { $response.ContentType = 'text/csv' }
            default { $response.ContentType = 'application/octet-stream' }
        }
        
        $response.OutputStream.Write($content, 0, $content.Length)
        $response.StatusCode = 200
    } else {
        $response.StatusCode = 404
        $buffer = [System.Text.Encoding]::UTF8.GetBytes('404 Not Found')
        $response.OutputStream.Write($buffer, 0, $buffer.Length)
    }
    
    $response.Close()
}

$listener.Stop()
