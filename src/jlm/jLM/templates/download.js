var linkData = window.atob('{{ blob }}');

var filename = '{{ filename }}';
var blob = new Blob([linkData], { type: 'application/octet-stream' });
if (navigator.msSaveBlob) { // IE 10+
    navigator.msSaveBlob(blob, filename);
} else {
    var link = document.createElement("a");
    if (link.download !== undefined) { // feature detection
        // Browsers that support HTML5 download attribute
        var url = URL.createObjectURL(blob);
        link.setAttribute("href", url);
        link.setAttribute("download", filename);
        link.style.visibility = 'hidden';
        document.body.appendChild(link);
        link.click();
        link.setAttribute("href", "#");
        document.body.removeChild(link);
    }
}
