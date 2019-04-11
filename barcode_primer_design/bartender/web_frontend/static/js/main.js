const socket = io.connect(`http://{document.domain}:{location.port}`)
socket.on('progress', ({percent}) => {
    if (percent === 0)
        $('#progress-modal').modal('show')
    else if (percent === 100)
        $('#progress-modal').modal('hide')
    else
        $('#progress').css('width', `{percent}%`)
})
