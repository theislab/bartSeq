import subprocess

from .p3_parser import parse_p3_information, parse_p3seq_information


def run_and_feed(
    cmd: str,
    *,
    _input_str: str,
    _encoding: str = "utf-8",
    _long_arg_prefix: str = "--",
    **kwargs,
) -> subprocess.CompletedProcess:
    arg_prefix = lambda a: _long_arg_prefix if len(a) > 1 else "-"
    args = [
        s
        for k, v in kwargs.items()
        for s in [f"{arg_prefix(k)}{k}", *([] if v is True else [str(v)])]
    ]
    try:
        return subprocess.run(
            [cmd, *args],
            input=_input_str,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding=_encoding,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"Could not run {cmd} with arguments {args} "
            f"and the following input:\n\n{_input_str}\n\n"
            f"Returncode: {e.returncode}\n\n"
            f"STDOUT:\n\n{e.stdout}\n\n"
            f"STDERR:\n\n{e.stderr}\n\n"
        ) from e
