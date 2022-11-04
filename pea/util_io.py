import asyncio
import sys
from asyncio.subprocess import PIPE, STDOUT
from pathlib import Path

def parse_fastqjoin_group_dict(file_group_rows):
    file_group_dict={}
    for _, row in file_group_rows.iterrows():
        beg, end = row.range.split('~')
        files = [Path(row.data_dir)/f'{i}.fastqjoin' for i in range(int(beg), int(end)+1, 1)]
        file_group_dict[row.key] = files
        
    return file_group_dict


def gen_input_fastqjoins(input_cells, group_dict):
    for val in input_cells:
        if val in group_dict:
            for input_file in group_dict[val]:
                yield input_file
        else:
            input_file = Path(val)
            yield input_file

            
@asyncio.coroutine
def run_async_command(cmd):
    tmp=cmd.split()
    exe, args = tmp[0], tmp[1:]
    p = yield from asyncio.create_subprocess_exec(exe, *args,
            stdin=PIPE, stdout=PIPE, stderr=STDOUT, cwd="./")
    return (yield from p.communicate())[0].splitlines()


def get_async_eventloop():
    if sys.platform.startswith('win'):
        loop = asyncio.ProactorEventLoop() # for subprocess' pipes on Windows
        asyncio.set_event_loop(loop)
    else:
        loop = asyncio.get_event_loop()
    
    return loop


def run_async_commands(loop, cmds):
    coros=[run_async_command(cmd) for cmd in cmds]
        
    outputs = loop.run_until_complete(asyncio.gather(*coros))

