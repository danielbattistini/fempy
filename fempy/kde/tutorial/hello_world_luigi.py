# run with python3 -m luigi --module hello_world_luigi HelloWorldTask  --local-scheduler

import luigi
import sys

class HelloTask(luigi.Task):
    def run(self):
        with open('hello.txt', 'w') as hello_file:
            hello_file.write('hello')
            hello_file.close()

    def output(self):
        return luigi.LocalTarget('hello.txt')


class WorldTask(luigi.Task):
    def run(self):
        with open('world.txt', 'w') as world_file:
            world_file.write('world')
            world_file.close()

    def output(self):
        return luigi.LocalTarget('world.txt')

class HelloWorldTask(luigi.Task):
    def run(self):
        with open('hello.txt', 'r') as hello_file:
            hello = hello_file.read()
        with open('world.txt', 'r') as world_file:
            world = world_file.read()
        with open('hello_world.txt', 'w') as output_file:
            content = f'{hello} {world}'
            output_file.write(content)
            output_file.close()

    def requires(self):
        return[HelloTask(), WorldTask()]

    def output(self):
        return luigi.LocalTarget('hello_world.txt')



if __name__ == '__main__':
    luigi.run()