"""main entry for microapp command-line interface"""


def main():
    from microapp import MicroappProject
    ret, _ = MicroappProject().run_command()
    return ret

if __name__ == "__main__":
    main()
