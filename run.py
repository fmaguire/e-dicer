import eDicer.eDicer as ed

if __name__=='__main__':

    parser = ed.get_parser()
    args = parser.parse_args()
    ed.main(args.input_file, args.output_file, k=args.k)
