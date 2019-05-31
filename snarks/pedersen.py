import argparse
from pycrypto.zokrates.pedersen import (
    pedersen_hash_bytes,
    pedersen_hash_bits,
    pedersen_hash_scalars,
    pedersen_hash_bits_table,
)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("preimage")
    args = parser.parse_args()
    preimage = bytes.fromhex(args.preimage)
    result = pedersen_hash_bytes(b"test", preimage)
    sign = 0
    if (result[0]) > 0: sign = 1
    print(sign, result[1]) 
    
if __name__ == "__main__":
    main()
