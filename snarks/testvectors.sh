#!/usr/bin/env sh
# Requires seth to be installed, see: https://dapp.tools/
set -ex

OUT=out/

PYCRYPTO="python3 pycrypto/cli.py" 

BUYER1="0000000000000000001"
BUYER1_KEY=$($PYCRYPTO keygen)
BUYER1_PUB=$(echo $BUYER1_KEY | cut -d' ' -f2)
BUYER1_PRIV=$(echo $BUYER1_KEY | cut -d' ' -f1)
BUYER1_RATING="64" # 100 in decimal
BUYER1_LEAF="$BUYER1$BUYER1_PUB$BUYER_RATING"
BUYER1_HASH=$(echo -n $BUYER1_LEAF | xxd -r -p | sha256sum | cut -d' ' -f1)

BUYER2="000000000000000000002"
BUYER2_KEY=$($PYCRYPTO keygen)
BUYER2_PUB=$(echo $BUYER2_KEY | cut -d' ' -f2)
BUYER2_PRIV=$(echo $BUYER2_KEY | cut -d' ' -f1)
BUYER2_RATING="0a" # 10 in decimal
BUYER2_HASH=$(echo -n "$BUYER2$BUYER2_PUB$BUYER_RATING" | xxd -r -p | sha256sum | cut -d' ' -f1)

BUYER3="000000000000000000003"
BUYER3_KEY=$($PYCRYPTO keygen)
BUYER3_PUB=$(echo $BUYER3_KEY | cut -d' ' -f2)
BUYER3_PRIV=$(echo $BUYER3_KEY | cut -d' ' -f1)
BUYER3_RATING="0a" # 10 in decimal
BUYER3_HASH=$(echo -n "$BUYER3$BUYER3_PUB$BUYER_RATING" | xxd -r -p | sha256sum | cut -d' ' -f1)

BUYER4="000000000000000000004"
BUYER4_KEY=$($PYCRYPTO keygen)
BUYER4_PUB=$(echo $BUYER4_KEY | cut -d' ' -f2)
BUYER4_PRIV=$(echo $BUYER4_KEY | cut -d' ' -f1)
BUYER4_RATING="60" # 96 in decimal
BUYER4_HASH=$(echo -n "$BUYER4$BUYER4_PUB$BUYER_RATING" | xxd -r -p | sha256sum | cut -d' ' -f1)

BUYER_NODE1=$($PYCRYPTO hash $BUYER1_HASH$BUYER2_HASH)
BUYER_NODE2=$($PYCRYPTO hash $BUYER3_HASH$BUYER4_HASH)
BUYER_ROOT=$($PYCRYPTO hash $BUYER_NODE1$BUYER_NODE2)

NFT_AMOUNT=$(printf "%064x" 800) # 800 in hex


DOC_INVOICE_AMOUNT=$(printf "%064x" 1000) # 1000 in hex
DOC_INVOICE_AMOUNT_SALT=$(dd if=/dev/urandom bs=32 count=1 | xxd -ps -c 200 | tr -d '\n')
DOC_INVOICE_AMOUNT_PROPERTY=00000000000000000
DOC_INVOICE_AMOUNT_LEAF="$DOC_INVOICE_AMOUNT_PROPERTY$DOC_INVOICE_AMOUNT_SALT"
DOC_INVOICE_AMOUNT_HASH=$(echo -n "$DOC_INVOICE_AMOUNT_LEAF" |  xxd -r -p | sha256sum | cut -d' ' -f1)

DOC_BUYER_SALT=$(dd if=/dev/urandom bs=32 count=1 | xxd -ps -c 200 | tr -d '\n')
DOC_BUYER_PROPERTY=0000000000
DOC_BUYER_LEAF="$DOC_BUYER_PROPERTY$BUYER1$DOC_BUYER_SALT"
DOC_BUYER_HASH=$(echo -n "$DOC_BUYER_LEAF" |  xxd -r -p | sha256sum | cut -d' ' -f1)

DOC_EXTRA_LEAF2=$(dd if=/dev/urandom bs=32 count=1 | xxd -ps -c 200 | tr -d '\n')
DOC_EXTRA_LEAF3=$(dd if=/dev/urandom bs=32 count=1 | xxd -ps -c 200 | tr -d '\n')
DOC_EXTRA_LEAF4=$(dd if=/dev/urandom bs=32 count=1 | xxd -ps -c 200 | tr -d '\n')
DOC_EXTRA_LEAF5=$(dd if=/dev/urandom bs=32 count=1 | xxd -ps -c 200 | tr -d '\n')
DOC_EXTRA_LEAF6=$(dd if=/dev/urandom bs=32 count=1 | xxd -ps -c 200 | tr -d '\n')
DOC_EXTRA_LEAF7=$(dd if=/dev/urandom bs=32 count=1 | xxd -ps -c 200 | tr -d '\n')

DOC_NODE1=$($PYCRYPTO hash $DOC_INVOICE_AMOUNT_HASH$DOC_BUYER_HASH)
DOC_NODE2=$($PYCRYPTO hash $DOC_EXTRA_LEAF2$DOC_EXTRA_LEAF3)
DOC_NODE3=$($PYCRYPTO hash $DOC_EXTRA_LEAF4$DOC_EXTRA_LEAF5)
DOC_NODE4=$($PYCRYPTO hash $DOC_EXTRA_LEAF6$DOC_EXTRA_LEAF7)

DOC_NODE4=$($PYCRYPTO hash $DOC_NODE1$DOC_NODE2)
DOC_NODE5=$($PYCRYPTO hash $DOC_NODE3$DOC_NODE4)

DOC_ROOT=$($PYCRYPTO hash $DOC_NODE4$DOC_NODE5)

BUYER1_SIG=$($PYCRYPTO sig-gen $BUYER1_PRIV $DOC_ROOT)
BUYER1_SIG=$(echo $BUYER1_SIG | cut -d' ' -f1)$(echo $BUYER1_SIG | cut -d' ' -f2)

cat > "$OUT/proof_data.json" << EOF
{ 
    "public": {
        "nft_amount": "$NFT_AMOUNT",
        "credit_rating_roothash": "$BUYER_ROOT",
        "rating": "$BUYER1_RATING",
	"document_roothash": "$DOC_ROOT"
    },
    "private": {
	"buyer_signature": "$BUYER1_SIG",
        "buyer_rating_proof": { 
            "hashes": ["$BUYER2_HASH", "$BUYER_LEAF2"],
            "right": [false, false],
            "value": "$BUYER1_LEAF"
        },
        "document_invoice_amount_proof": {
            "hashes": ["$DOC_BUYER_HASH", "$DOC_NODE2", "$DOC_NODE5"],
            "right": [false, false, false],
            "value": "$DOC_INVOICE_AMOUNT",
            "salt": "$DOC_INVOICE_AMOUNT_SALT",
            "property": "$DOC_INVOICE_AMOUNT_PROPERTY"
        "document_invoice_buyer_proof": {
            "hashes": ["$ODC_INVOICE_AMOUNT_HASH", "$DOC_NODE2", "$DOC_NODE5"],
            "right": [true, false, false],
            "value": "$BUYER1",
            "salt": "$DOC_BUYER_SALT",
            "property": "$DOC_BUYER_PROPERTY"
        }
    }
}
EOF
