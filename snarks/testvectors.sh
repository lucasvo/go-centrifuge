#!/usr/bin/env sh
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
BUYER2_RATING="a" # 100 in decimal
BUYER2_HASH=$(echo -n "$BUYER2$BUYER2_PUB$BUYER_RATING" | xxd -r -p | sha256sum | cut -d' ' -f1)

BUYER3="000000000000000000003"
BUYER3_KEY=$($PYCRYPTO keygen)
BUYER3_PUB=$(echo $BUYER3_KEY | cut -d' ' -f2)
BUYER3_PRIV=$(echo $BUYER3_KEY | cut -d' ' -f1)
BUYER3_RATING="a" # 100 in decimal
BUYER3_HASH=$(echo -n "$BUYER3$BUYER3_PUB$BUYER_RATING" | xxd -r -p | sha256sum | cut -d' ' -f1)

BUYER4="000000000000000000004"
BUYER4_KEY=$($PYCRYPTO keygen)
BUYER4_PUB=$(echo $BUYER4_KEY | cut -d' ' -f2)
BUYER4_PRIV=$(echo $BUYER4_KEY | cut -d' ' -f1)
BUYER4_RATING="a" # 100 in decimal
BUYER4_HASH=$(echo -n "$BUYER4$BUYER4_PUB$BUYER_RATING" | xxd -r -p | sha256sum | cut -d' ' -f1)

BUYER_LEAF1=$($PYCRYPTO "hash" $BUYER1_HASH$BUYER2_HASH)
BUYER_LEAF2=$($PYCRYPTO "hash" $BUYER3_HASH$BUYER4_HASH)
BUYER_ROOT=$($PYCRYPTO "hash" $BUYER_LEAF1$BUYER_LEAF2)

cat > "$OUT/proof_data.json" << EOF
{ 
    "public": {
        "credit_rating_roothash": "$BUYER_ROOT",
        "rating": $BUYER1_RATING
    },
    "private": {
        "buyer_rating_proof": { 
            "hashes": ["$BUYER2_HASH", "$BUYER_LEAF2"],
            "right": [false, false],
            "value": "$BUYER1_LEAF"
        }
    }
}
EOF
