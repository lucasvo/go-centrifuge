// +build integration

package receiver_test

import (
	"context"
	"flag"
	"math/big"
	"os"
	"testing"
	"time"

	"github.com/centrifuge/go-centrifuge/crypto/secp256k1"
	"github.com/ethereum/go-ethereum/common/hexutil"

	"github.com/centrifuge/centrifuge-protobufs/gen/go/coredocument"
	"github.com/centrifuge/go-centrifuge/crypto"

	"github.com/centrifuge/go-centrifuge/config/configstore"
	"github.com/ethereum/go-ethereum/common"

	"github.com/centrifuge/centrifuge-protobufs/gen/go/p2p"
	"github.com/centrifuge/go-centrifuge/anchors"
	"github.com/centrifuge/go-centrifuge/bootstrap"
	"github.com/centrifuge/go-centrifuge/bootstrap/bootstrappers/testingbootstrap"
	"github.com/centrifuge/go-centrifuge/config"
	"github.com/centrifuge/go-centrifuge/contextutil"
	cented25519 "github.com/centrifuge/go-centrifuge/crypto/ed25519"
	"github.com/centrifuge/go-centrifuge/documents"
	"github.com/centrifuge/go-centrifuge/documents/purchaseorder"
	"github.com/centrifuge/go-centrifuge/identity"
	"github.com/centrifuge/go-centrifuge/p2p/common"
	"github.com/centrifuge/go-centrifuge/p2p/receiver"
	"github.com/centrifuge/go-centrifuge/protobufs/gen/go/protocol"
	"github.com/centrifuge/go-centrifuge/storage"
	"github.com/centrifuge/go-centrifuge/testingutils/config"
	"github.com/centrifuge/go-centrifuge/testingutils/documents"
	"github.com/centrifuge/go-centrifuge/utils"
	"github.com/golang/protobuf/proto"
	"github.com/stretchr/testify/assert"
)

var (
	handler    *receiver.Handler
	anchorRepo anchors.AnchorRepository
	cfg        config.Configuration
	idService  identity.ServiceDID
	idFactory  identity.Factory
	cfgService config.Service
	docSrv     documents.Service
	defaultDID identity.DID
)

func TestMain(m *testing.M) {
	flag.Parse()
	ctx := testingbootstrap.TestFunctionalEthereumBootstrap()
	cfg = ctx[bootstrap.BootstrappedConfig].(config.Configuration)
	cfgService = ctx[config.BootstrappedConfigStorage].(config.Service)
	docSrv = ctx[documents.BootstrappedDocumentService].(documents.Service)
	anchorRepo = ctx[anchors.BootstrappedAnchorRepo].(anchors.AnchorRepository)
	idService = ctx[identity.BootstrappedDIDService].(identity.ServiceDID)
	idFactory = ctx[identity.BootstrappedDIDFactory].(identity.Factory)
	handler = receiver.New(cfgService, receiver.HandshakeValidator(cfg.GetNetworkID(), idService), docSrv, new(testingdocuments.MockRegistry), idService)
	defaultDID = createIdentity(&testing.T{})
	result := m.Run()
	testingbootstrap.TestFunctionalEthereumTearDown()
	os.Exit(result)
}

func TestHandler_GetDocument_nonexistentIdentifier(t *testing.T) {
	b := utils.RandomSlice(32)
	centrifugeId := createIdentity(t)
	req := &p2ppb.GetDocumentRequest{DocumentIdentifier: b}
	resp, err := handler.GetDocument(context.Background(), req, centrifugeId)
	assert.Error(t, err, "must return error")
	assert.Nil(t, resp, "must be nil")
}

func TestHandler_HandleInterceptorReqSignature(t *testing.T) {
	centID := createIdentity(t)
	tc, err := configstore.NewAccount("", cfg)
	assert.Nil(t, err)
	acc := tc.(*configstore.Account)
	acc.IdentityID = centID[:]
	ctxh, err := contextutil.New(context.Background(), acc)
	assert.Nil(t, err)
	_, err = cfgService.CreateAccount(acc)
	assert.NoError(t, err)
	_, cd := prepareDocumentForP2PHandler(t, nil)
	p2pEnv, err := p2pcommon.PrepareP2PEnvelope(ctxh, cfg.GetNetworkID(), p2pcommon.MessageTypeRequestSignature, &p2ppb.SignatureRequest{Document: &cd})

	pub, _ := acc.GetP2PKeyPair()
	publicKey, err := cented25519.GetPublicSigningKey(pub)
	assert.NoError(t, err)
	var bPk [32]byte
	copy(bPk[:], publicKey)
	peerID, err := cented25519.PublicKeyToP2PKey(bPk)
	assert.NoError(t, err)

	p2pResp, err := handler.HandleInterceptor(ctxh, peerID, p2pcommon.ProtocolForDID(&centID), p2pEnv)
	assert.Nil(t, err, "must be nil")
	assert.NotNil(t, p2pResp, "must be non nil")
	resp := resolveSignatureResponse(t, p2pResp)
	assert.NotNil(t, resp.Signature.Signature, "must be non nil")
	sig := resp.Signature
	assert.True(t, secp256k1.VerifySignatureWithAddress(common.BytesToAddress(sig.PublicKey).String(), hexutil.Encode(sig.Signature), cd.SigningRoot), "signature must be valid")
}

func TestHandler_RequestDocumentSignature_AlreadyExists(t *testing.T) {
	_, cd := prepareDocumentForP2PHandler(t, nil)
	ctxh := testingconfig.CreateAccountContext(t, cfg)
	resp, err := handler.RequestDocumentSignature(ctxh, &p2ppb.SignatureRequest{Document: &cd})
	assert.Nil(t, err, "must be nil")
	assert.NotNil(t, resp, "must be non nil")

	resp, err = handler.RequestDocumentSignature(ctxh, &p2ppb.SignatureRequest{Document: &cd})
	assert.NotNil(t, err, "must not be nil")
	assert.Contains(t, err.Error(), storage.ErrRepositoryModelCreateKeyExists.Error())
}

func TestHandler_RequestDocumentSignature_UpdateSucceeds(t *testing.T) {
	ctxh := testingconfig.CreateAccountContext(t, cfg)
	po, cd := prepareDocumentForP2PHandler(t, nil)
	resp, err := handler.RequestDocumentSignature(ctxh, &p2ppb.SignatureRequest{Document: &cd})
	assert.Nil(t, err, "must be nil")
	assert.NotNil(t, resp, "must be non nil")
	assert.NotNil(t, resp.Signature.Signature, "must be non nil")
	sig := resp.Signature
	assert.True(t, secp256k1.VerifySignatureWithAddress(common.BytesToAddress(sig.PublicKey).String(), hexutil.Encode(sig.Signature), cd.SigningRoot), "signature must be valid")
	//Update document
	po, cd = updateDocumentForP2Phandler(t, po)
	resp, err = handler.RequestDocumentSignature(ctxh, &p2ppb.SignatureRequest{Document: &cd})
	assert.Nil(t, err, "must be nil")
	assert.NotNil(t, resp, "must be non nil")
	assert.NotNil(t, resp.Signature.Signature, "must be non nil")
	sig = resp.Signature
	assert.True(t, secp256k1.VerifySignatureWithAddress(common.BytesToAddress(sig.PublicKey).String(), hexutil.Encode(sig.Signature), cd.SigningRoot), "signature must be valid")
}

func TestHandler_RequestDocumentSignatureFirstTimeOnUpdatedDocument(t *testing.T) {
	ctxh := testingconfig.CreateAccountContext(t, cfg)
	po, cd := prepareDocumentForP2PHandler(t, nil)
	po, cd = updateDocumentForP2Phandler(t, po)
	assert.NotEqual(t, cd.DocumentIdentifier, cd.CurrentVersion)
	resp, err := handler.RequestDocumentSignature(ctxh, &p2ppb.SignatureRequest{Document: &cd})
	assert.Nil(t, err, "must be nil")
	assert.NotNil(t, resp, "must be non nil")
	assert.NotNil(t, resp.Signature.Signature, "must be non nil")
	sig := resp.Signature
	assert.True(t, secp256k1.VerifySignatureWithAddress(common.BytesToAddress(sig.PublicKey).String(), hexutil.Encode(sig.Signature), cd.SigningRoot), "signature must be valid")
}

func TestHandler_RequestDocumentSignature(t *testing.T) {
	ctxh := testingconfig.CreateAccountContext(t, cfg)
	_, cd := prepareDocumentForP2PHandler(t, nil)
	resp, err := handler.RequestDocumentSignature(ctxh, &p2ppb.SignatureRequest{Document: &cd})
	assert.Nil(t, err, "must be nil")
	assert.NotNil(t, resp, "must be non nil")
	assert.NotNil(t, resp.Signature.Signature, "must be non nil")
	sig := resp.Signature
	assert.True(t, secp256k1.VerifySignatureWithAddress(common.BytesToAddress(sig.PublicKey).String(), hexutil.Encode(sig.Signature), cd.SigningRoot), "signature must be valid")
}

func TestHandler_SendAnchoredDocument_update_fail(t *testing.T) {
	_, cd := prepareDocumentForP2PHandler(t, nil)

	// Anchor document
	idConfig, err := identity.GetIdentityConfig(cfg)
	anchorIDTyped, _ := anchors.ToAnchorID(cd.CurrentVersion)
	docRootTyped, _ := anchors.ToDocumentRoot(cd.DocumentRoot)

	ctx := testingconfig.CreateAccountContext(t, cfg)
	anchorConfirmations, err := anchorRepo.CommitAnchor(ctx, anchorIDTyped, docRootTyped, [][anchors.DocumentProofLength]byte{utils.RandomByte32()})
	assert.Nil(t, err)

	watchCommittedAnchor := <-anchorConfirmations
	assert.True(t, watchCommittedAnchor, "No error should be thrown by context")

	anchorResp, err := handler.SendAnchoredDocument(ctx, &p2ppb.AnchorDocumentRequest{Document: &cd}, idConfig.ID[:])
	assert.Error(t, err)
	assert.Contains(t, err.Error(), storage.ErrRepositoryModelUpdateKeyNotFound.Error())
	assert.Nil(t, anchorResp)
}

func TestHandler_SendAnchoredDocument_EmptyDocument(t *testing.T) {
	ctxh := testingconfig.CreateAccountContext(t, cfg)
	id, err := cfg.GetIdentityID()
	assert.NoError(t, err)
	resp, err := handler.SendAnchoredDocument(ctxh, &p2ppb.AnchorDocumentRequest{}, id)
	assert.NotNil(t, err)
	assert.Nil(t, resp, "must be nil")
}

func TestHandler_SendAnchoredDocument(t *testing.T) {
	tc, err := configstore.NewAccount("", cfg)
	assert.Nil(t, err)
	centrifugeId := createIdentity(t)
	acc := tc.(*configstore.Account)
	acc.IdentityID = centrifugeId[:]

	ctxh, err := contextutil.New(context.Background(), tc)
	assert.Nil(t, err)

	po, cd := prepareDocumentForP2PHandler(t, nil)
	resp, err := handler.RequestDocumentSignature(ctxh, &p2ppb.SignatureRequest{Document: &cd})
	assert.Nil(t, err)
	assert.NotNil(t, resp)

	// Add signature received
	po.AppendSignatures(resp.Signature)
	tree, err := po.DocumentRootTree()
	po.Document.DocumentRoot = tree.RootHash()

	// Anchor document
	anchorIDTyped, _ := anchors.ToAnchorID(po.Document.CurrentVersion)
	docRootTyped, _ := anchors.ToDocumentRoot(po.Document.DocumentRoot)
	anchorConfirmations, err := anchorRepo.CommitAnchor(ctxh, anchorIDTyped, docRootTyped, [][anchors.DocumentProofLength]byte{utils.RandomByte32()})
	assert.Nil(t, err)

	watchCommittedAnchor := <-anchorConfirmations
	assert.True(t, watchCommittedAnchor, "No error should be thrown by context")
	cd, err = po.PackCoreDocument()
	assert.NoError(t, err)
	anchorResp, err := handler.SendAnchoredDocument(ctxh, &p2ppb.AnchorDocumentRequest{Document: &cd}, centrifugeId[:])
	assert.Nil(t, err)
	assert.NotNil(t, anchorResp, "must be non nil")
	assert.True(t, anchorResp.Accepted)
}

func createIdentity(t *testing.T) identity.DID {
	// Create Identity
	didAddr, err := idFactory.CalculateIdentityAddress(context.Background())
	assert.NoError(t, err)
	tc, err := configstore.NewAccount("", cfg)
	assert.Nil(t, err)
	acc := tc.(*configstore.Account)
	acc.IdentityID = didAddr.Bytes()

	ctx, err := contextutil.New(context.Background(), tc)
	assert.Nil(t, err)
	did, err := idFactory.CreateIdentity(ctx)
	assert.Nil(t, err, "should not error out when creating identity")
	assert.Equal(t, did.String(), didAddr.String(), "Resulting Identity should have the same ID as the input")

	idConfig, err := identity.GetIdentityConfig(cfg)
	assert.NoError(t, err)
	// Add Keys
	pk, err := utils.SliceToByte32(idConfig.Keys[identity.KeyPurposeP2P].PublicKey)
	assert.NoError(t, err)
	keyDID := identity.NewKey(pk, big.NewInt(identity.KeyPurposeP2P), big.NewInt(identity.KeyTypeECDSA))
	err = idService.AddKey(ctx, keyDID)
	assert.Nil(t, err, "should not error out when adding key to identity")

	sPk, err := utils.SliceToByte32(idConfig.Keys[identity.KeyPurposeSigning].PublicKey)
	assert.NoError(t, err)
	keyDID = identity.NewKey(sPk, big.NewInt(identity.KeyPurposeSigning), big.NewInt(identity.KeyTypeECDSA))
	err = idService.AddKey(ctx, keyDID)
	assert.Nil(t, err, "should not error out when adding key to identity")

	secPk, err := utils.SliceToByte32(idConfig.Keys[identity.KeyPurposeEthMsgAuth].PublicKey)
	assert.NoError(t, err)
	keyDID = identity.NewKey(secPk, big.NewInt(identity.KeyPurposeEthMsgAuth), big.NewInt(identity.KeyTypeECDSA))
	err = idService.AddKey(ctx, keyDID)
	assert.Nil(t, err, "should not error out when adding key to identity")

	return *did
}

func prepareDocumentForP2PHandler(t *testing.T, po *purchaseorder.PurchaseOrder) (*purchaseorder.PurchaseOrder, coredocumentpb.CoreDocument) {
	idConfig, err := identity.GetIdentityConfig(cfg)
	assert.Nil(t, err)
	idConfig.ID = defaultDID
	if po == nil {
		payload := testingdocuments.CreatePOPayload()
		po = new(purchaseorder.PurchaseOrder)
		err = po.InitPurchaseOrderInput(payload, idConfig.ID.String())
		assert.NoError(t, err)
	}
	_, err = po.CalculateDataRoot()
	assert.NoError(t, err)
	sr, err := po.CalculateSigningRoot()
	assert.NoError(t, err)
	s, err := crypto.SignMessage(idConfig.Keys[identity.KeyPurposeSigning].PrivateKey, sr, crypto.CurveSecp256K1)
	assert.NoError(t, err)
	sig := &coredocumentpb.Signature{
		EntityId:  idConfig.ID[:],
		PublicKey: idConfig.Keys[identity.KeyPurposeSigning].PublicKey,
		Signature: s,
		Timestamp: utils.ToTimestamp(time.Now().UTC()),
	}
	po.AppendSignatures(sig)
	_, err = po.CalculateDocumentRoot()
	assert.NoError(t, err)
	cd, err := po.PackCoreDocument()
	assert.NoError(t, err)
	return po, cd
}

func updateDocumentForP2Phandler(t *testing.T, po *purchaseorder.PurchaseOrder) (*purchaseorder.PurchaseOrder, coredocumentpb.CoreDocument) {
	cd, err := po.CoreDocument.PrepareNewVersion(nil, true)
	assert.NoError(t, err)
	po.CoreDocument = cd
	return prepareDocumentForP2PHandler(t, po)
}

func resolveSignatureResponse(t *testing.T, p2pEnv *protocolpb.P2PEnvelope) *p2ppb.SignatureResponse {
	signResp := new(p2ppb.SignatureResponse)
	dataEnv, err := p2pcommon.ResolveDataEnvelope(p2pEnv)
	assert.NoError(t, err)
	err = proto.Unmarshal(dataEnv.Body, signResp)
	assert.NoError(t, err)
	return signResp
}
