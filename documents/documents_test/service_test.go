// +build unit

package documents_test

import (
	"context"
	"os"
	"testing"

	"github.com/centrifuge/centrifuge-protobufs/gen/go/coredocument"
	"github.com/centrifuge/go-centrifuge/anchors"
	"github.com/centrifuge/go-centrifuge/bootstrap"
	"github.com/centrifuge/go-centrifuge/bootstrap/bootstrappers/testlogging"
	"github.com/centrifuge/go-centrifuge/config"
	"github.com/centrifuge/go-centrifuge/contextutil"
	"github.com/centrifuge/go-centrifuge/documents"
	"github.com/centrifuge/go-centrifuge/documents/invoice"
	"github.com/centrifuge/go-centrifuge/errors"
	"github.com/centrifuge/go-centrifuge/ethereum"
	"github.com/centrifuge/go-centrifuge/storage/leveldb"
	"github.com/centrifuge/go-centrifuge/testingutils/commons"
	"github.com/centrifuge/go-centrifuge/testingutils/config"
	"github.com/centrifuge/go-centrifuge/testingutils/documents"
	"github.com/centrifuge/go-centrifuge/testingutils/identity"
	"github.com/centrifuge/go-centrifuge/utils"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/mock"
)

var testRepoGlobal documents.Repository
var (
	cid         = testingidentity.GenerateRandomDID()
	centIDBytes = cid[:]
	tenantID    = cid[:]
	key1Pub     = [...]byte{230, 49, 10, 12, 200, 149, 43, 184, 145, 87, 163, 252, 114, 31, 91, 163, 24, 237, 36, 51, 165, 8, 34, 104, 97, 49, 114, 85, 255, 15, 195, 199}
	key1        = []byte{102, 109, 71, 239, 130, 229, 128, 189, 37, 96, 223, 5, 189, 91, 210, 47, 89, 4, 165, 6, 188, 53, 49, 250, 109, 151, 234, 139, 57, 205, 231, 253, 230, 49, 10, 12, 200, 149, 43, 184, 145, 87, 163, 252, 114, 31, 91, 163, 24, 237, 36, 51, 165, 8, 34, 104, 97, 49, 114, 85, 255, 15, 195, 199}
)

var ctx = map[string]interface{}{}
var cfg config.Configuration

func TestMain(m *testing.M) {
	ethClient := &testingcommons.MockEthClient{}
	ethClient.On("GetEthClient").Return(nil)
	ctx[ethereum.BootstrappedEthereumClient] = ethClient
	ibootstappers := []bootstrap.TestBootstrapper{
		&testlogging.TestLoggingBootstrapper{},
		&config.Bootstrapper{},
	}
	bootstrap.RunTestBootstrappers(ibootstappers, ctx)
	cfg = ctx[bootstrap.BootstrappedConfig].(config.Configuration)
	cfg.Set("identityId", cid.String())
	result := m.Run()
	bootstrap.RunTestTeardown(ibootstappers)
	os.Exit(result)
}

func TestService_ReceiveAnchoredDocumentFailed(t *testing.T) {
	poSrv := documents.DefaultService(nil, nil, documents.NewServiceRegistry(), nil)

	// self failed
	err := poSrv.ReceiveAnchoredDocument(context.Background(), nil, nil)
	assert.Error(t, err)
}

func getServiceWithMockedLayers() (documents.Service, testingcommons.MockIdentityService) {
	repo := testRepo()
	idService := testingcommons.MockIdentityService{}
	idService.On("ValidateSignature", mock.Anything, mock.Anything).Return(nil).Once()
	mockAnchor = &mockAnchorRepo{}
	return documents.DefaultService(repo, mockAnchor, documents.NewServiceRegistry(), &idService), idService
}

type mockAnchorRepo struct {
	mock.Mock
	anchors.AnchorRepository
}

var mockAnchor *mockAnchorRepo

func (r *mockAnchorRepo) GetDocumentRootOf(anchorID anchors.AnchorID) (anchors.DocumentRoot, error) {
	args := r.Called(anchorID)
	docRoot, _ := args.Get(0).(anchors.DocumentRoot)
	return docRoot, args.Error(1)
}

// Functions returns service mocks
func mockSignatureCheck(t *testing.T, i *invoice.Invoice, idService testingcommons.MockIdentityService) testingcommons.MockIdentityService {
	anchorID, _ := anchors.ToAnchorID(i.ID())
	dr, err := i.CalculateDocumentRoot()
	assert.NoError(t, err)
	docRoot, err := anchors.ToDocumentRoot(dr)
	assert.NoError(t, err)
	mockAnchor.On("GetDocumentRootOf", anchorID).Return(docRoot, nil).Once()
	return idService
}

func TestService_CreateProofs(t *testing.T) {
	service, idService := getServiceWithMockedLayers()
	ctxh := testingconfig.CreateAccountContext(t, cfg)
	i, _ := createCDWithEmbeddedInvoice(t, ctxh, false)
	idService = mockSignatureCheck(t, i.(*invoice.Invoice), idService)
	proof, err := service.CreateProofs(ctxh, i.ID(), []string{"invoice.invoice_number"})
	assert.Nil(t, err)
	assert.Equal(t, i.ID(), proof.DocumentID)
	assert.Equal(t, i.CurrentVersion(), proof.VersionID)
	assert.Equal(t, len(proof.FieldProofs), 1)
	assert.Equal(t, proof.FieldProofs[0].GetCompactName(), []byte{0x0, 0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x1})
}
func TestService_CreateProofsValidationFails(t *testing.T) {
	service, idService := getServiceWithMockedLayers()
	ctxh := testingconfig.CreateAccountContext(t, cfg)
	i, _ := createCDWithEmbeddedInvoice(t, ctxh, false)
	idService = mockSignatureCheck(t, i.(*invoice.Invoice), idService)
	i.(*invoice.Invoice).Document.DataRoot = nil
	i.(*invoice.Invoice).Document.SigningRoot = nil
	assert.Nil(t, testRepo().Update(tenantID, i.CurrentVersion(), i))
	_, err := service.CreateProofs(ctxh, i.ID(), []string{"invoice.invoice_number"})
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "failed to get signing root")
}

func TestService_CreateProofsInvalidField(t *testing.T) {
	service, idService := getServiceWithMockedLayers()
	ctxh := testingconfig.CreateAccountContext(t, cfg)
	i, _ := createCDWithEmbeddedInvoice(t, ctxh, false)
	idService = mockSignatureCheck(t, i.(*invoice.Invoice), idService)
	_, err := service.CreateProofs(ctxh, i.CurrentVersion(), []string{"invalid_field"})
	assert.Error(t, err)
	assert.True(t, errors.IsOfType(documents.ErrDocumentProof, err))
}

func TestService_CreateProofsDocumentDoesntExist(t *testing.T) {
	service, _ := getServiceWithMockedLayers()
	ctxh := testingconfig.CreateAccountContext(t, cfg)
	_, err := service.CreateProofs(ctxh, utils.RandomSlice(32), []string{"invoice.invoice_number"})
	assert.Error(t, err)
	assert.True(t, errors.IsOfType(documents.ErrDocumentNotFound, err))
}

func TestService_CreateProofsForVersion(t *testing.T) {
	service, idService := getServiceWithMockedLayers()
	ctxh := testingconfig.CreateAccountContext(t, cfg)
	i, _ := createCDWithEmbeddedInvoice(t, ctxh, false)
	idService = mockSignatureCheck(t, i.(*invoice.Invoice), idService)
	proof, err := service.CreateProofsForVersion(ctxh, i.ID(), i.CurrentVersion(), []string{"invoice.invoice_number"})
	assert.Nil(t, err)
	assert.Equal(t, i.ID(), proof.DocumentID)
	assert.Equal(t, i.CurrentVersion(), proof.VersionID)
	assert.Equal(t, len(proof.FieldProofs), 1)
	assert.Equal(t, proof.FieldProofs[0].GetCompactName(), []byte{0x0, 0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x1})
}

func TestService_RequestDocumentSignature_SigningRootNil(t *testing.T) {
	service, idService := getServiceWithMockedLayers()
	ctxh := testingconfig.CreateAccountContext(t, cfg)
	i, _ := createCDWithEmbeddedInvoice(t, ctxh, true)
	idService = mockSignatureCheck(t, i.(*invoice.Invoice), idService)
	i.(*invoice.Invoice).Document.DataRoot = nil
	i.(*invoice.Invoice).Document.SigningRoot = nil
	signature, err := service.RequestDocumentSignature(ctxh, i)
	assert.NotNil(t, err)
	assert.True(t, errors.IsOfType(documents.ErrDocumentInvalid, err))
	assert.Nil(t, signature)
}

func TestService_CreateProofsForVersionDocumentDoesntExist(t *testing.T) {
	ctxh := testingconfig.CreateAccountContext(t, cfg)
	i, _ := createCDWithEmbeddedInvoice(t, ctxh, false)
	s, _ := getServiceWithMockedLayers()
	_, err := s.CreateProofsForVersion(ctxh, i.ID(), utils.RandomSlice(32), []string{"invoice.invoice_number"})
	assert.Error(t, err)
	assert.True(t, errors.IsOfType(documents.ErrDocumentVersionNotFound, err))
}

func TestService_GetCurrentVersion_successful(t *testing.T) {
	service, _ := getServiceWithMockedLayers()
	documentIdentifier := utils.RandomSlice(32)
	const amountVersions = 10

	version := documentIdentifier
	var currentVersion []byte

	nonExistingVersion := utils.RandomSlice(32)

	for i := 0; i < amountVersions; i++ {

		var next []byte
		if i != amountVersions-1 {
			next = utils.RandomSlice(32)
		} else {
			next = nonExistingVersion
		}

		cd := coredocumentpb.CoreDocument{
			DocumentIdentifier: documentIdentifier,
			CurrentVersion:     version,
			NextVersion:        next,
		}

		inv := &invoice.Invoice{
			GrossAmount:  int64(i + 1),
			CoreDocument: documents.NewCoreDocumentFromProtobuf(cd),
		}

		err := testRepo().Create(tenantID, version, inv)
		currentVersion = version
		version = next
		assert.Nil(t, err)

	}

	ctxh := testingconfig.CreateAccountContext(t, cfg)
	model, err := service.GetCurrentVersion(ctxh, documentIdentifier)
	assert.Nil(t, err)
	assert.Equal(t, model.CurrentVersion(), currentVersion, "should return latest version")
	assert.Equal(t, model.NextVersion(), nonExistingVersion, "latest version should have a non existing id as nextVersion ")

}

func TestService_GetVersion_successful(t *testing.T) {
	service, _ := getServiceWithMockedLayers()
	documentIdentifier := utils.RandomSlice(32)
	currentVersion := utils.RandomSlice(32)
	cd := coredocumentpb.CoreDocument{
		DocumentIdentifier: documentIdentifier,
		CurrentVersion:     currentVersion,
	}
	inv := &invoice.Invoice{
		GrossAmount:  60,
		CoreDocument: documents.NewCoreDocumentFromProtobuf(cd),
	}

	ctxh := testingconfig.CreateAccountContext(t, cfg)
	err := testRepo().Create(tenantID, currentVersion, inv)
	assert.Nil(t, err)

	mod, err := service.GetVersion(ctxh, documentIdentifier, currentVersion)
	assert.Nil(t, err)

	assert.Equal(t, documentIdentifier, mod.ID(), "should be same document Identifier")
	assert.Equal(t, currentVersion, mod.CurrentVersion(), "should be same version")
}

func TestService_GetCurrentVersion_error(t *testing.T) {
	service, _ := getServiceWithMockedLayers()
	documentIdentifier := utils.RandomSlice(32)

	//document is not existing
	ctxh := testingconfig.CreateAccountContext(t, cfg)
	_, err := service.GetCurrentVersion(ctxh, documentIdentifier)
	assert.True(t, errors.IsOfType(documents.ErrDocumentVersionNotFound, err))

	cd := coredocumentpb.CoreDocument{
		DocumentIdentifier: documentIdentifier,
		CurrentVersion:     documentIdentifier,
	}

	inv := &invoice.Invoice{
		GrossAmount:  60,
		CoreDocument: documents.NewCoreDocumentFromProtobuf(cd),
	}

	err = testRepo().Create(tenantID, documentIdentifier, inv)
	assert.Nil(t, err)

	_, err = service.GetCurrentVersion(ctxh, documentIdentifier)
	assert.Nil(t, err)

}

func TestService_GetVersion_error(t *testing.T) {
	service, _ := getServiceWithMockedLayers()

	documentIdentifier := utils.RandomSlice(32)
	currentVersion := utils.RandomSlice(32)

	//document is not existing
	ctxh := testingconfig.CreateAccountContext(t, cfg)
	_, err := service.GetVersion(ctxh, documentIdentifier, currentVersion)
	assert.True(t, errors.IsOfType(documents.ErrDocumentVersionNotFound, err))

	cd := coredocumentpb.CoreDocument{
		DocumentIdentifier: documentIdentifier,
		CurrentVersion:     currentVersion,
	}
	inv := &invoice.Invoice{
		GrossAmount:  60,
		CoreDocument: documents.NewCoreDocumentFromProtobuf(cd),
	}
	err = testRepo().Create(tenantID, currentVersion, inv)
	assert.Nil(t, err)

	//random version
	_, err = service.GetVersion(ctxh, documentIdentifier, utils.RandomSlice(32))
	assert.True(t, errors.IsOfType(documents.ErrDocumentVersionNotFound, err))

	//random document id
	_, err = service.GetVersion(ctxh, utils.RandomSlice(32), documentIdentifier)
	assert.True(t, errors.IsOfType(documents.ErrDocumentVersionNotFound, err))
}

func testRepo() documents.Repository {
	if testRepoGlobal == nil {
		ldb, err := leveldb.NewLevelDBStorage(leveldb.GetRandomTestStoragePath())
		if err != nil {
			panic(err)
		}
		testRepoGlobal = documents.NewDBRepository(leveldb.NewLevelDBRepository(ldb))
		testRepoGlobal.Register(&invoice.Invoice{})
	}
	return testRepoGlobal
}

func TestService_Exists(t *testing.T) {
	service, _ := getServiceWithMockedLayers()
	documentIdentifier := utils.RandomSlice(32)
	ctxh := testingconfig.CreateAccountContext(t, cfg)

	//document is not existing
	_, err := service.GetCurrentVersion(ctxh, documentIdentifier)
	assert.True(t, errors.IsOfType(documents.ErrDocumentVersionNotFound, err))

	cd := coredocumentpb.CoreDocument{
		DocumentIdentifier: documentIdentifier,
		CurrentVersion:     documentIdentifier,
	}
	inv := &invoice.Invoice{
		GrossAmount:  60,
		CoreDocument: documents.NewCoreDocumentFromProtobuf(cd),
	}

	err = testRepo().Create(tenantID, documentIdentifier, inv)

	exists := service.Exists(ctxh, documentIdentifier)
	assert.True(t, exists, "document should exist")

	exists = service.Exists(ctxh, utils.RandomSlice(32))
	assert.False(t, exists, "document should not exist")

}

func createCDWithEmbeddedInvoice(t *testing.T, ctx context.Context, skipSave bool) (documents.Model, coredocumentpb.CoreDocument) {
	i := new(invoice.Invoice)
	err := i.InitInvoiceInput(testingdocuments.CreateInvoicePayload(), cid.String())
	assert.NoError(t, err)
	_, err = i.CalculateDataRoot()
	assert.NoError(t, err)
	sr, err := i.CalculateSigningRoot()
	assert.NoError(t, err)

	acc, err := contextutil.Account(ctx)
	assert.NoError(t, err)

	sig, err := acc.SignMsg(sr)
	assert.NoError(t, err)

	i.AppendSignatures(sig)
	_, err = i.CalculateDocumentRoot()
	assert.NoError(t, err)
	cd, err := i.PackCoreDocument()
	assert.NoError(t, err)

	if !skipSave {
		err = testRepo().Create(tenantID, i.CurrentVersion(), i)
		assert.NoError(t, err)
	}
	return i, cd
}
