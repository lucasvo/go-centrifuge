// +build integration

package anchors_test

import (
	"context"
	"os"
	"testing"

	"github.com/centrifuge/go-centrifuge/anchors"
	"github.com/centrifuge/go-centrifuge/bootstrap"
	cc "github.com/centrifuge/go-centrifuge/bootstrap/bootstrappers/testingbootstrap"
	"github.com/centrifuge/go-centrifuge/config"
	"github.com/centrifuge/go-centrifuge/ethereum"
	"github.com/centrifuge/go-centrifuge/testingutils/config"
	"github.com/centrifuge/go-centrifuge/utils"
	"github.com/stretchr/testify/assert"
)

var (
	anchorRepo anchors.AnchorRepository
	cfg        config.Configuration
)

func TestMain(m *testing.M) {
	ctx := cc.TestFunctionalEthereumBootstrap()
	anchorRepo = ctx[anchors.BootstrappedAnchorRepo].(anchors.AnchorRepository)
	cfg = ctx[bootstrap.BootstrappedConfig].(config.Configuration)
	result := m.Run()
	cc.TestFunctionalEthereumTearDown()
	os.Exit(result)
}

func TestCommitAnchor_Integration(t *testing.T) {
	anchorID := utils.RandomSlice(32)
	documentRoot := utils.RandomSlice(32)

	anchorIDTyped, _ := anchors.ToAnchorID(anchorID)
	docRootTyped, _ := anchors.ToDocumentRoot(documentRoot)
	commitAnchor(t, anchorID, documentRoot, [][anchors.DocumentProofLength]byte{utils.RandomByte32()})
	gotDocRoot, err := anchorRepo.GetDocumentRootOf(anchorIDTyped)
	assert.Nil(t, err)
	assert.Equal(t, docRootTyped, gotDocRoot)
}

func commitAnchor(t *testing.T, anchorID, documentRoot []byte, documentProofs [][32]byte) {
	anchorIDTyped, _ := anchors.ToAnchorID(anchorID)
	docRootTyped, _ := anchors.ToDocumentRoot(documentRoot)

	ctx := testingconfig.CreateAccountContext(t, cfg)
	done, err := anchorRepo.CommitAnchor(ctx, anchorIDTyped, docRootTyped, documentProofs)

	isDone := <-done

	assert.True(t, isDone, "isDone should be true")

	assert.Nil(t, err)

}

func TestCommitAnchor_Integration_Concurrent(t *testing.T) {
	var commitDataList [5]*anchors.CommitData
	var doneList [5]chan bool

	for ix := 0; ix < 5; ix++ {
		currentAnchorId := utils.RandomByte32()
		currentDocumentRoot := utils.RandomByte32()
		documentProofs := [][anchors.DocumentProofLength]byte{utils.RandomByte32()}
		h, err := ethereum.GetClient().GetEthClient().HeaderByNumber(context.Background(), nil)
		assert.Nil(t, err, " error must be nil")
		commitDataList[ix] = anchors.NewCommitData(h.Number.Uint64(), currentAnchorId, currentDocumentRoot, documentProofs)
		ctx := testingconfig.CreateAccountContext(t, cfg)
		doneList[ix], err = anchorRepo.CommitAnchor(ctx, currentAnchorId, currentDocumentRoot, documentProofs)
		if err != nil {
			t.Fatalf("Error commit Anchor %v", err)
		}
		assert.Nil(t, err, "should not error out upon anchor registration")
	}

	for ix := 0; ix < 5; ix++ {
		isDone := <-doneList[ix]
		assert.True(t, isDone)
		anchorID := commitDataList[ix].AnchorID
		docRoot := commitDataList[ix].DocumentRoot
		gotDocRoot, err := anchorRepo.GetDocumentRootOf(anchorID)
		assert.Nil(t, err)
		assert.Equal(t, docRoot, gotDocRoot)
	}
}
