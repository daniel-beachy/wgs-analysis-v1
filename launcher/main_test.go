package main

import (
	"fmt"
	"net"
	"net/http"
	"testing"
)

// A running dashboard on the preferred port is reused only when it serves the same data folder: the public demo
// must never open on top of a personal genome (or the other way round).
func TestListenReusesOnlySameData(t *testing.T) {
	busy, err := net.Listen("tcp", "127.0.0.1:0")
	if err != nil {
		t.Fatal(err)
	}
	defer busy.Close()
	mux := http.NewServeMux()
	mux.HandleFunc("/api/ping", func(w http.ResponseWriter, _ *http.Request) {
		w.Header().Set(dataHeader, dataID("/drive/wgs-data/releases"))
		fmt.Fprint(w, "genome-dashboard test")
	})
	go http.Serve(busy, mux)
	port := busy.Addr().(*net.TCPAddr).Port

	if _, _, existing := listen(port, dataID("/drive/wgs-data/releases")); existing == "" {
		t.Fatal("same data folder: expected the running instance to be reused")
	}
	ln, url, existing := listen(port, dataID("/drive/wgs-data-demo/releases"))
	if existing != "" || ln == nil {
		t.Fatalf("other data folder: expected a new port, got existing=%q", existing)
	}
	defer ln.Close()
	if url == fmt.Sprintf("http://127.0.0.1:%d/", port) {
		t.Fatal("other data folder: must not share the busy port")
	}
}
