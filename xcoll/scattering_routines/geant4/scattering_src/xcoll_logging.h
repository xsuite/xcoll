// root_logging.h
#pragma once

#include "TSystem.h"
#include "TError.h"
#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include <fstream>
#include <memory>
#include <mutex>
#include <unistd.h>
#include <fcntl.h>


#include <filesystem>
namespace fs = std::filesystem;

static fs::path resolve_path(const fs::path& workdir, const fs::path& name){
    if (workdir.empty() || name.is_absolute()) return name;
    return workdir / name;
}


// ------------------------------ //
// Logging functionality for Root //
// ------------------------------ //

namespace rootlog {
  inline std::unique_ptr<std::ofstream> info;
  inline std::unique_ptr<std::ofstream> err;
  inline ErrorHandlerFunc_t prev = nullptr;
  inline std::mutex mtx;

  inline void Handler(int level, Bool_t abort, const char* loc, const char* msg) {
    std::lock_guard<std::mutex> lock(mtx);
    if (!info || !err) return;
    std::ostream& os = (level >= kError) ? *err : *info;
    os << loc << ": " << msg << '\n';
    os.flush();
    // if (abort) ::abort();
  }

  inline void init(const std::string& infoPath,
                   const std::string& errPath,
                   bool append = true) {
    std::lock_guard<std::mutex> lock(mtx);

    auto mode = append ? std::ios::app : std::ios::trunc;
    info = std::make_unique<std::ofstream>(infoPath, mode);
    err  = std::make_unique<std::ofstream>(errPath.empty() ? infoPath : errPath, mode);
    if (!info->is_open() || !err->is_open())
      throw std::runtime_error("Failed to open ROOT log files");

    if (!prev) prev = GetErrorHandler();    // save once
    SetErrorHandler(&Handler);              // (re)install handler
  }

  inline void restore() {
    std::lock_guard<std::mutex> lock(mtx);
    if (prev) SetErrorHandler(prev);
    info.reset();
    err.reset();
  }
}



// -------------------------------- //
// Logging functionality for Geant4 //
// -------------------------------- //

class XcollLogFileSession : public G4UIsession {
public:
  XcollLogFileSession(const std::string& out_path,
                      const std::string& err_path)
  : logfile(out_path, std::ios::app),
    errfile(err_path, std::ios::app)
  {}

  ~XcollLogFileSession() override { logfile.close(); errfile.close(); }

  G4int ReceiveG4cout(const G4String& msg) override {
    logfile << msg;
    logfile.flush();
    return 0;
  }
  G4int ReceiveG4cerr(const G4String& msg) override {
    errfile << msg;
    errfile.flush();
    return 0;
  }

private:
  std::ofstream logfile;
  std::ofstream errfile;
};

inline void RedirectGeant4(const std::string& out_path,
                           const std::string& err_path) {
    auto ui = G4UImanager::GetUIpointer();
    if (!ui) {return;}
    static std::unique_ptr<XcollLogFileSession> session;
    session = std::make_unique<XcollLogFileSession>(out_path, err_path);
    ui->SetCoutDestination(session.get());
}



// ---------------------------------------------- //
// Logging functionality for others (e.g. Geant4) //
// ---------------------------------------------- //

class FDRedirect {
    int saved_out{-1}, saved_err{-1}, out_fd{-1}, err_fd{-1};
public:
    FDRedirect(const char* out_path, const char* err_path=nullptr) {
        saved_out = dup(STDOUT_FILENO);
        saved_err = dup(STDERR_FILENO);
        out_fd = ::open(out_path, O_WRONLY|O_CREAT|O_APPEND, 0644);
        if (!err_path) err_path = out_path;
        err_fd = ::open(err_path, O_WRONLY|O_CREAT|O_APPEND, 0644);
        dup2(out_fd, STDOUT_FILENO);
        dup2(err_fd, STDERR_FILENO);
    }
    ~FDRedirect() {
        fsync(STDOUT_FILENO);
        fsync(STDERR_FILENO);
        dup2(saved_out, STDOUT_FILENO); close(saved_out); close(out_fd);
        dup2(saved_err, STDERR_FILENO); close(saved_err); close(err_fd);
    }
};
