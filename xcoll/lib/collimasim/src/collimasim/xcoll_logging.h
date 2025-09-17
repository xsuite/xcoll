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

// ------------------------------ //
// Logging functionality for Root //
// ------------------------------ //

namespace rootlog {
  inline std::unique_ptr<std::ofstream> info;
  inline std::unique_ptr<std::ofstream> err;
  inline ErrorHandlerFunc_t prev = nullptr;
  inline std::once_flag init_once;

  inline void Handler(int level, Bool_t abort, const char* loc, const char* msg) {
    std::ostream& os = (level >= kError) ? *err : *info;
    os << loc << ": " << msg << '\n';
    os.flush();                 // keep logs timely
    if (abort) ::abort();       // mimic ROOT's behaviour
  }

  inline void init(const std::string& infoPath,
                   const std::string& errPath,
                   bool append = true)
  {
    std::call_once(init_once, [&]{
      auto mode = append ? std::ios::app : std::ios::trunc;
      info = std::make_unique<std::ofstream>(infoPath, mode);
      err  = std::make_unique<std::ofstream>(errPath.empty() ? infoPath : errPath, mode);
      if (!info->is_open() || !err->is_open())
        throw std::runtime_error("Failed to open ROOT log files");
      prev = GetErrorHandler();
      SetErrorHandler(&Handler);
    });
  }

  inline void restore() {
    if (prev) SetErrorHandler(prev);   // optional: restore if you ever want to undo
  }
} // namespace rootlog


// -------------------------------- //
// Logging functionality for Geant4 //
// -------------------------------- //

class XcollLogFileSession : public G4UIsession {
public:
  XcollLogFileSession()
  : logfile("geant4.out"), errfile("geant4.err") {}

  virtual ~XcollLogFileSession() { logfile.close(); errfile.close();}

  G4int ReceiveG4cout(const G4String& msg) override {
    logfile << msg;
    return 0;
  }
  G4int ReceiveG4cerr(const G4String& msg) override {
    errfile << msg;
    return 0;
  }

private:
  std::ofstream logfile;
  std::ofstream errfile;
};

void RedirectGeant4() {
  auto ui = G4UImanager::GetUIpointer();
  static XcollLogFileSession session;
  ui->SetCoutDestination(&session);
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
