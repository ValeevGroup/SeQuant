// sequant-schedule: validate, diff and render JSON schedule documents.
//
//   sequant-schedule validate <doc>
//   sequant-schedule diff <a> <b>
//   sequant-schedule view <doc> [--run <run.json>] [-o out.html]
//                             [--template <schedule-view.html>]
//
// exit codes: 0 ok, 1 differences / invalid input, 2 usage

#include <SeQuant/core/eval/schedule_document.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#ifndef SEQUANT_SCHEDULE_VIEW_TEMPLATE
#define SEQUANT_SCHEDULE_VIEW_TEMPLATE ""
#endif

namespace {

namespace sched = sequant::eval::schedule;

int usage() {
  std::cerr << "usage:\n"
               "  sequant-schedule validate <doc>\n"
               "  sequant-schedule diff <a> <b>\n"
               "  sequant-schedule view <doc> [--run <run.json>] [-o out.html]"
               " [--template <file>]\n";
  return 2;
}

bool slurp(std::string const& path, std::string& out) {
  std::ifstream in(path, std::ios::binary);
  if (!in) return false;
  out.assign(std::istreambuf_iterator<char>(in), {});
  return true;
}

// reader errors are prefixed with the path, so `diff a b` names the file
sched::Document load(std::string const& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw sequant::Exception("cannot open " + path);
  try {
    return sched::read(in);
  } catch (std::exception const& e) {
    throw sequant::Exception(path + ": " + e.what());
  }
}

int validate(std::string const& path) {
  std::cout << sched::summary_line(load(path)) << "\n";
  return 0;
}

int diff(std::string const& pa, std::string const& pb) {
  auto const a = load(pa);
  auto const b = load(pb);
  auto const diffs = sched::compare(a, b);
  for (auto const& d : diffs)
    std::cout << d.path << ": " << d.left << " vs " << d.right << "\n";
  return diffs.empty() ? 0 : 1;
}

int view(std::vector<std::string> const& args) {
  std::string doc_path, run_path, out_path;
  std::string tmpl_path = SEQUANT_SCHEDULE_VIEW_TEMPLATE;
  for (std::size_t i = 0; i < args.size(); ++i) {
    auto const& a = args[i];
    auto next = [&](std::string& dst) {
      if (i + 1 >= args.size()) return false;
      dst = args[++i];
      return true;
    };
    if (a == "--run") {
      if (!next(run_path)) return usage();
    } else if (a == "-o") {
      if (!next(out_path)) return usage();
    } else if (a == "--template") {
      if (!next(tmpl_path)) return usage();
    } else if (doc_path.empty() && !a.empty() && a[0] != '-') {
      doc_path = a;
    } else {
      return usage();
    }
  }
  if (doc_path.empty()) return usage();
  if (out_path.empty()) out_path = doc_path + ".html";

  auto const doc = load(doc_path);  // validates
  auto const doc_json = sched::to_json(doc);

  // the run file is a JSON `run` section; it is validated and re-serialized
  std::string run_json = "null";
  if (!run_path.empty()) {
    std::string text;
    if (!slurp(run_path, text)) {
      std::cerr << "cannot open " << run_path << "\n";
      return 1;
    }
    try {
      run_json = sched::to_json(sched::run_from_json(text));
    } catch (std::exception const& e) {
      throw sequant::Exception(run_path + ": " + e.what());
    }
  }

  std::string page;
  if (!slurp(tmpl_path, page)) {
    std::cerr << "viewer template not found: " << tmpl_path << "\n";
    return 1;
  }
  // Both texts go inside <script> elements, hence the escaping. The slots
  // are located in the template before anything is inserted, so neither
  // inserted text can pass for the other's slot, and they are filled from the
  // back so the earlier offset stays valid.
  constexpr std::string_view doc_slot = "/*__SCHEDULE_JSON__*/";
  constexpr std::string_view run_slot = "/*__RUN_JSON__*/";
  auto const doc_pos = page.find(doc_slot);
  auto const run_pos = page.find(run_slot);
  if (doc_pos == std::string::npos) {
    std::cerr << "template lacks the " << doc_slot << " slot: " << tmpl_path
              << "\n";
    return 1;
  }
  if (run_pos == std::string::npos && !run_path.empty()) {
    std::cerr << "template lacks the " << run_slot
              << " slot, so --run cannot be embedded: " << tmpl_path << "\n";
    return 1;
  }
  std::vector<std::pair<std::size_t, std::pair<std::string_view, std::string>>>
      fills{{doc_pos, {doc_slot, sched::escape_for_html_script(doc_json)}}};
  if (run_pos != std::string::npos)
    fills.emplace_back(
        run_pos, std::pair{run_slot, sched::escape_for_html_script(run_json)});
  std::sort(fills.begin(), fills.end(),
            [](auto const& a, auto const& b) { return a.first > b.first; });
  for (auto const& [pos, fill] : fills)
    page.replace(pos, fill.first.size(), fill.second);
  std::ofstream out(out_path, std::ios::binary);
  if (!out) {
    std::cerr << "cannot write " << out_path << "\n";
    return 1;
  }
  out << page;
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> args(argv + 1, argv + argc);
  if (args.empty()) return usage();
  try {
    auto const& verb = args[0];
    if (verb == "validate" && args.size() == 2) return validate(args[1]);
    if (verb == "diff" && args.size() == 3) return diff(args[1], args[2]);
    if (verb == "view")
      return view(std::vector<std::string>(args.begin() + 1, args.end()));
    return usage();
  } catch (std::exception const& e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
