//
// Created by alex on 7/4/2024.
//

#ifndef EVENT_HANDLER_H
#define EVENT_HANDLER_H
#include <memory>

namespace naomi
{
class spacecraft;
namespace events
{
class event_handler
{
public:
  virtual ~event_handler() = default;
  virtual void handle_event(const std::shared_ptr<spacecraft>& sc, double t) = 0;
};
}
}
#endif //EVENT_HANDLER_H
