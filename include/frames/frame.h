//
// Created by alex_ on 7/24/2024.
//

#ifndef FRAME_H
#define FRAME_H

namespace naomi::frames
{
class frame
{
public:
  virtual state_type get_transform_to(const frame& other) = 0;
};
}
#endif //FRAME_H
