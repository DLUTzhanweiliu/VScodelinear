# import matplotlib.pyplot as plt
# import matplotlib.cbook as cbook
# import numpy as np


# # Data are 256x256 16 bit integers.
# with cbook.get_sample_data('s1045.ima.gz') as dfile:
#     im = np.frombuffer(dfile.read(), np.uint16).reshape((256, 256))
#     for i in range(256):
#         for j in range(256):
#             print(im[i,j])

# fig, ax = plt.subplots(num="MRI_demo")
# ax.imshow(im, cmap="gray")
# # ax.axis('off')
# ax.plot([1,2,100],[1,2,100],color='r')

# delta = 2
# x = np.arange(0.0, 200.0, delta)
# y = np.arange(0.0, 200.0, delta)
# X, Y = np.meshgrid(x, y)
# Z1 = np.exp(-X**2 - Y**2)
# Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2)
# Z = (Z1 - Z2) * 2

# # fig, ax = plt.subplots()
# ax.contour(X, Y, Z, zorder=10)
# # ax.clabel(CS, inline=1, fontsize=10)
# # ax.set_title('Simplest default with labels')

# plt.show()
if True:
    raise Exception("你好")
    print(1)
