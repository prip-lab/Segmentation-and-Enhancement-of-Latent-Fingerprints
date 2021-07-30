import glob
import numpy as np
import matplotlib.pyplot as plt
from skimage.restoration import denoise_nl_means, estimate_sigma
import skimage.io
from skimage import exposure
from skimage.morphology import disk
from skimage.filters import rank
from scipy.ndimage import gaussian_filter
from skimage import data
from skimage import img_as_float
from skimage.morphology import reconstruction
import math
import cv2
# import sklearn.linear_model
import sys
sys.path.append('OF/')
import pdb
#import get_maps
from skimage.filters import gaussian
import os
from skimage import io
from skimage.transform import rescale
def nextpow2(x):
    return int(math.ceil(math.log(x, 2)))


def LowpassFiltering(img,L):

    h,w = img.shape
    h2,w2 = L.shape

    img = cv2.copyMakeBorder(img, 0, h2-h, 0, w2-w, cv2.BORDER_CONSTANT, value=0)

    img_fft = np.fft.fft2(img)
    img_fft = np.fft.fftshift(img_fft)
    #img_fft = fftshift(fft2(img,h2,w2));

    img_fft = img_fft * L

    rec_img = np.fft.ifft2(np.fft.fftshift(img_fft))
    rec_img = np.real(rec_img)
    rec_img = rec_img[:h,:w]
    #plt.imshow(rec_img, cmap='gray')
    #plt.show()
    return rec_img





def FastCartoonTexture(img,sigma=2.5,show=False):
    img = img.astype(np.float32)
    h, w = img.shape
    h2 = 2 ** nextpow2(h)
    w2 = 2 ** nextpow2(w)

    FFTsize = np.max([h2, w2])
    x, y = np.meshgrid(range(-int(FFTsize / 2), int(FFTsize / 2)), range(-int(FFTsize / 2), int(FFTsize / 2)))
    r = np.sqrt(x * x + y * y) + 0.0001
    r = r/FFTsize

    L = 1. / (1 + (2 * math.pi * r * sigma)** 4)
    img_low = LowpassFiltering(img, L)

    gradim1=  compute_gradient_norm(img)
    gradim1 = LowpassFiltering(gradim1,L)

    gradim2=  compute_gradient_norm(img_low)
    gradim2 = LowpassFiltering(gradim2,L)

    diff = gradim1-gradim2
    ar1 = np.abs(gradim1)
    diff[ar1>1] = diff[ar1>1]/ar1[ar1>1]
    diff[ar1 <= 1] = 0

    cmin = 0.3
    cmax = 0.7

    weight = (diff-cmin)/(cmax-cmin)
    weight[diff<cmin] = 0
    weight[diff>cmax] = 1


    u = weight * img_low + (1-weight)* img

    temp = img - u

    lim = 20

    temp1 = (temp + lim) * 255 / (2 * lim)

    temp1[temp1 < 0] = 0
    temp1[temp1 >255] = 255
    v = temp1
    if show:
        plt.imshow(v,cmap='gray')
        plt.show()
    return u, v, temp

def compute_gradient_norm(input):
    input = input.astype(np.float32)
    h,w = input.shape

    Gx, Gy = np.gradient(input)
    out = np.sqrt(Gx * Gx + Gy * Gy) + 0.000001
    #input = double(input);
    #[Gx, Gy] = gradient(input);
    #out = sqrt(Gx.*Gx + Gy.*Gy);
    return out



if __name__=='__main__':
    # construct ridge structure dictionary for quality estimation or ridge spacing estimation
    # img = io.imread(')

    output_path = './results/'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    # file = '/home/kaicao/Research/Data/Fingerprints/Rolled/NIST4/f0062.bmp'
    file = 'Female_right_palm_print_-_1.jpg'

    fname = os.path.basename(file)[:-4]
    img = io.imread(file)
    img = img.astype(np.float32)

    img = rescale(img, 0.25, anti_aliasing=False)

    # 
    if len(img.shape) == 3:
        img = np.mean(img,axis=2)

    # 
    u, v, tmp = FastCartoonTexture(img, sigma=3)
    u = np.clip(u, 0, 255)
    v = np.clip(v, 0, 255)
    
    tmp = tmp - np.min(tmp)
    tmp = np.clip(tmp, 0, 255)

    cartoon = u.astype(np.uint8)
    texture = v.astype(np.uint8)
    enhance_texture = tmp.astype(np.uint8)

    cartoon_file = output_path + fname + '_cartoon.bmp'
    original_file = output_path + fname + '.bmp'
    texture_file = output_path + fname + '_texture.bmp'

    enh_texture_file = 'results/' + fname + '_enh_texture.bmp'

    img = img.astype(np.uint8)
    io.imsave(original_file, img)
    io.imsave(cartoon_file, cartoon)
    io.imsave(texture_file, texture)
    io.imsave(enh_texture_file, enhance_texture)

    # plt.subplot(1,3,1)
    # plt.imshow(img, cmap='gray')
    # plt.title('Original image')
    # plt.subplot(1,3,2)
    # plt.imshow(u, cmap='gray')
    # plt.title('Cartoon image')
    # plt.subplot(1,3,3)
    # plt.imshow(v, cmap='gray')
    # plt.title('Texture image')
    # plt.show()

    # root_path = '/home/kaicao/Research/Data/Rolled/'
    # img_ext = ['jpg', 'bmp', 'png']

    # for root, dirs, files in os.walk(root_path, topdown=False):
    #     # pdb.set_trace()
    #     files = [file for file in files if file[-3:] in img_ext]
    #     out_path = root.replace('Rolled', 'Rolled_Texture')
    #     if not os.path.exists(out_path):
    #         os.makedirs(out_path)
    #     for file in files:
    #         img = io.imread(os.path.join(root,file))
    #         img = FastCartoonTexture(img,sigma=2.5,show=False)
    #         out_file = os.path.join(out_path,file) 
    #         io.imsave(out_file, img)
    #         # pdb.set_trace()

     

