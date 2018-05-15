# -*- coding:utf-8 -*-

"""PyHdust *hdrpil* module: third-part functions for HDR process of images

import pyhdust.hdrpil asd hdr
proj = hdr.hdr(case='img04', img_type='png', cur_dir='./')
proj.get_hdr(strength=[1.],naturalness=[1.0,1.1,1.2,1.5])

:license: ?
"""
from __future__ import print_function
import os
from copy import copy
import warnings as _warn

try:
    from PIL import Image
    import pylab
    pylab_loaded = 1
except ImportError:
    pylab_loaded = 0
    _warn.warn('matplotlib+pylab and/or PIL/Pillow not installed!!!')

__author__ = "bpowah"
__email__ = "bpowah@gmail.com"


class hdr:

    """
    a collection of images to merege into HDR

    blend an arbitrary number of photos into a single HDR image
    or several images with various combinations of HDR parameters

    it is assumed that project folders contain all the images you want to merge
    case = folder name
    cur_dir = folder where the project folders are located
    images are auto-sorted by degree of exposure (image brightness)
    """

    def get_imgs(self):
        """
        load a list of images from folder
        sort them from lightest to darkest
        """
        drks = []
        imgs = [Image.open(os.path.join(self.indir, fn)) for fn in self.fns]
        for img in imgs:
            samp = img.resize((50, 50))  # crop(box=((10,10,20,20)))
            drks.append(sum(samp.getdata(band=0)))
        if self.resize is not None:
            newsize = tuple([int(x * self.resize) for x in img.size])
            imgs = [img.resize(newsize) for img in imgs]

        newdrks = sorted(drks)
        idxs = [drks.index(x) for x in newdrks]
        imgs = [imgs[x] for x in idxs]
        self.fns = [self.fns[x] for x in idxs]

        self.drks = drks
        self.idxs = idxs

        print('got', len(imgs), 'images')
        return imgs

    def __init__(
        self,
        case='',
        resize=None,
        img_type='.jpg',
        cur_dir=r'./' ):
        """
        load a project
        all images of [img_type] are loaded from folder [case]
        and resized by a factor [resize]
        """
        if not case:
            raise ValueError('case is required')
        self.resize = resize
        self.ext = img_type.strip('.')
        self.case = case
        self.cur_dir = os.path.normpath(cur_dir)
        self.indir = os.path.join(self.cur_dir, case)
        self.outdir = os.path.join(self.indir, 'out')

        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)
        f = lambda x: x.upper().endswith(self.ext.upper())
        self.fns = filter(f, os.listdir(self.indir))
        self.imgs = self.get_imgs()
        self.fin_set = {}

    def get_masks(self, imgs):
        """
        create a set of masks from a list of images
        (one mask for every adjacent pair of images
        """
        masks = []
        mask_ct = len(imgs) - 1
        imgs = [self.bal(img.convert(mode='L')) for img in imgs]
        for i in range(mask_ct):
            blend_fraction = .5  # 1. - (float(i)+.5)/float(mask_ct)
            m = Image.blend(imgs[i], imgs[i + 1], blend_fraction)
            masks.append(m)
            # print blend_fraction
        print('blending using', mask_ct, 'masks')
        return masks

    def lev(self, im):
        # im = ImageEnhance.Brightness(im).enhance(self.bri)
        # im = ImageEnhance.Contrast(im).enhance(self.con)
        return self.bal(im, self.str)

    def bal(self, im):
        """
        adjust the balance of the mask
        (re-distribute the histogram so that there are more
        extreme blacks and whites)
        like increasing the contrast, but without clipping
        and maintains overall average image brightness
        """
        h = im.histogram()
        ln = range(len(h))
        up = [sum(h[0: i]) for i in ln]
        lo = [sum(h[i:-1]) for i in ln]
        ct = sum(h)
        st = int(self.cur_str * 255.)

        lut = [i + st * up[i] * lo[i] * (up[i] - lo[i]) / ct**3 for i in ln]
        for i in ln:
            if lut[i] < 1:
                lut[i] = 1
            if lut[i] > 255:
                lut[i] = 255
        return im.point(lut)

    def save_im(self, im, name):
        """
        save an image
        """
        print('saving', name)
        im.save(os.path.join(self.outdir, name + '.jpg'), format='JPEG')

        if pylab_loaded:
            pylab.cla()
            h = im.histogram()
            rgb = h[0:256], h[256:256 * 2], h[256 * 2:]
            x = range(256)
            for b in rgb:
                pylab.plot(x, b)
            pylab.xlim(xmax=255)
            pylab.savefig(os.path.join(self.outdir, name + '_histogram.png'))

    def merge(self, imgs):
        """
        combine a set images into a smaller set by combinding all
        adjacent images
        """
        masks = self.get_masks(imgs)
        imx = lambda i: Image.composite(imgs[i], imgs[i + 1], masks[i])
        return [imx(i) for i in range(len(masks))]

    def merge_all(self, imgs):
        """
        iteratively merge a set of images until only one remains
        """
        while len(imgs) > 1:
            imgs = self.merge(imgs)
        return imgs[0]

    def get_hdr(self, strength=[0.0, 1.0, 2.0], naturalness=[0.8, 0.9, 1.0]):
        """
        process the hdr image(s)
        strength - a list or a float that defines how strong the hdr
                   effect should be
                 - a value of zero will combine images by using a
                   greyscale image average
                 - a value greater than zero will use higher contrast
                   versions of those greyscale images
        naturalness- values between zero and one
                 - zero will be a very high-contrast image
                 - 1.0 will be a very flat image
                 - 0.7 to 0.9 tend to give the best results
        """
        contrast = 1.0
        brightness = 1.0
        self.con = contrast  # not used
        self.bri = brightness  # not used
        self.nat = self.to_list(naturalness)
        self.str = self.to_list(strength)
        self.fin_set = {}

        for s in self.str:
            self.cur_str = s
            print('getting saturation image, strength', str(s))
            imgs = copy(self.imgs)
            sat_img = self.merge_all(imgs)

            print('getting contrast image')
            imgs.reverse()
            con_img = self.merge_all(imgs)

            self.final_blend(con_img, sat_img)

        self.print_set(self.fin_set)
        ori_set = dict(
            ('orig_' + str(i), self.imgs[i]) for i in range(len(self.imgs)))
        self.print_set(ori_set)
        print(self.indir)

    def final_blend(self, im1, im2):
        """
        combines a saturated image with a contrast image
        and puts them in a dictionary of completed images
        """
        for nat in self.nat:
            n_s = '_nat' + str(round(nat, 2))
            s_s = '_str' + str(round(self.cur_str, 2))
            s = self.case + s_s + n_s
            self.fin_set[s] = Image.blend(im1, im2, nat)

    def print_set(self, im_dict):
        """
        print all rendered images
        """
        for k in im_dict.keys():
            self.save_im(im_dict[k], k)

    def to_list(self, val):
        # if type(val) != type(list()):
        if not isinstance(val, list):
            val = [float(val)]
        else:
            val = [float(v) for v in val]
        return val

# MAIN ###
if __name__ == "__main__":
    pass
