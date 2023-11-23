#!/usr/bin/env python
"""Visualization of the images."""

import numpy as np
import matplotlib.pyplot as plt

def image_background(
    img_region=None,
    img=None,
    cmap: str='binary_r',
    alpha=1,
    linewidths=1,
    colors='cyan',
    kws_region={},
    show_scalebar=None,
    scalebar_color: str='w',
    show_cbar: bool=True,
    test=False,
    ax=None,
    **kws_img,
    ) -> plt.Axes:
    """Plot the image, to be used as a background to the annotations.

    Args:
        img_region (_type_, optional): segmentation image. Defaults to None.
        img (_type_, optional): image with intensity values. Defaults to None.
        cmap (str, optional): colormap name. Defaults to 'binary_r'.
        alpha (int, optional): transparency. Defaults to 1.
        linewidths (int, optional): segmentation contour line width. Defaults to 1.
        colors (str, optional): color of the segmentation line. Defaults to 'cyan'.
        kws_region (dict, optional): parameters provided to the segmentation plot. Defaults to {}.
        show_scalebar (_type_, optional): show scale bar. Defaults to None.
        scalebar_color (str, optional): color of the scale bar. Defaults to 'w'.
        show_cbar (bool, optional): show colorbar. Defaults to True.
        test (bool, optional): test-mode. Defaults to False.
        ax (_type_, optional): subplot object. Defaults to None.

    Keyword Args:
        parameters provided to the `plt.imshow`.
    Returns:
        plt.Axes
    """
    if cmap=='gfp':
        from roux.viz.colors import make_cmap
        cmap=make_cmap(
            ["#000000",'#83f52c'],
            N=50,
            )        
    ax=plt.subplot(111) if ax is None else ax
    # TODO preprocess
    # if not rotation is None:
    #     assert rotation%90 ==0:
    # img=numpy.rot90(img,3)
    # img_region=numpy.rot90(img_region,3)

    if not img is None:
        ax_img=ax.imshow(
            img,
            cmap=cmap,
            alpha=alpha,
            **kws_img,
            )
    if not img_region is None:
        ax.contour(
            img_region, [0.5], 
            linewidths=linewidths,
            linestyles='dashed',
            colors=colors,
            **kws_region,
            )
    if not show_scalebar is None:
        ax=annotate_scalebar(
            ax,
            img=img,
            pixels=list(show_scalebar.values())[0],
            label=list(show_scalebar.keys())[0],
            )
    if not img is None and show_cbar:
        cax = ax.figure.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height if ax.figure.get_size_inches()[1]<5 else (2/ax.figure.get_size_inches()[1])])
        cbar = ax.figure.colorbar(ax_img,cax=cax)
        cbar.ax.set_ylabel('Intensity')
        # cbar.ax.hist(img.ravel()) # show histogram in the colorbar
    ax.grid(test)
    ax.set_aspect('equal')
    return ax

def annot_cells(
    label_image,
    show_boxes: bool=False,
    ax: plt.Axes=None,
    )-> plt.Axes:
    """Annotate the cells on an image.

    Args:
        label_image (_type_): image with the labeled regions 
        show_boxes (bool, optional): show boxes around regions. Defaults to False.
        ax (plt.Axes, optional): plt.Axes. Defaults to None.

    Returns:
        plt.Axes
    """
    from skimage.measure import regionprops
    import matplotlib.patches as mpatches
    for region in regionprops(label_image):
        minr, minc, maxr, maxc = region.bbox
        ax.text(x=np.mean([minc,maxc]),y=np.mean([minr,maxr]),s=f"{region.label:.0f}",ha='center',va='center',color='r',zorder=10,size=10)
        if show_boxes:
            # draw rectangle around segmented coins
            rect = mpatches.Rectangle((minc, minr), maxc - minc, maxr - minr,
                                      fill=False, edgecolor='red', linewidth=2)
            ax.add_patch(rect)
    return ax

def image_regions_annotated(
    img_region,
    img,
    show_boxes: bool=False,
    **kws_img,
    ) -> plt.Axes:
    """
    Image with the annotated regions.
    Usage: for QC of the segmentation.

    Args:
        img_region (_type_): image with segmentated regions.
        img (_type_): image with intensity.
        show_boxes (bool, optional): whether to show the boxes around the regions. Defaults to False.

    Keyword Args:
        parameters provided to the `image_background` function.

    Returns:
        plt.Axes
    """
    ax=image_background(
        img_region=img_region,
        img=img,
        **kws_img,
        )
    annot_cells(
        label_image=img_region,
        show_boxes=show_boxes,
        ax=ax,
        )
    return ax