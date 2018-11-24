import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

class Figure:
    default_figsize_inches=(8,8)
    motion_notify_throttle_timer_interval = 100
    left_button = 1
    middle_button = 2
    right_button = 3

    def __init__(self):
        self.mouse_down_handler = None
        self.mouse_up_handler = None
        self.mouse_over_handler = None
        self.mouse_wheel_handler = None
        self.key_down_handler = None
        self.key_up_handler = None

        fig,ax = plt.subplots(figsize=Figure.default_figsize_inches)
        self.fig = fig
        self.ax = ax

        fig.canvas.mpl_connect('button_press_event', self.on_button_press)
        fig.canvas.mpl_connect('button_release_event', self.on_button_release)
        fig.canvas.mpl_connect('motion_notify_event', self.on_motion_notify)
        fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        fig.canvas.mpl_connect('key_release_event', self.on_key_release)

        ax.set_axis_off()

        motion_notify_throttle_timer = fig.canvas.new_timer(interval=self.motion_notify_throttle_timer_interval)
        motion_notify_throttle_timer.add_callback(self.on_motion_notify_throttled)
        motion_notify_throttle_timer.single_shot = True

        self.motion_notify_throttle_timer = motion_notify_throttle_timer
        self.motion_notify_throttle_timer_started = False
        self.motion_notify_args = None

        self.views = dict()
        self.last_used_view_id = 0

    def set_title(self, title):
        self.fig.suptitle(title)

    def set_image(self, img):
        self.img = self.ax.imshow(img, cmap='gray')
        self.ax.autoscale(False)

    def set_mouse_down_handler(self, handler):
        self.mouse_down_handler = handler

    def set_mouse_up_handler(self, handler):
        self.mouse_up_handler = handler

    def set_mouse_over_handler(self, handler):
        self.mouse_over_handler = handler

    def set_mouse_wheel_handler(self, handler):
        self.mouse_wheel_handler = handler

    def set_key_down_handler(self, handler):
        self.key_down_handler = handler

    def set_key_up_handler(self, handler):
        self.key_up_handler = handler

    def create_mask_view(self, colors, fill_value):
        levels = range(len(colors))
        cmap, norm = mpl.colors.from_levels_and_colors(levels, colors, extend='max')
        mask_size = self.img.get_size()
        mask_data = np.full(mask_size, fill_value, dtype='uint8')
        mask_view = self.ax.imshow(mask_data, cmap=cmap, norm=norm)
        mask_view_id = self.register_view(mask_view)
        return mask_view_id

    def set_mask_view_data(self, mask_view_id, mask_data, colors = None):
        if mask_view_id in self.views:
            mask_view = self.views[mask_view_id]
            mask_view.set_data(mask_data)
            if colors!=None:
                levels = range(len(colors))
                cmap, norm = mpl.colors.from_levels_and_colors(levels, colors, extend='max')
                mask_view.set_cmap(cmap)


    def create_circle_view(self, x, y, radius, color):
        circle_view = plt.Circle((x, y), radius, color=color, fill=False)
        self.ax.add_artist(circle_view)
        circle_view_id = self.register_view(circle_view)
        return circle_view_id

    def set_circle_view_data(self, circle_view_id, x, y, radius):
        if circle_view_id in self.views:
            circle_view = self.views[circle_view_id]
            circle_view.center = (x,y)
            circle_view.radius = radius

    def create_polyline_view(self, color):
        polyline_view = plt.Line2D((), (), color=color)
        self.ax.add_artist(polyline_view)
        polyline_view_id = self.register_view(polyline_view)
        return polyline_view_id

    def set_polyline_view_data(self, polyline_view_id, xs, ys):
        if polyline_view_id in self.views:
            polyline_view = self.views[polyline_view_id]
            polyline_view.set_data(xs, ys)

    def remove_view(self, view_id):
        if view_id in self.views:
            view = self.views[view_id]
            view.remove()
            del self.views[view_id]

    def remove_views(self, view_ids):
        for view_id in view_ids:
            if view_id in self.views:
                view = self.views[view_id]
                view.remove()
                del self.views[view_id]

    def draw(self):
        self.fig.canvas.draw()

    def show(self):
        plt.show()

    def register_view(self, view):
        new_view_id = self.last_used_view_id + 1
        self.views[new_view_id] = view
        self.last_used_view_id = new_view_id
        return new_view_id

    def raise_mouse_down(self, e):
        if not self.mouse_down_handler is None:
            self.mouse_down_handler(e)

    def raise_mouse_up(self, e):
        if not self.mouse_up_handler is None:
            self.mouse_up_handler(e)

    def raise_mouse_over(self, e):
        if not self.mouse_over_handler is None:
            self.mouse_over_handler(e)

    def raise_mouse_wheel(self, e):
        if not self.mouse_wheel_handler is None:
            self.mouse_wheel_handler(e)

    def raise_key_down(self, e):
        if not self.key_down_handler is None:
            self.key_down_handler(e)

    def raise_key_up(self, e):
        if not self.key_up_handler is None:
            self.key_up_handler(e)

    def is_over_image(self, e):
        return (not self.img is None) and e.inaxes == self.img.axes

    def on_button_press(self, e):
        if (self.is_over_image(e)):
            self.raise_mouse_down(e)

    def on_button_release(self, e):
        if (self.is_over_image(e)):
            self.raise_mouse_up(e)

    def on_motion_notify(self, e):
        if self.is_over_image(e):
            self.motion_notify_args = e
            if not self.motion_notify_throttle_timer_started:
                self.motion_notify_throttle_timer_started = True
                self.motion_notify_throttle_timer.start()

    def on_motion_notify_throttled(self):
        self.raise_mouse_over(self.motion_notify_args)
        self.motion_notify_throttle_timer_started = False

    def on_scroll(self, e):
        if (self.is_over_image(e)):
            self.raise_mouse_wheel(e)

    def on_key_press(self, e):
        self.raise_key_down(e)

    def on_key_release(self, e):
        self.raise_key_up(e)

class MaskBuilder:
    def __init__(self, num_rows, num_cols, fill_value):
        self.mask = np.full((num_rows, num_cols), fill_value, dtype='uint8')
        self.meshgrid = np.ogrid[:num_rows, :num_cols]

    def add_point(self, x, y, fill_value):
        self.mask[y, x] = fill_value

    def add_disk(self, x, y, radius, fill_value):
        mask = self.mask
        ys, xs = self.meshgrid
        disk_mask = np.less_equal(np.square(xs - x) + np.square(ys - y), radius * radius)
        mask[disk_mask] = fill_value

    def get_value_at(self, x, y):
        return self.mask[y,x]

    def get_mask(self):
        return self.mask

class GraphCutsPresenter:
    bgr_seed_color = (1, 0, 0, 1)
    obj_seed_color = (0, 0, 1, 1)
    bgr_label_color = (1, 0, 0, 0.5)
    obj_label_color = (0, 0, 1, 0.5)
    none_color = (0, 0, 0, 0)
    circle_around_cursor_color = (0, 0, 0, 1)

    disk_radius = 10
    max_disk_radius = 25
    min_disk_radius = 3

    def __init__(self, img, alg):
        self.img = img
        self.alg = alg
        self.bgr_value = alg.bgr_value
        self.obj_value = alg.obj_value
        self.none_value = alg.none_value
        self.circle_view_ids = []

        num_rows = self.img.shape[0]
        num_cols = self.img.shape[1]

        self.label_mask_builder = MaskBuilder(num_rows, num_cols, self.none_value)
        self.seed_mask_builder = MaskBuilder(num_rows, num_cols, self.none_value)

    def connect_figure(self, fig):
        self.fig = fig
        fig.set_title('Graph Cuts')
        fig.set_image(self.img)
        fig.set_mouse_down_handler(self.on_mouse_down)
        fig.set_mouse_up_handler(self.on_mouse_up)
        fig.set_mouse_over_handler(self.on_mouse_over)
        fig.set_mouse_wheel_handler(self.on_mouse_wheel)

        label_mask_colors = (self.bgr_label_color, self.obj_label_color, self.none_color)
        self.label_mask_view_id = fig.create_mask_view(label_mask_colors, self.none_value)

        seed_mask_colors = (self.bgr_seed_color, self.obj_seed_color, self.none_color)
        self.seed_mask_view_id = fig.create_mask_view(seed_mask_colors, self.none_value)

        self.mask_value = self.none_value

        self.circle_around_cursor_view = fig.create_circle_view(-self.disk_radius, -self.disk_radius, self.disk_radius, self.circle_around_cursor_color)

        self.fig.draw()

    def on_mouse_down(self, e):
        sx = int(e.xdata)
        sy = int(e.ydata)

        if (e.button == Figure.left_button):
            self.mask_value = self.bgr_value
            self.add_disk_to_seed_mask(sx, sy, self.mask_value)
        elif (e.button == Figure.right_button):
            self.mask_value = self.obj_value
            self.add_disk_to_seed_mask(sx, sy, self.mask_value)

    def on_mouse_up(self, e):
        mask_value = self.mask_value
        if mask_value == self.bgr_value or mask_value == self.obj_value:
            seed_mask = self.seed_mask_builder.get_mask()
            label_mask = self.alg.compute_labels(seed_mask)

            self.fig.remove_views(self.circle_view_ids)
            self.circle_view_ids = []

            self.fig.set_mask_view_data(self.label_mask_view_id, label_mask)
            self.fig.set_mask_view_data(self.seed_mask_view_id, seed_mask)
            self.mask_value = self.none_value
            self.fig.draw()

    def on_mouse_over(self, e):
        mask_value = self.mask_value
        sx = int(e.xdata)
        sy = int(e.ydata)
        
        if mask_value == self.bgr_value or mask_value == self.obj_value:
            self.fig.set_circle_view_data(self.circle_around_cursor_view, -100, -100, 10)
            self.add_disk_to_seed_mask(sx, sy, self.mask_value)
        else:
            self.set_circle_around_cursor(sx, sy)

    def on_mouse_wheel(self, e):
        sx = int(e.xdata)
        sy = int(e.ydata)

        if e.button == 'up':
            if self.disk_radius < self.max_disk_radius:
                self.disk_radius = self.disk_radius + 1
                self.set_circle_around_cursor(sx, sy)
                
        elif e.button == 'down':
            if self.disk_radius > self.min_disk_radius:
                self.disk_radius = self.disk_radius - 1
                self.set_circle_around_cursor(sx, sy)
            
    def add_disk_to_seed_mask(self, sx, sy, mask_value):
        self.seed_mask_builder.add_disk(sx, sy, self.disk_radius, mask_value)
        circle_view_id = self.fig.create_circle_view(sx, sy, self.disk_radius, self.get_circle_color(mask_value))
        self.circle_view_ids.append(circle_view_id)
        self.fig.draw()

    def get_circle_color(self, mask_value):
        if mask_value == self.bgr_value:
            return self.bgr_seed_color
        elif mask_value == self.obj_value:
            return self.obj_seed_color
        else:
            return self.none_color

    def set_circle_around_cursor(self, sx, sy):
        self.fig.set_circle_view_data(self.circle_around_cursor_view, sx, sy, self.disk_radius)
        self.fig.draw()

class PolylineBuilder:
    def __init__(self):
        self.xs = []
        self.ys = []

    def add_points(self, xs, ys):
        self.xs.extend(xs)
        self.ys.extend(ys)

    def empty(self):
        return len(self.xs) == 0

    def get_starting_point(self):
        sx = self.xs[0]
        sy = self.ys[0]
        return sx,sy

    def get_ending_point(self):
        sx = self.xs[-1]
        sy = self.ys[-1]
        return sx,sy

    def get_polyline_xs(self):
        return self.xs

    def get_polyline_ys(self):
        return self.ys

class LiveWirePresenter:
    contour_polyline_color = 'blue'
    path_polyline_color = 'red'

    def __init__(self, img, alg):
        self.img = img
        self.alg = alg
        
        self.contour_builder = PolylineBuilder()

    def connect_figure(self, fig):
        self.fig = fig
        fig.set_title('Live Wire')
        fig.set_image(self.img)
        fig.set_mouse_down_handler(self.on_mouse_down)
        fig.set_mouse_over_handler(self.on_mouse_over)

        self.contour_polyline_view_id = fig.create_polyline_view(self.contour_polyline_color)
        self.path_polyline_view_id = fig.create_polyline_view(self.path_polyline_color)
        self.contour_polyline_closed = False

        self.fig.draw()

    def on_mouse_down(self, e):
        if self.contour_polyline_closed:
            pass
        elif e.button == Figure.left_button:
            sx = int(e.xdata)
            sy = int(e.ydata)

            if self.alg.paths_computed():
                xs, ys = self.alg.get_path_to(sx, sy)
                self.add_points_to_contour_polyline(reversed(xs), reversed(ys))
                self.fig.set_polyline_view_data(self.path_polyline_view_id, (), ())
                self.fig.draw()

            self.alg.compute_paths_starting_at(sx, sy)
        elif e.button == Figure.right_button and not self.contour_builder.empty():
            sx = int(e.xdata)
            sy = int(e.ydata)

            xs,ys = self.alg.get_path_to(sx, sy)
            self.add_points_to_contour_polyline(reversed(xs), reversed(ys))
        
            self.alg.compute_paths_starting_at(sx, sy)

            tx,ty = self.contour_builder.get_starting_point()
            xs,ys = self.alg.get_path_to(tx, ty)
            self.add_points_to_contour_polyline(reversed(xs), reversed(ys))
            self.fig.set_polyline_view_data(self.path_polyline_view_id, (), ())
            self.fig.draw()
            self.contour_polyline_closed = True

    def on_mouse_over(self, e):
        if not self.contour_polyline_closed and self.alg.paths_computed():
            tx = int(e.xdata)
            ty = int(e.ydata)
            xs, ys = self.alg.get_path_to(tx, ty)
            self.fig.set_polyline_view_data(self.path_polyline_view_id, xs, ys)
            self.fig.draw()

    def add_points_to_contour_polyline(self, xs, ys):
        self.contour_builder.add_points(xs, ys)
        contour_polyline_xs = self.contour_builder.get_polyline_xs()
        contour_polyline_ys = self.contour_builder.get_polyline_ys()
        self.fig.set_polyline_view_data(self.contour_polyline_view_id, contour_polyline_xs, contour_polyline_ys)

class RegionGrowingPresenter:
    def __init__(self, img, alg):
        self.img = img
        self.alg = alg
        self.colors = ('none','red','green','blue','orange','pink','lightblue','lightgreen')
        self.none_value = alg.none_value

    def connect_figure(self, fig):
        self.fig = fig
        fig.set_title('Region Growing')
        fig.set_image(self.img)
        fig.set_mouse_down_handler(self.on_mouse_down)

        self.region_mask_view_id = fig.create_mask_view(self.colors, self.none_value)
        self.fig.draw()

    def on_mouse_down(self, e):
        if e.button == Figure.left_button:
            sx = int(e.xdata)
            sy = int(e.ydata)
            self.alg.grow_new_region_starting_at(sx, sy)
            region_mask = self.alg.get_region_mask()
            self.fig.set_mask_view_data(self.region_mask_view_id, region_mask)
            self.fig.draw()

from matplotlib.colors import colorConverter, cnames

class KmeansPresenter:
    def __init__(self, img, alg):
        self.img = img
        self.alg = alg
        self.colors = ['none']     # example of a color list is ['red','green','blue',(1.0,0.0,0.0,0.3)]
        self.palette = 'none'  # other possible values are 'transparent' 'solid', and 'mean'
        self.set_colors('transparent')
 
    def connect_figure(self, fig):
        self.fig = fig
        fig.set_title('K-means')
        fig.set_image(self.img)

        fig.set_key_down_handler(self.on_key_down)
        self.region_mask_view_id = fig.create_mask_view(self.colors, self.alg.no_label)
        self.fig.draw()
        
    def set_colors(self, palette):
        self.palette = palette
        all_colors = cnames.keys()
        if palette=='transparent':
            self.colors = list(colorConverter.to_rgba_array(np.random.choice(all_colors,self.alg.k),0.7)) + ['none']
        if palette=='solid':
            self.colors = list(np.random.choice(all_colors,self.alg.k)) + ['none']
        if palette=='mean':
            self.colors = list(np.minimum(1.0,self.alg.means[c,:3]/255.0) for c in range(self.alg.k)) + ['none']
        
    def on_key_down(self, e): # do the following when certain "key" is pressed
        if e.key == 'i':    # make 1 iteration of K-means
            self.alg.compute_k_means_clusters()
            region_mask = self.alg.get_region_mask()
            if (self.palette=='mean'):
                self.set_colors('mean')
            self.fig.set_mask_view_data(self.region_mask_view_id, region_mask, self.colors)
            self.fig.draw() 
        if e.key == 'c':    # run K-means to convergence (when SSE (energy) improvement is less than given threshold)
            delta = np.infty
            while (delta>0.001):
                delta = self.alg.compute_k_means_clusters()
            region_mask = self.alg.get_region_mask()
            if (self.palette=='mean'):
                self.set_colors('mean')
            self.fig.set_mask_view_data(self.region_mask_view_id, region_mask, self.colors)
            self.fig.draw() 
        if e.key == 'v':    # run K-means to convergence with visualization of each iteration
            delta = float('inf')
            while (delta>0.001):
                delta = self.alg.compute_k_means_clusters()
                region_mask = self.alg.get_region_mask()
                if (self.palette=='mean'):
                    self.set_colors('mean')
                self.fig.set_mask_view_data(self.region_mask_view_id, region_mask, self.colors)
                self.fig.draw() 
        if e.key == 'r':    # restart K-means (with random clusters)
            self.alg.init_means()
            self.alg.compute_k_means_clusters()
            region_mask = self.alg.get_region_mask()
            if (self.palette=='mean'):
                self.set_colors('mean')
            self.fig.set_mask_view_data(self.region_mask_view_id, region_mask, self.colors)
            self.fig.draw() 
        if e.key == 's':    # changes colors to some Solid colors (randomly selected)
            self.set_colors('solid')
            region_mask = self.alg.get_region_mask()
            self.fig.set_mask_view_data(self.region_mask_view_id, region_mask, self.colors)
            self.fig.draw() 
        if e.key == 't':    # change colors to some Transparent colors (randomly selected)
            region_mask = self.alg.get_region_mask()
            self.set_colors('transparent')
            self.fig.set_mask_view_data(self.region_mask_view_id, region_mask, self.colors)
            self.fig.draw() 
        if e.key == 'm':    # changes colors to mean RGB in each segment
            region_mask = self.alg.get_region_mask()
            self.set_colors('mean')
            self.fig.set_mask_view_data(self.region_mask_view_id, region_mask, self.colors)
            self.fig.draw() 