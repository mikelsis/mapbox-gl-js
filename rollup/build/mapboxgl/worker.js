define(['./shared'], (function (ref_properties) { 'use strict';

//      
function stringify(obj) {
    if (typeof obj === 'number' || typeof obj === 'boolean' || typeof obj === 'string' || obj === undefined || obj === null)
        return JSON.stringify(obj);
    if (Array.isArray(obj)) {
        let str = '[';
        for (const val of obj) {
            str += `${ stringify(val) },`;
        }
        return `${ str }]`;
    }
    let str = '{';
    for (const key of Object.keys(obj).sort()) {
        str += `${ key }:${ stringify(obj[key]) },`;
    }
    return `${ str }}`;
}
function getKey(layer) {
    let key = '';
    for (const k of ref_properties.refProperties) {
        key += `/${ stringify(layer[k]) }`;
    }
    return key;
}
/**
 * Given an array of layers, return an array of arrays of layers where all
 * layers in each group have identical layout-affecting properties. These
 * are the properties that were formerly used by explicit `ref` mechanism
 * for layers: 'type', 'source', 'source-layer', 'minzoom', 'maxzoom',
 * 'filter', and 'layout'.
 *
 * The input is not modified. The output layers are references to the
 * input layers.
 *
 * @private
 * @param {Array<Layer>} layers
 * @param {Object} [cachedKeys] - an object to keep already calculated keys.
 * @returns {Array<Array<Layer>>}
 */
function groupByLayout(layers, cachedKeys) {
    const groups = {};
    for (let i = 0; i < layers.length; i++) {
        const k = cachedKeys && cachedKeys[layers[i].id] || getKey(layers[i]);
        // update the cache if there is one
        if (cachedKeys)
            cachedKeys[layers[i].id] = k;
        let group = groups[k];
        if (!group) {
            group = groups[k] = [];
        }
        group.push(layers[i]);
    }
    const result = [];
    for (const k in groups) {
        result.push(groups[k]);
    }
    return result;
}

//      
class StyleLayerIndex {
    constructor(layerConfigs) {
        this.keyCache = {};
        if (layerConfigs) {
            this.replace(layerConfigs);
        }
    }
    replace(layerConfigs) {
        this._layerConfigs = {};
        this._layers = {};
        this.update(layerConfigs, []);
    }
    update(layerConfigs, removedIds) {
        for (const layerConfig of layerConfigs) {
            this._layerConfigs[layerConfig.id] = layerConfig;
            const layer = this._layers[layerConfig.id] = ref_properties.createStyleLayer(layerConfig);
            layer.compileFilter();
            if (this.keyCache[layerConfig.id])
                delete this.keyCache[layerConfig.id];
        }
        for (const id of removedIds) {
            delete this.keyCache[id];
            delete this._layerConfigs[id];
            delete this._layers[id];
        }
        this.familiesBySource = {};
        const groups = groupByLayout(ref_properties.values(this._layerConfigs), this.keyCache);
        for (const layerConfigs of groups) {
            const layers = layerConfigs.map(layerConfig => this._layers[layerConfig.id]);
            const layer = layers[0];
            if (layer.visibility === 'none') {
                continue;
            }
            const sourceId = layer.source || '';
            let sourceGroup = this.familiesBySource[sourceId];
            if (!sourceGroup) {
                sourceGroup = this.familiesBySource[sourceId] = {};
            }
            const sourceLayerId = layer.sourceLayer || '_geojsonTileLayer';
            let sourceLayerFamilies = sourceGroup[sourceLayerId];
            if (!sourceLayerFamilies) {
                sourceLayerFamilies = sourceGroup[sourceLayerId] = [];
            }
            sourceLayerFamilies.push(layers);
        }
    }
}

//      
class RasterDEMTileWorkerSource {
    loadTile(params, callback) {
        const {uid, encoding, rawImageData, padding, buildQuadTree} = params;
        // Main thread will transfer ImageBitmap if offscreen decode with OffscreenCanvas is supported, else it will transfer an already decoded image.
        // Flow struggles to refine ImageBitmap type, likely due to the JSDom shim
        const imagePixels = ref_properties.window.ImageBitmap && rawImageData instanceof ref_properties.window.ImageBitmap ? this.getImageData(rawImageData, padding) : rawImageData;
        const dem = new ref_properties.DEMData(uid, imagePixels, encoding, padding < 1, buildQuadTree);
        callback(null, dem);
    }
    getImageData(imgBitmap, padding) {
        // Lazily initialize OffscreenCanvas
        if (!this.offscreenCanvas || !this.offscreenCanvasContext) {
            // Dem tiles are typically 256x256
            this.offscreenCanvas = new OffscreenCanvas(imgBitmap.width, imgBitmap.height);
            this.offscreenCanvasContext = this.offscreenCanvas.getContext('2d');
        }
        this.offscreenCanvas.width = imgBitmap.width;
        this.offscreenCanvas.height = imgBitmap.height;
        this.offscreenCanvasContext.drawImage(imgBitmap, 0, 0, imgBitmap.width, imgBitmap.height);
        // Insert or remove defined padding around the image to allow backfilling for neighboring data.
        const imgData = this.offscreenCanvasContext.getImageData(-padding, -padding, imgBitmap.width + 2 * padding, imgBitmap.height + 2 * padding);
        this.offscreenCanvasContext.clearRect(0, 0, this.offscreenCanvas.width, this.offscreenCanvas.height);
        return imgData;
    }
}

var geojsonRewind = rewind$1;
function rewind$1(gj, outer) {
    var type = gj && gj.type, i;
    if (type === 'FeatureCollection') {
        for (i = 0; i < gj.features.length; i++)
            rewind$1(gj.features[i], outer);
    } else if (type === 'GeometryCollection') {
        for (i = 0; i < gj.geometries.length; i++)
            rewind$1(gj.geometries[i], outer);
    } else if (type === 'Feature') {
        rewind$1(gj.geometry, outer);
    } else if (type === 'Polygon') {
        rewindRings(gj.coordinates, outer);
    } else if (type === 'MultiPolygon') {
        for (i = 0; i < gj.coordinates.length; i++)
            rewindRings(gj.coordinates[i], outer);
    }
    return gj;
}
function rewindRings(rings, outer) {
    if (rings.length === 0)
        return;
    rewindRing(rings[0], outer);
    for (var i = 1; i < rings.length; i++) {
        rewindRing(rings[i], !outer);
    }
}
function rewindRing(ring, dir) {
    var area = 0, err = 0;
    for (var i = 0, len = ring.length, j = len - 1; i < len; j = i++) {
        var k = (ring[i][0] - ring[j][0]) * (ring[j][1] + ring[i][1]);
        var m = area + k;
        err += Math.abs(area) >= Math.abs(k) ? area - m + k : k - m + area;
        area = m;
    }
    if (area + err >= 0 !== !!dir)
        ring.reverse();
}

//      
const toGeoJSON = ref_properties.vectorTile.VectorTileFeature.prototype.toGeoJSON;
// The feature type used by geojson-vt and supercluster. Should be extracted to
// global type and used in module definitions for those two modules.
class FeatureWrapper$1 {
    constructor(feature) {
        this._feature = feature;
        this.extent = ref_properties.EXTENT;
        this.type = feature.type;
        this.properties = feature.tags;
        // If the feature has a top-level `id` property, copy it over, but only
        // if it can be coerced to an integer, because this wrapper is used for
        // serializing geojson feature data into vector tile PBF data, and the
        // vector tile spec only supports integer values for feature ids --
        // allowing non-integer values here results in a non-compliant PBF
        // that causes an exception when it is parsed with vector-tile-js
        if ('id' in feature && !isNaN(feature.id)) {
            this.id = parseInt(feature.id, 10);
        }
    }
    loadGeometry() {
        if (this._feature.type === 1) {
            const geometry = [];
            for (const point of this._feature.geometry) {
                geometry.push([new ref_properties.pointGeometry(point[0], point[1])]);
            }
            return geometry;
        } else {
            const geometry = [];
            for (const ring of this._feature.geometry) {
                const newRing = [];
                for (const point of ring) {
                    newRing.push(new ref_properties.pointGeometry(point[0], point[1]));
                }
                geometry.push(newRing);
            }
            return geometry;
        }
    }
    toGeoJSON(x, y, z) {
        return toGeoJSON.call(this, x, y, z);
    }
}
class GeoJSONWrapper$2 {
    constructor(features) {
        this.layers = { '_geojsonTileLayer': this };
        this.name = '_geojsonTileLayer';
        this.extent = ref_properties.EXTENT;
        this.length = features.length;
        this._features = features;
    }
    feature(i) {
        return new FeatureWrapper$1(this._features[i]);
    }
}

var vtPbf = {exports: {}};

var Point = ref_properties.pointGeometry;
var VectorTileFeature = ref_properties.vectorTile.VectorTileFeature;
var geojson_wrapper = GeoJSONWrapper$1;
// conform to vectortile api
function GeoJSONWrapper$1(features, options) {
    this.options = options || {};
    this.features = features;
    this.length = features.length;
}
GeoJSONWrapper$1.prototype.feature = function (i) {
    return new FeatureWrapper(this.features[i], this.options.extent);
};
function FeatureWrapper(feature, extent) {
    this.id = typeof feature.id === 'number' ? feature.id : undefined;
    this.type = feature.type;
    this.rawGeometry = feature.type === 1 ? [feature.geometry] : feature.geometry;
    this.properties = feature.tags;
    this.extent = extent || 4096;
}
FeatureWrapper.prototype.loadGeometry = function () {
    var rings = this.rawGeometry;
    this.geometry = [];
    for (var i = 0; i < rings.length; i++) {
        var ring = rings[i];
        var newRing = [];
        for (var j = 0; j < ring.length; j++) {
            newRing.push(new Point(ring[j][0], ring[j][1]));
        }
        this.geometry.push(newRing);
    }
    return this.geometry;
};
FeatureWrapper.prototype.bbox = function () {
    if (!this.geometry)
        this.loadGeometry();
    var rings = this.geometry;
    var x1 = Infinity;
    var x2 = -Infinity;
    var y1 = Infinity;
    var y2 = -Infinity;
    for (var i = 0; i < rings.length; i++) {
        var ring = rings[i];
        for (var j = 0; j < ring.length; j++) {
            var coord = ring[j];
            x1 = Math.min(x1, coord.x);
            x2 = Math.max(x2, coord.x);
            y1 = Math.min(y1, coord.y);
            y2 = Math.max(y2, coord.y);
        }
    }
    return [
        x1,
        y1,
        x2,
        y2
    ];
};
FeatureWrapper.prototype.toGeoJSON = VectorTileFeature.prototype.toGeoJSON;

var Pbf = ref_properties.pbf;
var GeoJSONWrapper = geojson_wrapper;
vtPbf.exports = fromVectorTileJs;
vtPbf.exports.fromVectorTileJs = fromVectorTileJs;
vtPbf.exports.fromGeojsonVt = fromGeojsonVt;
vtPbf.exports.GeoJSONWrapper = GeoJSONWrapper;
/**
 * Serialize a vector-tile-js-created tile to pbf
 *
 * @param {Object} tile
 * @return {Buffer} uncompressed, pbf-serialized tile data
 */
function fromVectorTileJs(tile) {
    var out = new Pbf();
    writeTile(tile, out);
    return out.finish();
}
/**
 * Serialized a geojson-vt-created tile to pbf.
 *
 * @param {Object} layers - An object mapping layer names to geojson-vt-created vector tile objects
 * @param {Object} [options] - An object specifying the vector-tile specification version and extent that were used to create `layers`.
 * @param {Number} [options.version=1] - Version of vector-tile spec used
 * @param {Number} [options.extent=4096] - Extent of the vector tile
 * @return {Buffer} uncompressed, pbf-serialized tile data
 */
function fromGeojsonVt(layers, options) {
    options = options || {};
    var l = {};
    for (var k in layers) {
        l[k] = new GeoJSONWrapper(layers[k].features, options);
        l[k].name = k;
        l[k].version = options.version;
        l[k].extent = options.extent;
    }
    return fromVectorTileJs({ layers: l });
}
function writeTile(tile, pbf) {
    for (var key in tile.layers) {
        pbf.writeMessage(3, writeLayer, tile.layers[key]);
    }
}
function writeLayer(layer, pbf) {
    pbf.writeVarintField(15, layer.version || 1);
    pbf.writeStringField(1, layer.name || '');
    pbf.writeVarintField(5, layer.extent || 4096);
    var i;
    var context = {
        keys: [],
        values: [],
        keycache: {},
        valuecache: {}
    };
    for (i = 0; i < layer.length; i++) {
        context.feature = layer.feature(i);
        pbf.writeMessage(2, writeFeature, context);
    }
    var keys = context.keys;
    for (i = 0; i < keys.length; i++) {
        pbf.writeStringField(3, keys[i]);
    }
    var values = context.values;
    for (i = 0; i < values.length; i++) {
        pbf.writeMessage(4, writeValue, values[i]);
    }
}
function writeFeature(context, pbf) {
    var feature = context.feature;
    if (feature.id !== undefined) {
        pbf.writeVarintField(1, feature.id);
    }
    pbf.writeMessage(2, writeProperties, context);
    pbf.writeVarintField(3, feature.type);
    pbf.writeMessage(4, writeGeometry, feature);
}
function writeProperties(context, pbf) {
    var feature = context.feature;
    var keys = context.keys;
    var values = context.values;
    var keycache = context.keycache;
    var valuecache = context.valuecache;
    for (var key in feature.properties) {
        var value = feature.properties[key];
        var keyIndex = keycache[key];
        if (value === null)
            continue;
        // don't encode null value properties
        if (typeof keyIndex === 'undefined') {
            keys.push(key);
            keyIndex = keys.length - 1;
            keycache[key] = keyIndex;
        }
        pbf.writeVarint(keyIndex);
        var type = typeof value;
        if (type !== 'string' && type !== 'boolean' && type !== 'number') {
            value = JSON.stringify(value);
        }
        var valueKey = type + ':' + value;
        var valueIndex = valuecache[valueKey];
        if (typeof valueIndex === 'undefined') {
            values.push(value);
            valueIndex = values.length - 1;
            valuecache[valueKey] = valueIndex;
        }
        pbf.writeVarint(valueIndex);
    }
}
function command(cmd, length) {
    return (length << 3) + (cmd & 7);
}
function zigzag(num) {
    return num << 1 ^ num >> 31;
}
function writeGeometry(feature, pbf) {
    var geometry = feature.loadGeometry();
    var type = feature.type;
    var x = 0;
    var y = 0;
    var rings = geometry.length;
    for (var r = 0; r < rings; r++) {
        var ring = geometry[r];
        var count = 1;
        if (type === 1) {
            count = ring.length;
        }
        pbf.writeVarint(command(1, count));
        // moveto
        // do not write polygon closing path as lineto
        var lineCount = type === 3 ? ring.length - 1 : ring.length;
        for (var i = 0; i < lineCount; i++) {
            if (i === 1 && type !== 1) {
                pbf.writeVarint(command(2, lineCount - 1))    // lineto
;
            }
            var dx = ring[i].x - x;
            var dy = ring[i].y - y;
            pbf.writeVarint(zigzag(dx));
            pbf.writeVarint(zigzag(dy));
            x += dx;
            y += dy;
        }
        if (type === 3) {
            pbf.writeVarint(command(7, 1))    // closepath
;
        }
    }
}
function writeValue(value, pbf) {
    var type = typeof value;
    if (type === 'string') {
        pbf.writeStringField(1, value);
    } else if (type === 'boolean') {
        pbf.writeBooleanField(7, value);
    } else if (type === 'number') {
        if (value % 1 !== 0) {
            pbf.writeDoubleField(3, value);
        } else if (value < 0) {
            pbf.writeSVarintField(6, value);
        } else {
            pbf.writeVarintField(5, value);
        }
    }
}

function sortKD(ids, coords, nodeSize, left, right, depth) {
    if (right - left <= nodeSize)
        return;
    const m = left + right >> 1;
    select(ids, coords, m, left, right, depth % 2);
    sortKD(ids, coords, nodeSize, left, m - 1, depth + 1);
    sortKD(ids, coords, nodeSize, m + 1, right, depth + 1);
}
function select(ids, coords, k, left, right, inc) {
    while (right > left) {
        if (right - left > 600) {
            const n = right - left + 1;
            const m = k - left + 1;
            const z = Math.log(n);
            const s = 0.5 * Math.exp(2 * z / 3);
            const sd = 0.5 * Math.sqrt(z * s * (n - s) / n) * (m - n / 2 < 0 ? -1 : 1);
            const newLeft = Math.max(left, Math.floor(k - m * s / n + sd));
            const newRight = Math.min(right, Math.floor(k + (n - m) * s / n + sd));
            select(ids, coords, k, newLeft, newRight, inc);
        }
        const t = coords[2 * k + inc];
        let i = left;
        let j = right;
        swapItem(ids, coords, left, k);
        if (coords[2 * right + inc] > t)
            swapItem(ids, coords, left, right);
        while (i < j) {
            swapItem(ids, coords, i, j);
            i++;
            j--;
            while (coords[2 * i + inc] < t)
                i++;
            while (coords[2 * j + inc] > t)
                j--;
        }
        if (coords[2 * left + inc] === t)
            swapItem(ids, coords, left, j);
        else {
            j++;
            swapItem(ids, coords, j, right);
        }
        if (j <= k)
            left = j + 1;
        if (k <= j)
            right = j - 1;
    }
}
function swapItem(ids, coords, i, j) {
    swap(ids, i, j);
    swap(coords, 2 * i, 2 * j);
    swap(coords, 2 * i + 1, 2 * j + 1);
}
function swap(arr, i, j) {
    const tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
}

function range(ids, coords, minX, minY, maxX, maxY, nodeSize) {
    const stack = [
        0,
        ids.length - 1,
        0
    ];
    const result = [];
    let x, y;
    while (stack.length) {
        const axis = stack.pop();
        const right = stack.pop();
        const left = stack.pop();
        if (right - left <= nodeSize) {
            for (let i = left; i <= right; i++) {
                x = coords[2 * i];
                y = coords[2 * i + 1];
                if (x >= minX && x <= maxX && y >= minY && y <= maxY)
                    result.push(ids[i]);
            }
            continue;
        }
        const m = Math.floor((left + right) / 2);
        x = coords[2 * m];
        y = coords[2 * m + 1];
        if (x >= minX && x <= maxX && y >= minY && y <= maxY)
            result.push(ids[m]);
        const nextAxis = (axis + 1) % 2;
        if (axis === 0 ? minX <= x : minY <= y) {
            stack.push(left);
            stack.push(m - 1);
            stack.push(nextAxis);
        }
        if (axis === 0 ? maxX >= x : maxY >= y) {
            stack.push(m + 1);
            stack.push(right);
            stack.push(nextAxis);
        }
    }
    return result;
}

function within(ids, coords, qx, qy, r, nodeSize) {
    const stack = [
        0,
        ids.length - 1,
        0
    ];
    const result = [];
    const r2 = r * r;
    while (stack.length) {
        const axis = stack.pop();
        const right = stack.pop();
        const left = stack.pop();
        if (right - left <= nodeSize) {
            for (let i = left; i <= right; i++) {
                if (sqDist(coords[2 * i], coords[2 * i + 1], qx, qy) <= r2)
                    result.push(ids[i]);
            }
            continue;
        }
        const m = Math.floor((left + right) / 2);
        const x = coords[2 * m];
        const y = coords[2 * m + 1];
        if (sqDist(x, y, qx, qy) <= r2)
            result.push(ids[m]);
        const nextAxis = (axis + 1) % 2;
        if (axis === 0 ? qx - r <= x : qy - r <= y) {
            stack.push(left);
            stack.push(m - 1);
            stack.push(nextAxis);
        }
        if (axis === 0 ? qx + r >= x : qy + r >= y) {
            stack.push(m + 1);
            stack.push(right);
            stack.push(nextAxis);
        }
    }
    return result;
}
function sqDist(ax, ay, bx, by) {
    const dx = ax - bx;
    const dy = ay - by;
    return dx * dx + dy * dy;
}

const defaultGetX = p => p[0];
const defaultGetY = p => p[1];
class KDBush {
    constructor(points, getX = defaultGetX, getY = defaultGetY, nodeSize = 64, ArrayType = Float64Array) {
        this.nodeSize = nodeSize;
        this.points = points;
        const IndexArrayType = points.length < 65536 ? Uint16Array : Uint32Array;
        const ids = this.ids = new IndexArrayType(points.length);
        const coords = this.coords = new ArrayType(points.length * 2);
        for (let i = 0; i < points.length; i++) {
            ids[i] = i;
            coords[2 * i] = getX(points[i]);
            coords[2 * i + 1] = getY(points[i]);
        }
        sortKD(ids, coords, nodeSize, 0, ids.length - 1, 0);
    }
    range(minX, minY, maxX, maxY) {
        return range(this.ids, this.coords, minX, minY, maxX, maxY, this.nodeSize);
    }
    within(x, y, r) {
        return within(this.ids, this.coords, x, y, r, this.nodeSize);
    }
}

const defaultOptions = {
    minZoom: 0,
    // min zoom to generate clusters on
    maxZoom: 16,
    // max zoom level to cluster the points on
    minPoints: 2,
    // minimum points to form a cluster
    radius: 40,
    // cluster radius in pixels
    extent: 512,
    // tile extent (radius is calculated relative to it)
    nodeSize: 64,
    // size of the KD-tree leaf node, affects performance
    log: false,
    // whether to log timing info
    // whether to generate numeric ids for input features (in vector tiles)
    generateId: false,
    // a reduce function for calculating custom cluster properties
    reduce: null,
    // (accumulated, props) => { accumulated.sum += props.sum; }
    // properties to use for individual points when running the reducer
    map: props => props    // props => ({sum: props.my_value})
};
const fround = Math.fround || (tmp => x => {
    tmp[0] = +x;
    return tmp[0];
})(new Float32Array(1));
class Supercluster {
    constructor(options) {
        this.options = extend$1(Object.create(defaultOptions), options);
        this.trees = new Array(this.options.maxZoom + 1);
    }
    load(points) {
        const {log, minZoom, maxZoom, nodeSize} = this.options;
        if (log)
            console.time('total time');
        const timerId = `prepare ${ points.length } points`;
        if (log)
            console.time(timerId);
        this.points = points;
        // generate a cluster object for each point and index input points into a KD-tree
        let clusters = [];
        for (let i = 0; i < points.length; i++) {
            if (!points[i].geometry)
                continue;
            clusters.push(createPointCluster(points[i], i));
        }
        this.trees[maxZoom + 1] = new KDBush(clusters, getX, getY, nodeSize, Float32Array);
        if (log)
            console.timeEnd(timerId);
        // cluster points on max zoom, then cluster the results on previous zoom, etc.;
        // results in a cluster hierarchy across zoom levels
        for (let z = maxZoom; z >= minZoom; z--) {
            const now = +Date.now();
            // create a new set of clusters for the zoom and index them with a KD-tree
            clusters = this._cluster(clusters, z);
            this.trees[z] = new KDBush(clusters, getX, getY, nodeSize, Float32Array);
            if (log)
                console.log('z%d: %d clusters in %dms', z, clusters.length, +Date.now() - now);
        }
        if (log)
            console.timeEnd('total time');
        return this;
    }
    getClusters(bbox, zoom) {
        let minLng = ((bbox[0] + 180) % 360 + 360) % 360 - 180;
        const minLat = Math.max(-90, Math.min(90, bbox[1]));
        let maxLng = bbox[2] === 180 ? 180 : ((bbox[2] + 180) % 360 + 360) % 360 - 180;
        const maxLat = Math.max(-90, Math.min(90, bbox[3]));
        if (bbox[2] - bbox[0] >= 360) {
            minLng = -180;
            maxLng = 180;
        } else if (minLng > maxLng) {
            const easternHem = this.getClusters([
                minLng,
                minLat,
                180,
                maxLat
            ], zoom);
            const westernHem = this.getClusters([
                -180,
                minLat,
                maxLng,
                maxLat
            ], zoom);
            return easternHem.concat(westernHem);
        }
        const tree = this.trees[this._limitZoom(zoom)];
        const ids = tree.range(lngX(minLng), latY(maxLat), lngX(maxLng), latY(minLat));
        const clusters = [];
        for (const id of ids) {
            const c = tree.points[id];
            clusters.push(c.numPoints ? getClusterJSON(c) : this.points[c.index]);
        }
        return clusters;
    }
    getChildren(clusterId) {
        const originId = this._getOriginId(clusterId);
        const originZoom = this._getOriginZoom(clusterId);
        const errorMsg = 'No cluster with the specified id.';
        const index = this.trees[originZoom];
        if (!index)
            throw new Error(errorMsg);
        const origin = index.points[originId];
        if (!origin)
            throw new Error(errorMsg);
        const r = this.options.radius / (this.options.extent * Math.pow(2, originZoom - 1));
        const ids = index.within(origin.x, origin.y, r);
        const children = [];
        for (const id of ids) {
            const c = index.points[id];
            if (c.parentId === clusterId) {
                children.push(c.numPoints ? getClusterJSON(c) : this.points[c.index]);
            }
        }
        if (children.length === 0)
            throw new Error(errorMsg);
        return children;
    }
    getLeaves(clusterId, limit, offset) {
        limit = limit || 10;
        offset = offset || 0;
        const leaves = [];
        this._appendLeaves(leaves, clusterId, limit, offset, 0);
        return leaves;
    }
    getTile(z, x, y) {
        const tree = this.trees[this._limitZoom(z)];
        const z2 = Math.pow(2, z);
        const {extent, radius} = this.options;
        const p = radius / extent;
        const top = (y - p) / z2;
        const bottom = (y + 1 + p) / z2;
        const tile = { features: [] };
        this._addTileFeatures(tree.range((x - p) / z2, top, (x + 1 + p) / z2, bottom), tree.points, x, y, z2, tile);
        if (x === 0) {
            this._addTileFeatures(tree.range(1 - p / z2, top, 1, bottom), tree.points, z2, y, z2, tile);
        }
        if (x === z2 - 1) {
            this._addTileFeatures(tree.range(0, top, p / z2, bottom), tree.points, -1, y, z2, tile);
        }
        return tile.features.length ? tile : null;
    }
    getClusterExpansionZoom(clusterId) {
        let expansionZoom = this._getOriginZoom(clusterId) - 1;
        while (expansionZoom <= this.options.maxZoom) {
            const children = this.getChildren(clusterId);
            expansionZoom++;
            if (children.length !== 1)
                break;
            clusterId = children[0].properties.cluster_id;
        }
        return expansionZoom;
    }
    _appendLeaves(result, clusterId, limit, offset, skipped) {
        const children = this.getChildren(clusterId);
        for (const child of children) {
            const props = child.properties;
            if (props && props.cluster) {
                if (skipped + props.point_count <= offset) {
                    // skip the whole cluster
                    skipped += props.point_count;
                } else {
                    // enter the cluster
                    skipped = this._appendLeaves(result, props.cluster_id, limit, offset, skipped);    // exit the cluster
                }
            } else if (skipped < offset) {
                // skip a single point
                skipped++;
            } else {
                // add a single point
                result.push(child);
            }
            if (result.length === limit)
                break;
        }
        return skipped;
    }
    _addTileFeatures(ids, points, x, y, z2, tile) {
        for (const i of ids) {
            const c = points[i];
            const isCluster = c.numPoints;
            let tags, px, py;
            if (isCluster) {
                tags = getClusterProperties(c);
                px = c.x;
                py = c.y;
            } else {
                const p = this.points[c.index];
                tags = p.properties;
                px = lngX(p.geometry.coordinates[0]);
                py = latY(p.geometry.coordinates[1]);
            }
            const f = {
                type: 1,
                geometry: [[
                        Math.round(this.options.extent * (px * z2 - x)),
                        Math.round(this.options.extent * (py * z2 - y))
                    ]],
                tags
            };
            // assign id
            let id;
            if (isCluster) {
                id = c.id;
            } else if (this.options.generateId) {
                // optionally generate id
                id = c.index;
            } else if (this.points[c.index].id) {
                // keep id if already assigned
                id = this.points[c.index].id;
            }
            if (id !== undefined)
                f.id = id;
            tile.features.push(f);
        }
    }
    _limitZoom(z) {
        return Math.max(this.options.minZoom, Math.min(Math.floor(+z), this.options.maxZoom + 1));
    }
    _cluster(points, zoom) {
        const clusters = [];
        const {radius, extent, reduce, minPoints} = this.options;
        const r = radius / (extent * Math.pow(2, zoom));
        // loop through each point
        for (let i = 0; i < points.length; i++) {
            const p = points[i];
            // if we've already visited the point at this zoom level, skip it
            if (p.zoom <= zoom)
                continue;
            p.zoom = zoom;
            // find all nearby points
            const tree = this.trees[zoom + 1];
            const neighborIds = tree.within(p.x, p.y, r);
            const numPointsOrigin = p.numPoints || 1;
            let numPoints = numPointsOrigin;
            // count the number of points in a potential cluster
            for (const neighborId of neighborIds) {
                const b = tree.points[neighborId];
                // filter out neighbors that are already processed
                if (b.zoom > zoom)
                    numPoints += b.numPoints || 1;
            }
            // if there were neighbors to merge, and there are enough points to form a cluster
            if (numPoints > numPointsOrigin && numPoints >= minPoints) {
                let wx = p.x * numPointsOrigin;
                let wy = p.y * numPointsOrigin;
                let clusterProperties = reduce && numPointsOrigin > 1 ? this._map(p, true) : null;
                // encode both zoom and point index on which the cluster originated -- offset by total length of features
                const id = (i << 5) + (zoom + 1) + this.points.length;
                for (const neighborId of neighborIds) {
                    const b = tree.points[neighborId];
                    if (b.zoom <= zoom)
                        continue;
                    b.zoom = zoom;
                    // save the zoom (so it doesn't get processed twice)
                    const numPoints2 = b.numPoints || 1;
                    wx += b.x * numPoints2;
                    // accumulate coordinates for calculating weighted center
                    wy += b.y * numPoints2;
                    b.parentId = id;
                    if (reduce) {
                        if (!clusterProperties)
                            clusterProperties = this._map(p, true);
                        reduce(clusterProperties, this._map(b));
                    }
                }
                p.parentId = id;
                clusters.push(createCluster(wx / numPoints, wy / numPoints, id, numPoints, clusterProperties));
            } else {
                // left points as unclustered
                clusters.push(p);
                if (numPoints > 1) {
                    for (const neighborId of neighborIds) {
                        const b = tree.points[neighborId];
                        if (b.zoom <= zoom)
                            continue;
                        b.zoom = zoom;
                        clusters.push(b);
                    }
                }
            }
        }
        return clusters;
    }
    // get index of the point from which the cluster originated
    _getOriginId(clusterId) {
        return clusterId - this.points.length >> 5;
    }
    // get zoom of the point from which the cluster originated
    _getOriginZoom(clusterId) {
        return (clusterId - this.points.length) % 32;
    }
    _map(point, clone) {
        if (point.numPoints) {
            return clone ? extend$1({}, point.properties) : point.properties;
        }
        const original = this.points[point.index].properties;
        const result = this.options.map(original);
        return clone && result === original ? extend$1({}, result) : result;
    }
}
function createCluster(x, y, id, numPoints, properties) {
    return {
        x: fround(x),
        // weighted cluster center; round for consistency with Float32Array index
        y: fround(y),
        zoom: Infinity,
        // the last zoom the cluster was processed at
        id,
        // encodes index of the first child of the cluster and its zoom level
        parentId: -1,
        // parent cluster id
        numPoints,
        properties
    };
}
function createPointCluster(p, id) {
    const [x, y] = p.geometry.coordinates;
    return {
        x: fround(lngX(x)),
        // projected point coordinates
        y: fround(latY(y)),
        zoom: Infinity,
        // the last zoom the point was processed at
        index: id,
        // index of the source feature in the original input array,
        parentId: -1    // parent cluster id
    };
}
function getClusterJSON(cluster) {
    return {
        type: 'Feature',
        id: cluster.id,
        properties: getClusterProperties(cluster),
        geometry: {
            type: 'Point',
            coordinates: [
                xLng(cluster.x),
                yLat(cluster.y)
            ]
        }
    };
}
function getClusterProperties(cluster) {
    const count = cluster.numPoints;
    const abbrev = count >= 10000 ? `${ Math.round(count / 1000) }k` : count >= 1000 ? `${ Math.round(count / 100) / 10 }k` : count;
    return extend$1(extend$1({}, cluster.properties), {
        cluster: true,
        cluster_id: cluster.id,
        point_count: count,
        point_count_abbreviated: abbrev
    });
}
// longitude/latitude to spherical mercator in [0..1] range
function lngX(lng) {
    return lng / 360 + 0.5;
}
function latY(lat) {
    const sin = Math.sin(lat * Math.PI / 180);
    const y = 0.5 - 0.25 * Math.log((1 + sin) / (1 - sin)) / Math.PI;
    return y < 0 ? 0 : y > 1 ? 1 : y;
}
// spherical mercator to longitude/latitude
function xLng(x) {
    return (x - 0.5) * 360;
}
function yLat(y) {
    const y2 = (180 - y * 360) * Math.PI / 180;
    return 360 * Math.atan(Math.exp(y2)) / Math.PI - 90;
}
function extend$1(dest, src) {
    for (const id in src)
        dest[id] = src[id];
    return dest;
}
function getX(p) {
    return p.x;
}
function getY(p) {
    return p.y;
}

// calculate simplification data using optimized Douglas-Peucker algorithm
function simplify(coords, first, last, sqTolerance) {
    var maxSqDist = sqTolerance;
    var mid = last - first >> 1;
    var minPosToMid = last - first;
    var index;
    var ax = coords[first];
    var ay = coords[first + 1];
    var bx = coords[last];
    var by = coords[last + 1];
    for (var i = first + 3; i < last; i += 3) {
        var d = getSqSegDist(coords[i], coords[i + 1], ax, ay, bx, by);
        if (d > maxSqDist) {
            index = i;
            maxSqDist = d;
        } else if (d === maxSqDist) {
            // a workaround to ensure we choose a pivot close to the middle of the list,
            // reducing recursion depth, for certain degenerate inputs
            // https://github.com/mapbox/geojson-vt/issues/104
            var posToMid = Math.abs(i - mid);
            if (posToMid < minPosToMid) {
                index = i;
                minPosToMid = posToMid;
            }
        }
    }
    if (maxSqDist > sqTolerance) {
        if (index - first > 3)
            simplify(coords, first, index, sqTolerance);
        coords[index + 2] = maxSqDist;
        if (last - index > 3)
            simplify(coords, index, last, sqTolerance);
    }
}
// square distance from a point to a segment
function getSqSegDist(px, py, x, y, bx, by) {
    var dx = bx - x;
    var dy = by - y;
    if (dx !== 0 || dy !== 0) {
        var t = ((px - x) * dx + (py - y) * dy) / (dx * dx + dy * dy);
        if (t > 1) {
            x = bx;
            y = by;
        } else if (t > 0) {
            x += dx * t;
            y += dy * t;
        }
    }
    dx = px - x;
    dy = py - y;
    return dx * dx + dy * dy;
}

function createFeature(id, type, geom, tags) {
    var feature = {
        id: typeof id === 'undefined' ? null : id,
        type: type,
        geometry: geom,
        tags: tags,
        minX: Infinity,
        minY: Infinity,
        maxX: -Infinity,
        maxY: -Infinity
    };
    calcBBox(feature);
    return feature;
}
function calcBBox(feature) {
    var geom = feature.geometry;
    var type = feature.type;
    if (type === 'Point' || type === 'MultiPoint' || type === 'LineString') {
        calcLineBBox(feature, geom);
    } else if (type === 'Polygon' || type === 'MultiLineString') {
        for (var i = 0; i < geom.length; i++) {
            calcLineBBox(feature, geom[i]);
        }
    } else if (type === 'MultiPolygon') {
        for (i = 0; i < geom.length; i++) {
            for (var j = 0; j < geom[i].length; j++) {
                calcLineBBox(feature, geom[i][j]);
            }
        }
    }
}
function calcLineBBox(feature, geom) {
    for (var i = 0; i < geom.length; i += 3) {
        feature.minX = Math.min(feature.minX, geom[i]);
        feature.minY = Math.min(feature.minY, geom[i + 1]);
        feature.maxX = Math.max(feature.maxX, geom[i]);
        feature.maxY = Math.max(feature.maxY, geom[i + 1]);
    }
}

// converts GeoJSON feature into an intermediate projected JSON vector format with simplification data
function convert(data, options) {
    var features = [];
    if (data.type === 'FeatureCollection') {
        for (var i = 0; i < data.features.length; i++) {
            convertFeature(features, data.features[i], options, i);
        }
    } else if (data.type === 'Feature') {
        convertFeature(features, data, options);
    } else {
        // single geometry or a geometry collection
        convertFeature(features, { geometry: data }, options);
    }
    return features;
}
function convertFeature(features, geojson, options, index) {
    if (!geojson.geometry)
        return;
    var coords = geojson.geometry.coordinates;
    var type = geojson.geometry.type;
    var tolerance = Math.pow(options.tolerance / ((1 << options.maxZoom) * options.extent), 2);
    var geometry = [];
    var id = geojson.id;
    if (options.promoteId) {
        id = geojson.properties[options.promoteId];
    } else if (options.generateId) {
        id = index || 0;
    }
    if (type === 'Point') {
        convertPoint(coords, geometry);
    } else if (type === 'MultiPoint') {
        for (var i = 0; i < coords.length; i++) {
            convertPoint(coords[i], geometry);
        }
    } else if (type === 'LineString') {
        convertLine(coords, geometry, tolerance, false);
    } else if (type === 'MultiLineString') {
        if (options.lineMetrics) {
            // explode into linestrings to be able to track metrics
            for (i = 0; i < coords.length; i++) {
                geometry = [];
                convertLine(coords[i], geometry, tolerance, false);
                features.push(createFeature(id, 'LineString', geometry, geojson.properties));
            }
            return;
        } else {
            convertLines(coords, geometry, tolerance, false);
        }
    } else if (type === 'Polygon') {
        convertLines(coords, geometry, tolerance, true);
    } else if (type === 'MultiPolygon') {
        for (i = 0; i < coords.length; i++) {
            var polygon = [];
            convertLines(coords[i], polygon, tolerance, true);
            geometry.push(polygon);
        }
    } else if (type === 'GeometryCollection') {
        for (i = 0; i < geojson.geometry.geometries.length; i++) {
            convertFeature(features, {
                id: id,
                geometry: geojson.geometry.geometries[i],
                properties: geojson.properties
            }, options, index);
        }
        return;
    } else {
        throw new Error('Input data is not a valid GeoJSON object.');
    }
    features.push(createFeature(id, type, geometry, geojson.properties));
}
function convertPoint(coords, out) {
    out.push(projectX(coords[0]));
    out.push(projectY(coords[1]));
    out.push(0);
}
function convertLine(ring, out, tolerance, isPolygon) {
    var x0, y0;
    var size = 0;
    for (var j = 0; j < ring.length; j++) {
        var x = projectX(ring[j][0]);
        var y = projectY(ring[j][1]);
        out.push(x);
        out.push(y);
        out.push(0);
        if (j > 0) {
            if (isPolygon) {
                size += (x0 * y - x * y0) / 2;    // area
            } else {
                size += Math.sqrt(Math.pow(x - x0, 2) + Math.pow(y - y0, 2));    // length
            }
        }
        x0 = x;
        y0 = y;
    }
    var last = out.length - 3;
    out[2] = 1;
    simplify(out, 0, last, tolerance);
    out[last + 2] = 1;
    out.size = Math.abs(size);
    out.start = 0;
    out.end = out.size;
}
function convertLines(rings, out, tolerance, isPolygon) {
    for (var i = 0; i < rings.length; i++) {
        var geom = [];
        convertLine(rings[i], geom, tolerance, isPolygon);
        out.push(geom);
    }
}
function projectX(x) {
    return x / 360 + 0.5;
}
function projectY(y) {
    var sin = Math.sin(y * Math.PI / 180);
    var y2 = 0.5 - 0.25 * Math.log((1 + sin) / (1 - sin)) / Math.PI;
    return y2 < 0 ? 0 : y2 > 1 ? 1 : y2;
}

/* clip features between two axis-parallel lines:
 *     |        |
 *  ___|___     |     /
 * /   |   \____|____/
 *     |        |
 */
function clip(features, scale, k1, k2, axis, minAll, maxAll, options) {
    k1 /= scale;
    k2 /= scale;
    if (minAll >= k1 && maxAll < k2)
        return features;    // trivial accept
    else if (maxAll < k1 || minAll >= k2)
        return null;
    // trivial reject
    var clipped = [];
    for (var i = 0; i < features.length; i++) {
        var feature = features[i];
        var geometry = feature.geometry;
        var type = feature.type;
        var min = axis === 0 ? feature.minX : feature.minY;
        var max = axis === 0 ? feature.maxX : feature.maxY;
        if (min >= k1 && max < k2) {
            // trivial accept
            clipped.push(feature);
            continue;
        } else if (max < k1 || min >= k2) {
            // trivial reject
            continue;
        }
        var newGeometry = [];
        if (type === 'Point' || type === 'MultiPoint') {
            clipPoints(geometry, newGeometry, k1, k2, axis);
        } else if (type === 'LineString') {
            clipLine(geometry, newGeometry, k1, k2, axis, false, options.lineMetrics);
        } else if (type === 'MultiLineString') {
            clipLines(geometry, newGeometry, k1, k2, axis, false);
        } else if (type === 'Polygon') {
            clipLines(geometry, newGeometry, k1, k2, axis, true);
        } else if (type === 'MultiPolygon') {
            for (var j = 0; j < geometry.length; j++) {
                var polygon = [];
                clipLines(geometry[j], polygon, k1, k2, axis, true);
                if (polygon.length) {
                    newGeometry.push(polygon);
                }
            }
        }
        if (newGeometry.length) {
            if (options.lineMetrics && type === 'LineString') {
                for (j = 0; j < newGeometry.length; j++) {
                    clipped.push(createFeature(feature.id, type, newGeometry[j], feature.tags));
                }
                continue;
            }
            if (type === 'LineString' || type === 'MultiLineString') {
                if (newGeometry.length === 1) {
                    type = 'LineString';
                    newGeometry = newGeometry[0];
                } else {
                    type = 'MultiLineString';
                }
            }
            if (type === 'Point' || type === 'MultiPoint') {
                type = newGeometry.length === 3 ? 'Point' : 'MultiPoint';
            }
            clipped.push(createFeature(feature.id, type, newGeometry, feature.tags));
        }
    }
    return clipped.length ? clipped : null;
}
function clipPoints(geom, newGeom, k1, k2, axis) {
    for (var i = 0; i < geom.length; i += 3) {
        var a = geom[i + axis];
        if (a >= k1 && a <= k2) {
            newGeom.push(geom[i]);
            newGeom.push(geom[i + 1]);
            newGeom.push(geom[i + 2]);
        }
    }
}
function clipLine(geom, newGeom, k1, k2, axis, isPolygon, trackMetrics) {
    var slice = newSlice(geom);
    var intersect = axis === 0 ? intersectX : intersectY;
    var len = geom.start;
    var segLen, t;
    for (var i = 0; i < geom.length - 3; i += 3) {
        var ax = geom[i];
        var ay = geom[i + 1];
        var az = geom[i + 2];
        var bx = geom[i + 3];
        var by = geom[i + 4];
        var a = axis === 0 ? ax : ay;
        var b = axis === 0 ? bx : by;
        var exited = false;
        if (trackMetrics)
            segLen = Math.sqrt(Math.pow(ax - bx, 2) + Math.pow(ay - by, 2));
        if (a < k1) {
            // ---|-->  | (line enters the clip region from the left)
            if (b > k1) {
                t = intersect(slice, ax, ay, bx, by, k1);
                if (trackMetrics)
                    slice.start = len + segLen * t;
            }
        } else if (a > k2) {
            // |  <--|--- (line enters the clip region from the right)
            if (b < k2) {
                t = intersect(slice, ax, ay, bx, by, k2);
                if (trackMetrics)
                    slice.start = len + segLen * t;
            }
        } else {
            addPoint(slice, ax, ay, az);
        }
        if (b < k1 && a >= k1) {
            // <--|---  | or <--|-----|--- (line exits the clip region on the left)
            t = intersect(slice, ax, ay, bx, by, k1);
            exited = true;
        }
        if (b > k2 && a <= k2) {
            // |  ---|--> or ---|-----|--> (line exits the clip region on the right)
            t = intersect(slice, ax, ay, bx, by, k2);
            exited = true;
        }
        if (!isPolygon && exited) {
            if (trackMetrics)
                slice.end = len + segLen * t;
            newGeom.push(slice);
            slice = newSlice(geom);
        }
        if (trackMetrics)
            len += segLen;
    }
    // add the last point
    var last = geom.length - 3;
    ax = geom[last];
    ay = geom[last + 1];
    az = geom[last + 2];
    a = axis === 0 ? ax : ay;
    if (a >= k1 && a <= k2)
        addPoint(slice, ax, ay, az);
    // close the polygon if its endpoints are not the same after clipping
    last = slice.length - 3;
    if (isPolygon && last >= 3 && (slice[last] !== slice[0] || slice[last + 1] !== slice[1])) {
        addPoint(slice, slice[0], slice[1], slice[2]);
    }
    // add the final slice
    if (slice.length) {
        newGeom.push(slice);
    }
}
function newSlice(line) {
    var slice = [];
    slice.size = line.size;
    slice.start = line.start;
    slice.end = line.end;
    return slice;
}
function clipLines(geom, newGeom, k1, k2, axis, isPolygon) {
    for (var i = 0; i < geom.length; i++) {
        clipLine(geom[i], newGeom, k1, k2, axis, isPolygon, false);
    }
}
function addPoint(out, x, y, z) {
    out.push(x);
    out.push(y);
    out.push(z);
}
function intersectX(out, ax, ay, bx, by, x) {
    var t = (x - ax) / (bx - ax);
    out.push(x);
    out.push(ay + (by - ay) * t);
    out.push(1);
    return t;
}
function intersectY(out, ax, ay, bx, by, y) {
    var t = (y - ay) / (by - ay);
    out.push(ax + (bx - ax) * t);
    out.push(y);
    out.push(1);
    return t;
}

function wrap(features, options) {
    var buffer = options.buffer / options.extent;
    var merged = features;
    var left = clip(features, 1, -1 - buffer, buffer, 0, -1, 2, options);
    // left world copy
    var right = clip(features, 1, 1 - buffer, 2 + buffer, 0, -1, 2, options);
    // right world copy
    if (left || right) {
        merged = clip(features, 1, -buffer, 1 + buffer, 0, -1, 2, options) || [];
        // center world copy
        if (left)
            merged = shiftFeatureCoords(left, 1).concat(merged);
        // merge left into center
        if (right)
            merged = merged.concat(shiftFeatureCoords(right, -1));    // merge right into center
    }
    return merged;
}
function shiftFeatureCoords(features, offset) {
    var newFeatures = [];
    for (var i = 0; i < features.length; i++) {
        var feature = features[i], type = feature.type;
        var newGeometry;
        if (type === 'Point' || type === 'MultiPoint' || type === 'LineString') {
            newGeometry = shiftCoords(feature.geometry, offset);
        } else if (type === 'MultiLineString' || type === 'Polygon') {
            newGeometry = [];
            for (var j = 0; j < feature.geometry.length; j++) {
                newGeometry.push(shiftCoords(feature.geometry[j], offset));
            }
        } else if (type === 'MultiPolygon') {
            newGeometry = [];
            for (j = 0; j < feature.geometry.length; j++) {
                var newPolygon = [];
                for (var k = 0; k < feature.geometry[j].length; k++) {
                    newPolygon.push(shiftCoords(feature.geometry[j][k], offset));
                }
                newGeometry.push(newPolygon);
            }
        }
        newFeatures.push(createFeature(feature.id, type, newGeometry, feature.tags));
    }
    return newFeatures;
}
function shiftCoords(points, offset) {
    var newPoints = [];
    newPoints.size = points.size;
    if (points.start !== undefined) {
        newPoints.start = points.start;
        newPoints.end = points.end;
    }
    for (var i = 0; i < points.length; i += 3) {
        newPoints.push(points[i] + offset, points[i + 1], points[i + 2]);
    }
    return newPoints;
}

// Transforms the coordinates of each feature in the given tile from
// mercator-projected space into (extent x extent) tile space.
function transformTile(tile, extent) {
    if (tile.transformed)
        return tile;
    var z2 = 1 << tile.z, tx = tile.x, ty = tile.y, i, j, k;
    for (i = 0; i < tile.features.length; i++) {
        var feature = tile.features[i], geom = feature.geometry, type = feature.type;
        feature.geometry = [];
        if (type === 1) {
            for (j = 0; j < geom.length; j += 2) {
                feature.geometry.push(transformPoint(geom[j], geom[j + 1], extent, z2, tx, ty));
            }
        } else {
            for (j = 0; j < geom.length; j++) {
                var ring = [];
                for (k = 0; k < geom[j].length; k += 2) {
                    ring.push(transformPoint(geom[j][k], geom[j][k + 1], extent, z2, tx, ty));
                }
                feature.geometry.push(ring);
            }
        }
    }
    tile.transformed = true;
    return tile;
}
function transformPoint(x, y, extent, z2, tx, ty) {
    return [
        Math.round(extent * (x * z2 - tx)),
        Math.round(extent * (y * z2 - ty))
    ];
}

function createTile(features, z, tx, ty, options) {
    var tolerance = z === options.maxZoom ? 0 : options.tolerance / ((1 << z) * options.extent);
    var tile = {
        features: [],
        numPoints: 0,
        numSimplified: 0,
        numFeatures: 0,
        source: null,
        x: tx,
        y: ty,
        z: z,
        transformed: false,
        minX: 2,
        minY: 1,
        maxX: -1,
        maxY: 0
    };
    for (var i = 0; i < features.length; i++) {
        tile.numFeatures++;
        addFeature(tile, features[i], tolerance, options);
        var minX = features[i].minX;
        var minY = features[i].minY;
        var maxX = features[i].maxX;
        var maxY = features[i].maxY;
        if (minX < tile.minX)
            tile.minX = minX;
        if (minY < tile.minY)
            tile.minY = minY;
        if (maxX > tile.maxX)
            tile.maxX = maxX;
        if (maxY > tile.maxY)
            tile.maxY = maxY;
    }
    return tile;
}
function addFeature(tile, feature, tolerance, options) {
    var geom = feature.geometry, type = feature.type, simplified = [];
    if (type === 'Point' || type === 'MultiPoint') {
        for (var i = 0; i < geom.length; i += 3) {
            simplified.push(geom[i]);
            simplified.push(geom[i + 1]);
            tile.numPoints++;
            tile.numSimplified++;
        }
    } else if (type === 'LineString') {
        addLine(simplified, geom, tile, tolerance, false, false);
    } else if (type === 'MultiLineString' || type === 'Polygon') {
        for (i = 0; i < geom.length; i++) {
            addLine(simplified, geom[i], tile, tolerance, type === 'Polygon', i === 0);
        }
    } else if (type === 'MultiPolygon') {
        for (var k = 0; k < geom.length; k++) {
            var polygon = geom[k];
            for (i = 0; i < polygon.length; i++) {
                addLine(simplified, polygon[i], tile, tolerance, true, i === 0);
            }
        }
    }
    if (simplified.length) {
        var tags = feature.tags || null;
        if (type === 'LineString' && options.lineMetrics) {
            tags = {};
            for (var key in feature.tags)
                tags[key] = feature.tags[key];
            tags['mapbox_clip_start'] = geom.start / geom.size;
            tags['mapbox_clip_end'] = geom.end / geom.size;
        }
        var tileFeature = {
            geometry: simplified,
            type: type === 'Polygon' || type === 'MultiPolygon' ? 3 : type === 'LineString' || type === 'MultiLineString' ? 2 : 1,
            tags: tags
        };
        if (feature.id !== null) {
            tileFeature.id = feature.id;
        }
        tile.features.push(tileFeature);
    }
}
function addLine(result, geom, tile, tolerance, isPolygon, isOuter) {
    var sqTolerance = tolerance * tolerance;
    if (tolerance > 0 && geom.size < (isPolygon ? sqTolerance : tolerance)) {
        tile.numPoints += geom.length / 3;
        return;
    }
    var ring = [];
    for (var i = 0; i < geom.length; i += 3) {
        if (tolerance === 0 || geom[i + 2] > sqTolerance) {
            tile.numSimplified++;
            ring.push(geom[i]);
            ring.push(geom[i + 1]);
        }
        tile.numPoints++;
    }
    if (isPolygon)
        rewind(ring, isOuter);
    result.push(ring);
}
function rewind(ring, clockwise) {
    var area = 0;
    for (var i = 0, len = ring.length, j = len - 2; i < len; j = i, i += 2) {
        area += (ring[i] - ring[j]) * (ring[i + 1] + ring[j + 1]);
    }
    if (area > 0 === clockwise) {
        for (i = 0, len = ring.length; i < len / 2; i += 2) {
            var x = ring[i];
            var y = ring[i + 1];
            ring[i] = ring[len - 2 - i];
            ring[i + 1] = ring[len - 1 - i];
            ring[len - 2 - i] = x;
            ring[len - 1 - i] = y;
        }
    }
}

// final simplified tile generation
function geojsonvt(data, options) {
    return new GeoJSONVT(data, options);
}
function GeoJSONVT(data, options) {
    options = this.options = extend(Object.create(this.options), options);
    var debug = options.debug;
    if (debug)
        console.time('preprocess data');
    if (options.maxZoom < 0 || options.maxZoom > 24)
        throw new Error('maxZoom should be in the 0-24 range');
    if (options.promoteId && options.generateId)
        throw new Error('promoteId and generateId cannot be used together.');
    var features = convert(data, options);
    this.tiles = {};
    this.tileCoords = [];
    if (debug) {
        console.timeEnd('preprocess data');
        console.log('index: maxZoom: %d, maxPoints: %d', options.indexMaxZoom, options.indexMaxPoints);
        console.time('generate tiles');
        this.stats = {};
        this.total = 0;
    }
    features = wrap(features, options);
    // start slicing from the top tile down
    if (features.length)
        this.splitTile(features, 0, 0, 0);
    if (debug) {
        if (features.length)
            console.log('features: %d, points: %d', this.tiles[0].numFeatures, this.tiles[0].numPoints);
        console.timeEnd('generate tiles');
        console.log('tiles generated:', this.total, JSON.stringify(this.stats));
    }
}
GeoJSONVT.prototype.options = {
    maxZoom: 14,
    // max zoom to preserve detail on
    indexMaxZoom: 5,
    // max zoom in the tile index
    indexMaxPoints: 100000,
    // max number of points per tile in the tile index
    tolerance: 3,
    // simplification tolerance (higher means simpler)
    extent: 4096,
    // tile extent
    buffer: 64,
    // tile buffer on each side
    lineMetrics: false,
    // whether to calculate line metrics
    promoteId: null,
    // name of a feature property to be promoted to feature.id
    generateId: false,
    // whether to generate feature ids. Cannot be used with promoteId
    debug: 0    // logging level (0, 1 or 2)
};
GeoJSONVT.prototype.splitTile = function (features, z, x, y, cz, cx, cy) {
    var stack = [
            features,
            z,
            x,
            y
        ], options = this.options, debug = options.debug;
    // avoid recursion by using a processing queue
    while (stack.length) {
        y = stack.pop();
        x = stack.pop();
        z = stack.pop();
        features = stack.pop();
        var z2 = 1 << z, id = toID(z, x, y), tile = this.tiles[id];
        if (!tile) {
            if (debug > 1)
                console.time('creation');
            tile = this.tiles[id] = createTile(features, z, x, y, options);
            this.tileCoords.push({
                z: z,
                x: x,
                y: y
            });
            if (debug) {
                if (debug > 1) {
                    console.log('tile z%d-%d-%d (features: %d, points: %d, simplified: %d)', z, x, y, tile.numFeatures, tile.numPoints, tile.numSimplified);
                    console.timeEnd('creation');
                }
                var key = 'z' + z;
                this.stats[key] = (this.stats[key] || 0) + 1;
                this.total++;
            }
        }
        // save reference to original geometry in tile so that we can drill down later if we stop now
        tile.source = features;
        // if it's the first-pass tiling
        if (!cz) {
            // stop tiling if we reached max zoom, or if the tile is too simple
            if (z === options.indexMaxZoom || tile.numPoints <= options.indexMaxPoints)
                continue;    // if a drilldown to a specific tile
        } else {
            // stop tiling if we reached base zoom or our target tile zoom
            if (z === options.maxZoom || z === cz)
                continue;
            // stop tiling if it's not an ancestor of the target tile
            var m = 1 << cz - z;
            if (x !== Math.floor(cx / m) || y !== Math.floor(cy / m))
                continue;
        }
        // if we slice further down, no need to keep source geometry
        tile.source = null;
        if (features.length === 0)
            continue;
        if (debug > 1)
            console.time('clipping');
        // values we'll use for clipping
        var k1 = 0.5 * options.buffer / options.extent, k2 = 0.5 - k1, k3 = 0.5 + k1, k4 = 1 + k1, tl, bl, tr, br, left, right;
        tl = bl = tr = br = null;
        left = clip(features, z2, x - k1, x + k3, 0, tile.minX, tile.maxX, options);
        right = clip(features, z2, x + k2, x + k4, 0, tile.minX, tile.maxX, options);
        features = null;
        if (left) {
            tl = clip(left, z2, y - k1, y + k3, 1, tile.minY, tile.maxY, options);
            bl = clip(left, z2, y + k2, y + k4, 1, tile.minY, tile.maxY, options);
            left = null;
        }
        if (right) {
            tr = clip(right, z2, y - k1, y + k3, 1, tile.minY, tile.maxY, options);
            br = clip(right, z2, y + k2, y + k4, 1, tile.minY, tile.maxY, options);
            right = null;
        }
        if (debug > 1)
            console.timeEnd('clipping');
        stack.push(tl || [], z + 1, x * 2, y * 2);
        stack.push(bl || [], z + 1, x * 2, y * 2 + 1);
        stack.push(tr || [], z + 1, x * 2 + 1, y * 2);
        stack.push(br || [], z + 1, x * 2 + 1, y * 2 + 1);
    }
};
GeoJSONVT.prototype.getTile = function (z, x, y) {
    var options = this.options, extent = options.extent, debug = options.debug;
    if (z < 0 || z > 24)
        return null;
    var z2 = 1 << z;
    x = (x % z2 + z2) % z2;
    // wrap tile x coordinate
    var id = toID(z, x, y);
    if (this.tiles[id])
        return transformTile(this.tiles[id], extent);
    if (debug > 1)
        console.log('drilling down to z%d-%d-%d', z, x, y);
    var z0 = z, x0 = x, y0 = y, parent;
    while (!parent && z0 > 0) {
        z0--;
        x0 = Math.floor(x0 / 2);
        y0 = Math.floor(y0 / 2);
        parent = this.tiles[toID(z0, x0, y0)];
    }
    if (!parent || !parent.source)
        return null;
    // if we found a parent tile containing the original geometry, we can drill down from it
    if (debug > 1)
        console.log('found parent tile z%d-%d-%d', z0, x0, y0);
    if (debug > 1)
        console.time('drilling down');
    this.splitTile(parent.source, z0, x0, y0, z, x, y);
    if (debug > 1)
        console.timeEnd('drilling down');
    return this.tiles[id] ? transformTile(this.tiles[id], extent) : null;
};
function toID(z, x, y) {
    return ((1 << z) * y + x) * 32 + z;
}
function extend(dest, src) {
    for (var i in src)
        dest[i] = src[i];
    return dest;
}

//      
function loadGeoJSONTile(params, callback) {
    const canonical = params.tileID.canonical;
    if (!this._geoJSONIndex) {
        return callback(null, null);    // we couldn't load the file
    }
    const geoJSONTile = this._geoJSONIndex.getTile(canonical.z, canonical.x, canonical.y);
    if (!geoJSONTile) {
        return callback(null, null);    // nothing in the given tile
    }
    const geojsonWrapper = new GeoJSONWrapper$2(geoJSONTile.features);
    // Encode the geojson-vt tile into binary vector tile form.  This
    // is a convenience that allows `FeatureIndex` to operate the same way
    // across `VectorTileSource` and `GeoJSONSource` data.
    let pbf = vtPbf.exports(geojsonWrapper);
    if (pbf.byteOffset !== 0 || pbf.byteLength !== pbf.buffer.byteLength) {
        // Compatibility with node Buffer (https://github.com/mapbox/pbf/issues/35)
        pbf = new Uint8Array(pbf);
    }
    callback(null, {
        vectorTile: geojsonWrapper,
        rawData: pbf.buffer
    });
}
/**
 * The {@link WorkerSource} implementation that supports {@link GeoJSONSource}.
 * This class is designed to be easily reused to support custom source types
 * for data formats that can be parsed/converted into an in-memory GeoJSON
 * representation.  To do so, create it with
 * `new GeoJSONWorkerSource(actor, layerIndex, customLoadGeoJSONFunction)`.
 * For a full example, see [mapbox-gl-topojson](https://github.com/developmentseed/mapbox-gl-topojson).
 *
 * @private
 */
class GeoJSONWorkerSource extends ref_properties.VectorTileWorkerSource {
    /**
     * @param [loadGeoJSON] Optional method for custom loading/parsing of
     * GeoJSON based on parameters passed from the main-thread Source.
     * See {@link GeoJSONWorkerSource#loadGeoJSON}.
     * @private
     */
    constructor(actor, layerIndex, availableImages, isSpriteLoaded, loadGeoJSON) {
        super(actor, layerIndex, availableImages, isSpriteLoaded, loadGeoJSONTile);
        if (loadGeoJSON) {
            this.loadGeoJSON = loadGeoJSON;
        }
    }
    /**
     * Fetches (if appropriate), parses, and index geojson data into tiles. This
     * preparatory method must be called before {@link GeoJSONWorkerSource#loadTile}
     * can correctly serve up tiles.
     *
     * Defers to {@link GeoJSONWorkerSource#loadGeoJSON} for the fetching/parsing,
     * expecting `callback(error, data)` to be called with either an error or a
     * parsed GeoJSON object.
     *
     * When `loadData` requests come in faster than they can be processed,
     * they are coalesced into a single request using the latest data.
     * See {@link GeoJSONWorkerSource#coalesce}
     *
     * @param params
     * @param callback
     * @private
     */
    loadData(params, callback) {
        const requestParam = params && params.request;
        const perf = requestParam && requestParam.collectResourceTiming;
        this.loadGeoJSON(params, (err, data) => {
            if (err || !data) {
                return callback(err);
            } else if (typeof data !== 'object') {
                return callback(new Error(`Input data given to '${ params.source }' is not a valid GeoJSON object.`));
            } else {
                geojsonRewind(data, true);
                try {
                    if (params.filter) {
                        const compiled = ref_properties.createExpression(params.filter, {
                            type: 'boolean',
                            'property-type': 'data-driven',
                            overridable: false,
                            transition: false
                        });
                        if (compiled.result === 'error')
                            throw new Error(compiled.value.map(err => `${ err.key }: ${ err.message }`).join(', '));
                        const features = data.features.filter(feature => compiled.value.evaluate({ zoom: 0 }, feature));
                        data = {
                            type: 'FeatureCollection',
                            features
                        };
                    }
                    this._geoJSONIndex = params.cluster ? new Supercluster(getSuperclusterOptions(params)).load(data.features) : geojsonvt(data, params.geojsonVtOptions);
                } catch (err) {
                    return callback(err);
                }
                this.loaded = {};
                const result = {};
                if (perf) {
                    const resourceTimingData = ref_properties.getPerformanceMeasurement(requestParam);
                    // it's necessary to eval the result of getEntriesByName() here via parse/stringify
                    // late evaluation in the main thread causes TypeError: illegal invocation
                    if (resourceTimingData) {
                        result.resourceTiming = {};
                        result.resourceTiming[params.source] = JSON.parse(JSON.stringify(resourceTimingData));
                    }
                }
                callback(null, result);
            }
        });
    }
    /**
    * Implements {@link WorkerSource#reloadTile}.
    *
    * If the tile is loaded, uses the implementation in VectorTileWorkerSource.
    * Otherwise, such as after a setData() call, we load the tile fresh.
    *
    * @param params
    * @param params.uid The UID for this tile.
    * @private
    */
    reloadTile(params, callback) {
        const loaded = this.loaded, uid = params.uid;
        if (loaded && loaded[uid]) {
            return super.reloadTile(params, callback);
        } else {
            return this.loadTile(params, callback);
        }
    }
    /**
     * Fetch and parse GeoJSON according to the given params.  Calls `callback`
     * with `(err, data)`, where `data` is a parsed GeoJSON object.
     *
     * GeoJSON is loaded and parsed from `params.url` if it exists, or else
     * expected as a literal (string or object) `params.data`.
     *
     * @param params
     * @param [params.url] A URL to the remote GeoJSON data.
     * @param [params.data] Literal GeoJSON data. Must be provided if `params.url` is not.
     * @private
     */
    loadGeoJSON(params, callback) {
        // Because of same origin issues, urls must either include an explicit
        // origin or absolute path.
        // ie: /foo/bar.json or http://example.com/bar.json
        // but not ../foo/bar.json
        if (params.request) {
            ref_properties.getJSON(params.request, callback);
        } else if (typeof params.data === 'string') {
            try {
                return callback(null, JSON.parse(params.data));
            } catch (e) {
                return callback(new Error(`Input data given to '${ params.source }' is not a valid GeoJSON object.`));
            }
        } else {
            return callback(new Error(`Input data given to '${ params.source }' is not a valid GeoJSON object.`));
        }
    }
    getClusterExpansionZoom(params, callback) {
        try {
            callback(null, this._geoJSONIndex.getClusterExpansionZoom(params.clusterId));
        } catch (e) {
            callback(e);
        }
    }
    getClusterChildren(params, callback) {
        try {
            callback(null, this._geoJSONIndex.getChildren(params.clusterId));
        } catch (e) {
            callback(e);
        }
    }
    getClusterLeaves(params, callback) {
        try {
            callback(null, this._geoJSONIndex.getLeaves(params.clusterId, params.limit, params.offset));
        } catch (e) {
            callback(e);
        }
    }
}
function getSuperclusterOptions({superclusterOptions, clusterProperties}) {
    if (!clusterProperties || !superclusterOptions)
        return superclusterOptions;
    const mapExpressions = {};
    const reduceExpressions = {};
    const globals = {
        accumulated: null,
        zoom: 0
    };
    const feature = { properties: null };
    const propertyNames = Object.keys(clusterProperties);
    for (const key of propertyNames) {
        const [operator, mapExpression] = clusterProperties[key];
        const mapExpressionParsed = ref_properties.createExpression(mapExpression);
        const reduceExpressionParsed = ref_properties.createExpression(typeof operator === 'string' ? [
            operator,
            ['accumulated'],
            [
                'get',
                key
            ]
        ] : operator);
        mapExpressions[key] = mapExpressionParsed.value;
        reduceExpressions[key] = reduceExpressionParsed.value;
    }
    superclusterOptions.map = pointProperties => {
        feature.properties = pointProperties;
        const properties = {};
        for (const key of propertyNames) {
            properties[key] = mapExpressions[key].evaluate(globals, feature);
        }
        return properties;
    };
    superclusterOptions.reduce = (accumulated, clusterProperties) => {
        feature.properties = clusterProperties;
        for (const key of propertyNames) {
            globals.accumulated = accumulated[key];
            accumulated[key] = reduceExpressions[key].evaluate(globals, feature);
        }
    };
    return superclusterOptions;
}

//      
/**
 * @private
 */
class Worker {
    constructor(self) {
        this.self = self;
        this.actor = new ref_properties.Actor(self, this);
        this.layerIndexes = {};
        this.availableImages = {};
        this.isSpriteLoaded = {};
        this.projections = {};
        this.defaultProjection = ref_properties.getProjection({ name: 'mercator' });
        this.workerSourceTypes = {
            vector: ref_properties.VectorTileWorkerSource,
            geojson: GeoJSONWorkerSource
        };
        // [mapId][sourceType][sourceName] => worker source instance
        this.workerSources = {};
        this.demWorkerSources = {};
        this.self.registerWorkerSource = (name, WorkerSource) => {
            if (this.workerSourceTypes[name]) {
                throw new Error(`Worker source with name "${ name }" already registered.`);
            }
            this.workerSourceTypes[name] = WorkerSource;
        };
        // This is invoked by the RTL text plugin when the download via the `importScripts` call has finished, and the code has been parsed.
        this.self.registerRTLTextPlugin = rtlTextPlugin => {
            if (ref_properties.plugin.isParsed()) {
                throw new Error('RTL text plugin already registered.');
            }
            ref_properties.plugin['applyArabicShaping'] = rtlTextPlugin.applyArabicShaping;
            ref_properties.plugin['processBidirectionalText'] = rtlTextPlugin.processBidirectionalText;
            ref_properties.plugin['processStyledBidirectionalText'] = rtlTextPlugin.processStyledBidirectionalText;
        };
    }
    clearCaches(mapId, unused, callback) {
        delete this.layerIndexes[mapId];
        delete this.availableImages[mapId];
        delete this.workerSources[mapId];
        delete this.demWorkerSources[mapId];
        callback();
    }
    checkIfReady(mapID, unused, callback) {
        // noop, used to check if a worker is fully set up and ready to receive messages
        callback();
    }
    setReferrer(mapID, referrer) {
        this.referrer = referrer;
    }
    spriteLoaded(mapId, bool) {
        this.isSpriteLoaded[mapId] = bool;
        for (const workerSource in this.workerSources[mapId]) {
            const ws = this.workerSources[mapId][workerSource];
            for (const source in ws) {
                if (ws[source] instanceof ref_properties.VectorTileWorkerSource) {
                    ws[source].isSpriteLoaded = bool;
                    ws[source].fire(new ref_properties.Event('isSpriteLoaded'));
                }
            }
        }
    }
    setImages(mapId, images, callback) {
        this.availableImages[mapId] = images;
        for (const workerSource in this.workerSources[mapId]) {
            const ws = this.workerSources[mapId][workerSource];
            for (const source in ws) {
                ws[source].availableImages = images;
            }
        }
        callback();
    }
    enableTerrain(mapId, enable, callback) {
        this.terrain = enable;
        callback();
    }
    setProjection(mapId, config) {
        this.projections[mapId] = ref_properties.getProjection(config);
    }
    setLayers(mapId, layers, callback) {
        this.getLayerIndex(mapId).replace(layers);
        callback();
    }
    updateLayers(mapId, params, callback) {
        this.getLayerIndex(mapId).update(params.layers, params.removedIds);
        callback();
    }
    loadTile(mapId, params, callback) {
        const p = this.enableTerrain ? ref_properties.extend({ enableTerrain: this.terrain }, params) : params;
        p.projection = this.projections[mapId] || this.defaultProjection;
        this.getWorkerSource(mapId, params.type, params.source).loadTile(p, callback);
    }
    loadDEMTile(mapId, params, callback) {
        const p = this.enableTerrain ? ref_properties.extend({ buildQuadTree: this.terrain }, params) : params;
        this.getDEMWorkerSource(mapId, params.source).loadTile(p, callback);
    }
    reloadTile(mapId, params, callback) {
        const p = this.enableTerrain ? ref_properties.extend({ enableTerrain: this.terrain }, params) : params;
        p.projection = this.projections[mapId] || this.defaultProjection;
        this.getWorkerSource(mapId, params.type, params.source).reloadTile(p, callback);
    }
    abortTile(mapId, params, callback) {
        this.getWorkerSource(mapId, params.type, params.source).abortTile(params, callback);
    }
    removeTile(mapId, params, callback) {
        this.getWorkerSource(mapId, params.type, params.source).removeTile(params, callback);
    }
    removeSource(mapId, params, callback) {
        if (!this.workerSources[mapId] || !this.workerSources[mapId][params.type] || !this.workerSources[mapId][params.type][params.source]) {
            return;
        }
        const worker = this.workerSources[mapId][params.type][params.source];
        delete this.workerSources[mapId][params.type][params.source];
        if (worker.removeSource !== undefined) {
            worker.removeSource(params, callback);
        } else {
            callback();
        }
    }
    /**
     * Load a {@link WorkerSource} script at params.url.  The script is run
     * (using importScripts) with `registerWorkerSource` in scope, which is a
     * function taking `(name, workerSourceObject)`.
     *  @private
     */
    loadWorkerSource(map, params, callback) {
        try {
            this.self.importScripts(params.url);
            callback();
        } catch (e) {
            callback(e.toString());
        }
    }
    syncRTLPluginState(map, state, callback) {
        try {
            ref_properties.plugin.setState(state);
            const pluginURL = ref_properties.plugin.getPluginURL();
            if (ref_properties.plugin.isLoaded() && !ref_properties.plugin.isParsed() && pluginURL != null    // Not possible when `isLoaded` is true, but keeps flow happy
) {
                this.self.importScripts(pluginURL);
                const complete = ref_properties.plugin.isParsed();
                const error = complete ? undefined : new Error(`RTL Text Plugin failed to import scripts from ${ pluginURL }`);
                callback(error, complete);
            }
        } catch (e) {
            callback(e.toString());
        }
    }
    getAvailableImages(mapId) {
        let availableImages = this.availableImages[mapId];
        if (!availableImages) {
            availableImages = [];
        }
        return availableImages;
    }
    getLayerIndex(mapId) {
        let layerIndexes = this.layerIndexes[mapId];
        if (!layerIndexes) {
            layerIndexes = this.layerIndexes[mapId] = new StyleLayerIndex();
        }
        return layerIndexes;
    }
    getWorkerSource(mapId, type, source) {
        if (!this.workerSources[mapId])
            this.workerSources[mapId] = {};
        if (!this.workerSources[mapId][type])
            this.workerSources[mapId][type] = {};
        if (!this.workerSources[mapId][type][source]) {
            // use a wrapped actor so that we can attach a target mapId param
            // to any messages invoked by the WorkerSource
            const actor = {
                send: (type, data, callback, _, mustQueue, metadata) => {
                    this.actor.send(type, data, callback, mapId, mustQueue, metadata);
                },
                scheduler: this.actor.scheduler
            };
            this.workerSources[mapId][type][source] = new this.workerSourceTypes[type](actor, this.getLayerIndex(mapId), this.getAvailableImages(mapId), this.isSpriteLoaded[mapId]);
        }
        return this.workerSources[mapId][type][source];
    }
    getDEMWorkerSource(mapId, source) {
        if (!this.demWorkerSources[mapId])
            this.demWorkerSources[mapId] = {};
        if (!this.demWorkerSources[mapId][source]) {
            this.demWorkerSources[mapId][source] = new RasterDEMTileWorkerSource();
        }
        return this.demWorkerSources[mapId][source];
    }
    enforceCacheSizeLimit(mapId, limit) {
        ref_properties.enforceCacheSizeLimit(limit);
    }
    getWorkerPerformanceMetrics(mapId, params, callback) {
        callback(undefined, void 0);
    }
}
/* global self, WorkerGlobalScope */
if (typeof WorkerGlobalScope !== 'undefined' && typeof self !== 'undefined' && self instanceof WorkerGlobalScope) {
    self.worker = new Worker(self);
}

return Worker;

}));
//# sourceMappingURL=data:application/json;charset=utf-8;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoid29ya2VyLmpzIiwic291cmNlcyI6WyIuLi8uLi8uLi9zcmMvc3R5bGUtc3BlYy9ncm91cF9ieV9sYXlvdXQuanMiLCIuLi8uLi8uLi9zcmMvc3R5bGUvc3R5bGVfbGF5ZXJfaW5kZXguanMiLCIuLi8uLi8uLi9zcmMvc291cmNlL3Jhc3Rlcl9kZW1fdGlsZV93b3JrZXJfc291cmNlLmpzIiwiLi4vLi4vLi4vbm9kZV9tb2R1bGVzL0BtYXBib3gvZ2VvanNvbi1yZXdpbmQvaW5kZXguanMiLCIuLi8uLi8uLi9zcmMvc291cmNlL2dlb2pzb25fd3JhcHBlci5qcyIsIi4uLy4uLy4uL25vZGVfbW9kdWxlcy92dC1wYmYvbGliL2dlb2pzb25fd3JhcHBlci5qcyIsIi4uLy4uLy4uL25vZGVfbW9kdWxlcy92dC1wYmYvaW5kZXguanMiLCIuLi8uLi8uLi9ub2RlX21vZHVsZXMva2RidXNoL3NyYy9zb3J0LmpzIiwiLi4vLi4vLi4vbm9kZV9tb2R1bGVzL2tkYnVzaC9zcmMvcmFuZ2UuanMiLCIuLi8uLi8uLi9ub2RlX21vZHVsZXMva2RidXNoL3NyYy93aXRoaW4uanMiLCIuLi8uLi8uLi9ub2RlX21vZHVsZXMva2RidXNoL3NyYy9pbmRleC5qcyIsIi4uLy4uLy4uL25vZGVfbW9kdWxlcy9zdXBlcmNsdXN0ZXIvaW5kZXguanMiLCIuLi8uLi8uLi9ub2RlX21vZHVsZXMvZ2VvanNvbi12dC9zcmMvc2ltcGxpZnkuanMiLCIuLi8uLi8uLi9ub2RlX21vZHVsZXMvZ2VvanNvbi12dC9zcmMvZmVhdHVyZS5qcyIsIi4uLy4uLy4uL25vZGVfbW9kdWxlcy9nZW9qc29uLXZ0L3NyYy9jb252ZXJ0LmpzIiwiLi4vLi4vLi4vbm9kZV9tb2R1bGVzL2dlb2pzb24tdnQvc3JjL2NsaXAuanMiLCIuLi8uLi8uLi9ub2RlX21vZHVsZXMvZ2VvanNvbi12dC9zcmMvd3JhcC5qcyIsIi4uLy4uLy4uL25vZGVfbW9kdWxlcy9nZW9qc29uLXZ0L3NyYy90cmFuc2Zvcm0uanMiLCIuLi8uLi8uLi9ub2RlX21vZHVsZXMvZ2VvanNvbi12dC9zcmMvdGlsZS5qcyIsIi4uLy4uLy4uL25vZGVfbW9kdWxlcy9nZW9qc29uLXZ0L3NyYy9pbmRleC5qcyIsIi4uLy4uLy4uL3NyYy9zb3VyY2UvZ2VvanNvbl93b3JrZXJfc291cmNlLmpzIiwiLi4vLi4vLi4vc3JjL3NvdXJjZS93b3JrZXIuanMiXSwic291cmNlc0NvbnRlbnQiOlsiLy8gQGZsb3dcblxuaW1wb3J0IHR5cGUge0xheWVyU3BlY2lmaWNhdGlvbn0gZnJvbSAnLi90eXBlcy5qcyc7XG5cbmltcG9ydCByZWZQcm9wZXJ0aWVzIGZyb20gJy4vdXRpbC9yZWZfcHJvcGVydGllcy5qcyc7XG5cbmZ1bmN0aW9uIHN0cmluZ2lmeShvYmopIHtcbiAgICBpZiAodHlwZW9mIG9iaiA9PT0gJ251bWJlcicgfHwgdHlwZW9mIG9iaiA9PT0gJ2Jvb2xlYW4nIHx8IHR5cGVvZiBvYmogPT09ICdzdHJpbmcnIHx8IG9iaiA9PT0gdW5kZWZpbmVkIHx8IG9iaiA9PT0gbnVsbClcbiAgICAgICAgcmV0dXJuIEpTT04uc3RyaW5naWZ5KG9iaik7XG5cbiAgICBpZiAoQXJyYXkuaXNBcnJheShvYmopKSB7XG4gICAgICAgIGxldCBzdHIgPSAnWyc7XG4gICAgICAgIGZvciAoY29uc3QgdmFsIG9mIG9iaikge1xuICAgICAgICAgICAgc3RyICs9IGAke3N0cmluZ2lmeSh2YWwpfSxgO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiBgJHtzdHJ9XWA7XG4gICAgfVxuXG4gICAgbGV0IHN0ciA9ICd7JztcbiAgICBmb3IgKGNvbnN0IGtleSBvZiBPYmplY3Qua2V5cyhvYmopLnNvcnQoKSkge1xuICAgICAgICBzdHIgKz0gYCR7a2V5fToke3N0cmluZ2lmeSgob2JqOiBhbnkpW2tleV0pfSxgO1xuICAgIH1cbiAgICByZXR1cm4gYCR7c3RyfX1gO1xufVxuXG5mdW5jdGlvbiBnZXRLZXkobGF5ZXIpIHtcbiAgICBsZXQga2V5ID0gJyc7XG4gICAgZm9yIChjb25zdCBrIG9mIHJlZlByb3BlcnRpZXMpIHtcbiAgICAgICAga2V5ICs9IGAvJHtzdHJpbmdpZnkoKGxheWVyOiBhbnkpW2tdKX1gO1xuICAgIH1cbiAgICByZXR1cm4ga2V5O1xufVxuXG4vKipcbiAqIEdpdmVuIGFuIGFycmF5IG9mIGxheWVycywgcmV0dXJuIGFuIGFycmF5IG9mIGFycmF5cyBvZiBsYXllcnMgd2hlcmUgYWxsXG4gKiBsYXllcnMgaW4gZWFjaCBncm91cCBoYXZlIGlkZW50aWNhbCBsYXlvdXQtYWZmZWN0aW5nIHByb3BlcnRpZXMuIFRoZXNlXG4gKiBhcmUgdGhlIHByb3BlcnRpZXMgdGhhdCB3ZXJlIGZvcm1lcmx5IHVzZWQgYnkgZXhwbGljaXQgYHJlZmAgbWVjaGFuaXNtXG4gKiBmb3IgbGF5ZXJzOiAndHlwZScsICdzb3VyY2UnLCAnc291cmNlLWxheWVyJywgJ21pbnpvb20nLCAnbWF4em9vbScsXG4gKiAnZmlsdGVyJywgYW5kICdsYXlvdXQnLlxuICpcbiAqIFRoZSBpbnB1dCBpcyBub3QgbW9kaWZpZWQuIFRoZSBvdXRwdXQgbGF5ZXJzIGFyZSByZWZlcmVuY2VzIHRvIHRoZVxuICogaW5wdXQgbGF5ZXJzLlxuICpcbiAqIEBwcml2YXRlXG4gKiBAcGFyYW0ge0FycmF5PExheWVyPn0gbGF5ZXJzXG4gKiBAcGFyYW0ge09iamVjdH0gW2NhY2hlZEtleXNdIC0gYW4gb2JqZWN0IHRvIGtlZXAgYWxyZWFkeSBjYWxjdWxhdGVkIGtleXMuXG4gKiBAcmV0dXJucyB7QXJyYXk8QXJyYXk8TGF5ZXI+Pn1cbiAqL1xuZXhwb3J0IGRlZmF1bHQgZnVuY3Rpb24gZ3JvdXBCeUxheW91dChsYXllcnM6IEFycmF5PExheWVyU3BlY2lmaWNhdGlvbj4sIGNhY2hlZEtleXM6IHtbaWQ6IHN0cmluZ106IHN0cmluZ30pOiBBcnJheTxBcnJheTxMYXllclNwZWNpZmljYXRpb24+PiB7XG4gICAgY29uc3QgZ3JvdXBzID0ge307XG5cbiAgICBmb3IgKGxldCBpID0gMDsgaSA8IGxheWVycy5sZW5ndGg7IGkrKykge1xuXG4gICAgICAgIGNvbnN0IGsgPSAoY2FjaGVkS2V5cyAmJiBjYWNoZWRLZXlzW2xheWVyc1tpXS5pZF0pIHx8IGdldEtleShsYXllcnNbaV0pO1xuICAgICAgICAvLyB1cGRhdGUgdGhlIGNhY2hlIGlmIHRoZXJlIGlzIG9uZVxuICAgICAgICBpZiAoY2FjaGVkS2V5cylcbiAgICAgICAgICAgIGNhY2hlZEtleXNbbGF5ZXJzW2ldLmlkXSA9IGs7XG5cbiAgICAgICAgbGV0IGdyb3VwID0gZ3JvdXBzW2tdO1xuICAgICAgICBpZiAoIWdyb3VwKSB7XG4gICAgICAgICAgICBncm91cCA9IGdyb3Vwc1trXSA9IFtdO1xuICAgICAgICB9XG4gICAgICAgIGdyb3VwLnB1c2gobGF5ZXJzW2ldKTtcbiAgICB9XG5cbiAgICBjb25zdCByZXN1bHQgPSBbXTtcblxuICAgIGZvciAoY29uc3QgayBpbiBncm91cHMpIHtcbiAgICAgICAgcmVzdWx0LnB1c2goZ3JvdXBzW2tdKTtcbiAgICB9XG5cbiAgICByZXR1cm4gcmVzdWx0O1xufVxuIiwiLy8gQGZsb3dcblxuaW1wb3J0IGNyZWF0ZVN0eWxlTGF5ZXIgZnJvbSAnLi9jcmVhdGVfc3R5bGVfbGF5ZXIuanMnO1xuXG5pbXBvcnQge3ZhbHVlc30gZnJvbSAnLi4vdXRpbC91dGlsLmpzJztcbmltcG9ydCBncm91cEJ5TGF5b3V0IGZyb20gJy4uL3N0eWxlLXNwZWMvZ3JvdXBfYnlfbGF5b3V0LmpzJztcblxuaW1wb3J0IHR5cGUge1R5cGVkU3R5bGVMYXllcn0gZnJvbSAnLi9zdHlsZV9sYXllci90eXBlZF9zdHlsZV9sYXllci5qcyc7XG5pbXBvcnQgdHlwZSB7TGF5ZXJTcGVjaWZpY2F0aW9ufSBmcm9tICcuLi9zdHlsZS1zcGVjL3R5cGVzLmpzJztcblxuZXhwb3J0IHR5cGUgTGF5ZXJDb25maWdzID0ge1tfOiBzdHJpbmddOiBMYXllclNwZWNpZmljYXRpb24gfTtcbmV4cG9ydCB0eXBlIEZhbWlseTxMYXllcjogVHlwZWRTdHlsZUxheWVyPiA9IEFycmF5PExheWVyPjtcblxuY2xhc3MgU3R5bGVMYXllckluZGV4IHtcbiAgICBmYW1pbGllc0J5U291cmNlOiB7IFtzb3VyY2U6IHN0cmluZ106IHsgW3NvdXJjZUxheWVyOiBzdHJpbmddOiBBcnJheTxGYW1pbHk8VHlwZWRTdHlsZUxheWVyPj4gfSB9O1xuICAgIGtleUNhY2hlOiB7IFtzb3VyY2U6IHN0cmluZ106IHN0cmluZyB9O1xuXG4gICAgX2xheWVyQ29uZmlnczogTGF5ZXJDb25maWdzO1xuICAgIF9sYXllcnM6IHtbXzogc3RyaW5nXTogVHlwZWRTdHlsZUxheWVyIH07XG5cbiAgICBjb25zdHJ1Y3RvcihsYXllckNvbmZpZ3M6ID9BcnJheTxMYXllclNwZWNpZmljYXRpb24+KSB7XG4gICAgICAgIHRoaXMua2V5Q2FjaGUgPSB7fTtcbiAgICAgICAgaWYgKGxheWVyQ29uZmlncykge1xuICAgICAgICAgICAgdGhpcy5yZXBsYWNlKGxheWVyQ29uZmlncyk7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICByZXBsYWNlKGxheWVyQ29uZmlnczogQXJyYXk8TGF5ZXJTcGVjaWZpY2F0aW9uPikge1xuICAgICAgICB0aGlzLl9sYXllckNvbmZpZ3MgPSB7fTtcbiAgICAgICAgdGhpcy5fbGF5ZXJzID0ge307XG4gICAgICAgIHRoaXMudXBkYXRlKGxheWVyQ29uZmlncywgW10pO1xuICAgIH1cblxuICAgIHVwZGF0ZShsYXllckNvbmZpZ3M6IEFycmF5PExheWVyU3BlY2lmaWNhdGlvbj4sIHJlbW92ZWRJZHM6IEFycmF5PHN0cmluZz4pIHtcbiAgICAgICAgZm9yIChjb25zdCBsYXllckNvbmZpZyBvZiBsYXllckNvbmZpZ3MpIHtcbiAgICAgICAgICAgIHRoaXMuX2xheWVyQ29uZmlnc1tsYXllckNvbmZpZy5pZF0gPSBsYXllckNvbmZpZztcblxuICAgICAgICAgICAgY29uc3QgbGF5ZXIgPSB0aGlzLl9sYXllcnNbbGF5ZXJDb25maWcuaWRdID0gKChjcmVhdGVTdHlsZUxheWVyKGxheWVyQ29uZmlnKTogYW55KTogVHlwZWRTdHlsZUxheWVyKTtcbiAgICAgICAgICAgIGxheWVyLmNvbXBpbGVGaWx0ZXIoKTtcbiAgICAgICAgICAgIGlmICh0aGlzLmtleUNhY2hlW2xheWVyQ29uZmlnLmlkXSlcbiAgICAgICAgICAgICAgICBkZWxldGUgdGhpcy5rZXlDYWNoZVtsYXllckNvbmZpZy5pZF07XG4gICAgICAgIH1cbiAgICAgICAgZm9yIChjb25zdCBpZCBvZiByZW1vdmVkSWRzKSB7XG4gICAgICAgICAgICBkZWxldGUgdGhpcy5rZXlDYWNoZVtpZF07XG4gICAgICAgICAgICBkZWxldGUgdGhpcy5fbGF5ZXJDb25maWdzW2lkXTtcbiAgICAgICAgICAgIGRlbGV0ZSB0aGlzLl9sYXllcnNbaWRdO1xuICAgICAgICB9XG5cbiAgICAgICAgdGhpcy5mYW1pbGllc0J5U291cmNlID0ge307XG5cbiAgICAgICAgY29uc3QgZ3JvdXBzID0gZ3JvdXBCeUxheW91dCh2YWx1ZXModGhpcy5fbGF5ZXJDb25maWdzKSwgdGhpcy5rZXlDYWNoZSk7XG5cbiAgICAgICAgZm9yIChjb25zdCBsYXllckNvbmZpZ3Mgb2YgZ3JvdXBzKSB7XG4gICAgICAgICAgICBjb25zdCBsYXllcnMgPSBsYXllckNvbmZpZ3MubWFwKChsYXllckNvbmZpZykgPT4gdGhpcy5fbGF5ZXJzW2xheWVyQ29uZmlnLmlkXSk7XG5cbiAgICAgICAgICAgIGNvbnN0IGxheWVyID0gbGF5ZXJzWzBdO1xuICAgICAgICAgICAgaWYgKGxheWVyLnZpc2liaWxpdHkgPT09ICdub25lJykge1xuICAgICAgICAgICAgICAgIGNvbnRpbnVlO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICBjb25zdCBzb3VyY2VJZCA9IGxheWVyLnNvdXJjZSB8fCAnJztcbiAgICAgICAgICAgIGxldCBzb3VyY2VHcm91cCA9IHRoaXMuZmFtaWxpZXNCeVNvdXJjZVtzb3VyY2VJZF07XG4gICAgICAgICAgICBpZiAoIXNvdXJjZUdyb3VwKSB7XG4gICAgICAgICAgICAgICAgc291cmNlR3JvdXAgPSB0aGlzLmZhbWlsaWVzQnlTb3VyY2Vbc291cmNlSWRdID0ge307XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgICAgIGNvbnN0IHNvdXJjZUxheWVySWQgPSBsYXllci5zb3VyY2VMYXllciB8fCAnX2dlb2pzb25UaWxlTGF5ZXInO1xuICAgICAgICAgICAgbGV0IHNvdXJjZUxheWVyRmFtaWxpZXMgPSBzb3VyY2VHcm91cFtzb3VyY2VMYXllcklkXTtcbiAgICAgICAgICAgIGlmICghc291cmNlTGF5ZXJGYW1pbGllcykge1xuICAgICAgICAgICAgICAgIHNvdXJjZUxheWVyRmFtaWxpZXMgPSBzb3VyY2VHcm91cFtzb3VyY2VMYXllcklkXSA9IFtdO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICBzb3VyY2VMYXllckZhbWlsaWVzLnB1c2gobGF5ZXJzKTtcbiAgICAgICAgfVxuICAgIH1cbn1cblxuZXhwb3J0IGRlZmF1bHQgU3R5bGVMYXllckluZGV4O1xuIiwiLy8gQGZsb3dcblxuaW1wb3J0IERFTURhdGEgZnJvbSAnLi4vZGF0YS9kZW1fZGF0YS5qcyc7XG5pbXBvcnQgd2luZG93IGZyb20gJy4uL3V0aWwvd2luZG93LmpzJztcblxuaW1wb3J0IHR5cGUgQWN0b3IgZnJvbSAnLi4vdXRpbC9hY3Rvci5qcyc7XG5pbXBvcnQgdHlwZSB7V29ya2VyREVNVGlsZVBhcmFtZXRlcnMsIFdvcmtlckRFTVRpbGVDYWxsYmFja30gZnJvbSAnLi93b3JrZXJfc291cmNlLmpzJztcblxuY2xhc3MgUmFzdGVyREVNVGlsZVdvcmtlclNvdXJjZSB7XG4gICAgYWN0b3I6IEFjdG9yO1xuICAgIG9mZnNjcmVlbkNhbnZhczogT2Zmc2NyZWVuQ2FudmFzO1xuICAgIG9mZnNjcmVlbkNhbnZhc0NvbnRleHQ6IENhbnZhc1JlbmRlcmluZ0NvbnRleHQyRDtcblxuICAgIGxvYWRUaWxlKHBhcmFtczogV29ya2VyREVNVGlsZVBhcmFtZXRlcnMsIGNhbGxiYWNrOiBXb3JrZXJERU1UaWxlQ2FsbGJhY2spIHtcbiAgICAgICAgY29uc3Qge3VpZCwgZW5jb2RpbmcsIHJhd0ltYWdlRGF0YSwgcGFkZGluZywgYnVpbGRRdWFkVHJlZX0gPSBwYXJhbXM7XG4gICAgICAgIC8vIE1haW4gdGhyZWFkIHdpbGwgdHJhbnNmZXIgSW1hZ2VCaXRtYXAgaWYgb2Zmc2NyZWVuIGRlY29kZSB3aXRoIE9mZnNjcmVlbkNhbnZhcyBpcyBzdXBwb3J0ZWQsIGVsc2UgaXQgd2lsbCB0cmFuc2ZlciBhbiBhbHJlYWR5IGRlY29kZWQgaW1hZ2UuXG4gICAgICAgIC8vIEZsb3cgc3RydWdnbGVzIHRvIHJlZmluZSBJbWFnZUJpdG1hcCB0eXBlLCBsaWtlbHkgZHVlIHRvIHRoZSBKU0RvbSBzaGltXG4gICAgICAgIGNvbnN0IGltYWdlUGl4ZWxzID0gd2luZG93LkltYWdlQml0bWFwICYmIHJhd0ltYWdlRGF0YSBpbnN0YW5jZW9mIHdpbmRvdy5JbWFnZUJpdG1hcCA/IHRoaXMuZ2V0SW1hZ2VEYXRhKHJhd0ltYWdlRGF0YSwgcGFkZGluZykgOiAoKHJhd0ltYWdlRGF0YTogYW55KTogSW1hZ2VEYXRhKTtcbiAgICAgICAgY29uc3QgZGVtID0gbmV3IERFTURhdGEodWlkLCBpbWFnZVBpeGVscywgZW5jb2RpbmcsIHBhZGRpbmcgPCAxLCBidWlsZFF1YWRUcmVlKTtcbiAgICAgICAgY2FsbGJhY2sobnVsbCwgZGVtKTtcbiAgICB9XG5cbiAgICBnZXRJbWFnZURhdGEoaW1nQml0bWFwOiBJbWFnZUJpdG1hcCwgcGFkZGluZzogbnVtYmVyKTogSW1hZ2VEYXRhIHtcbiAgICAgICAgLy8gTGF6aWx5IGluaXRpYWxpemUgT2Zmc2NyZWVuQ2FudmFzXG4gICAgICAgIGlmICghdGhpcy5vZmZzY3JlZW5DYW52YXMgfHwgIXRoaXMub2Zmc2NyZWVuQ2FudmFzQ29udGV4dCkge1xuICAgICAgICAgICAgLy8gRGVtIHRpbGVzIGFyZSB0eXBpY2FsbHkgMjU2eDI1NlxuICAgICAgICAgICAgdGhpcy5vZmZzY3JlZW5DYW52YXMgPSBuZXcgT2Zmc2NyZWVuQ2FudmFzKGltZ0JpdG1hcC53aWR0aCwgaW1nQml0bWFwLmhlaWdodCk7XG4gICAgICAgICAgICB0aGlzLm9mZnNjcmVlbkNhbnZhc0NvbnRleHQgPSB0aGlzLm9mZnNjcmVlbkNhbnZhcy5nZXRDb250ZXh0KCcyZCcpO1xuICAgICAgICB9XG5cbiAgICAgICAgdGhpcy5vZmZzY3JlZW5DYW52YXMud2lkdGggPSBpbWdCaXRtYXAud2lkdGg7XG4gICAgICAgIHRoaXMub2Zmc2NyZWVuQ2FudmFzLmhlaWdodCA9IGltZ0JpdG1hcC5oZWlnaHQ7XG5cbiAgICAgICAgdGhpcy5vZmZzY3JlZW5DYW52YXNDb250ZXh0LmRyYXdJbWFnZShpbWdCaXRtYXAsIDAsIDAsIGltZ0JpdG1hcC53aWR0aCwgaW1nQml0bWFwLmhlaWdodCk7XG4gICAgICAgIC8vIEluc2VydCBvciByZW1vdmUgZGVmaW5lZCBwYWRkaW5nIGFyb3VuZCB0aGUgaW1hZ2UgdG8gYWxsb3cgYmFja2ZpbGxpbmcgZm9yIG5laWdoYm9yaW5nIGRhdGEuXG4gICAgICAgIGNvbnN0IGltZ0RhdGEgPSB0aGlzLm9mZnNjcmVlbkNhbnZhc0NvbnRleHQuZ2V0SW1hZ2VEYXRhKC1wYWRkaW5nLCAtcGFkZGluZywgaW1nQml0bWFwLndpZHRoICsgMiAqIHBhZGRpbmcsIGltZ0JpdG1hcC5oZWlnaHQgKyAyICogcGFkZGluZyk7XG4gICAgICAgIHRoaXMub2Zmc2NyZWVuQ2FudmFzQ29udGV4dC5jbGVhclJlY3QoMCwgMCwgdGhpcy5vZmZzY3JlZW5DYW52YXMud2lkdGgsIHRoaXMub2Zmc2NyZWVuQ2FudmFzLmhlaWdodCk7XG4gICAgICAgIHJldHVybiBpbWdEYXRhO1xuICAgIH1cbn1cblxuZXhwb3J0IGRlZmF1bHQgUmFzdGVyREVNVGlsZVdvcmtlclNvdXJjZTtcbiIsIlxubW9kdWxlLmV4cG9ydHMgPSByZXdpbmQ7XG5cbmZ1bmN0aW9uIHJld2luZChnaiwgb3V0ZXIpIHtcbiAgICB2YXIgdHlwZSA9IGdqICYmIGdqLnR5cGUsIGk7XG5cbiAgICBpZiAodHlwZSA9PT0gJ0ZlYXR1cmVDb2xsZWN0aW9uJykge1xuICAgICAgICBmb3IgKGkgPSAwOyBpIDwgZ2ouZmVhdHVyZXMubGVuZ3RoOyBpKyspIHJld2luZChnai5mZWF0dXJlc1tpXSwgb3V0ZXIpO1xuXG4gICAgfSBlbHNlIGlmICh0eXBlID09PSAnR2VvbWV0cnlDb2xsZWN0aW9uJykge1xuICAgICAgICBmb3IgKGkgPSAwOyBpIDwgZ2ouZ2VvbWV0cmllcy5sZW5ndGg7IGkrKykgcmV3aW5kKGdqLmdlb21ldHJpZXNbaV0sIG91dGVyKTtcblxuICAgIH0gZWxzZSBpZiAodHlwZSA9PT0gJ0ZlYXR1cmUnKSB7XG4gICAgICAgIHJld2luZChnai5nZW9tZXRyeSwgb3V0ZXIpO1xuXG4gICAgfSBlbHNlIGlmICh0eXBlID09PSAnUG9seWdvbicpIHtcbiAgICAgICAgcmV3aW5kUmluZ3MoZ2ouY29vcmRpbmF0ZXMsIG91dGVyKTtcblxuICAgIH0gZWxzZSBpZiAodHlwZSA9PT0gJ011bHRpUG9seWdvbicpIHtcbiAgICAgICAgZm9yIChpID0gMDsgaSA8IGdqLmNvb3JkaW5hdGVzLmxlbmd0aDsgaSsrKSByZXdpbmRSaW5ncyhnai5jb29yZGluYXRlc1tpXSwgb3V0ZXIpO1xuICAgIH1cblxuICAgIHJldHVybiBnajtcbn1cblxuZnVuY3Rpb24gcmV3aW5kUmluZ3MocmluZ3MsIG91dGVyKSB7XG4gICAgaWYgKHJpbmdzLmxlbmd0aCA9PT0gMCkgcmV0dXJuO1xuXG4gICAgcmV3aW5kUmluZyhyaW5nc1swXSwgb3V0ZXIpO1xuICAgIGZvciAodmFyIGkgPSAxOyBpIDwgcmluZ3MubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgcmV3aW5kUmluZyhyaW5nc1tpXSwgIW91dGVyKTtcbiAgICB9XG59XG5cbmZ1bmN0aW9uIHJld2luZFJpbmcocmluZywgZGlyKSB7XG4gICAgdmFyIGFyZWEgPSAwLCBlcnIgPSAwO1xuICAgIGZvciAodmFyIGkgPSAwLCBsZW4gPSByaW5nLmxlbmd0aCwgaiA9IGxlbiAtIDE7IGkgPCBsZW47IGogPSBpKyspIHtcbiAgICAgICAgdmFyIGsgPSAocmluZ1tpXVswXSAtIHJpbmdbal1bMF0pICogKHJpbmdbal1bMV0gKyByaW5nW2ldWzFdKTtcbiAgICAgICAgdmFyIG0gPSBhcmVhICsgaztcbiAgICAgICAgZXJyICs9IE1hdGguYWJzKGFyZWEpID49IE1hdGguYWJzKGspID8gYXJlYSAtIG0gKyBrIDogayAtIG0gKyBhcmVhO1xuICAgICAgICBhcmVhID0gbTtcbiAgICB9XG4gICAgaWYgKGFyZWEgKyBlcnIgPj0gMCAhPT0gISFkaXIpIHJpbmcucmV2ZXJzZSgpO1xufVxuIiwiLy8gQGZsb3dcblxuaW1wb3J0IFBvaW50IGZyb20gJ0BtYXBib3gvcG9pbnQtZ2VvbWV0cnknO1xuXG5pbXBvcnQgbXZ0IGZyb20gJ0BtYXBib3gvdmVjdG9yLXRpbGUnO1xuY29uc3QgdG9HZW9KU09OID0gbXZ0LlZlY3RvclRpbGVGZWF0dXJlLnByb3RvdHlwZS50b0dlb0pTT047XG5pbXBvcnQgRVhURU5UIGZyb20gJy4uL2RhdGEvZXh0ZW50LmpzJztcblxuLy8gVGhlIGZlYXR1cmUgdHlwZSB1c2VkIGJ5IGdlb2pzb24tdnQgYW5kIHN1cGVyY2x1c3Rlci4gU2hvdWxkIGJlIGV4dHJhY3RlZCB0b1xuLy8gZ2xvYmFsIHR5cGUgYW5kIHVzZWQgaW4gbW9kdWxlIGRlZmluaXRpb25zIGZvciB0aG9zZSB0d28gbW9kdWxlcy5cbnR5cGUgRmVhdHVyZSA9IHtcbiAgICB0eXBlOiAxLFxuICAgIGlkOiBtaXhlZCxcbiAgICB0YWdzOiB7W186IHN0cmluZ106IHN0cmluZyB8IG51bWJlciB8IGJvb2xlYW59LFxuICAgIGdlb21ldHJ5OiBBcnJheTxbbnVtYmVyLCBudW1iZXJdPixcbn0gfCB7XG4gICAgdHlwZTogMiB8IDMsXG4gICAgaWQ6IG1peGVkLFxuICAgIHRhZ3M6IHtbXzogc3RyaW5nXTogc3RyaW5nIHwgbnVtYmVyIHwgYm9vbGVhbn0sXG4gICAgZ2VvbWV0cnk6IEFycmF5PEFycmF5PFtudW1iZXIsIG51bWJlcl0+Pixcbn1cblxuY2xhc3MgRmVhdHVyZVdyYXBwZXIgaW1wbGVtZW50cyBWZWN0b3JUaWxlRmVhdHVyZSB7XG4gICAgX2ZlYXR1cmU6IEZlYXR1cmU7XG5cbiAgICBleHRlbnQ6IG51bWJlcjtcbiAgICB0eXBlOiAxIHwgMiB8IDM7XG4gICAgaWQ6IG51bWJlcjtcbiAgICBwcm9wZXJ0aWVzOiB7W186IHN0cmluZ106IHN0cmluZyB8IG51bWJlciB8IGJvb2xlYW59O1xuXG4gICAgY29uc3RydWN0b3IoZmVhdHVyZTogRmVhdHVyZSkge1xuICAgICAgICB0aGlzLl9mZWF0dXJlID0gZmVhdHVyZTtcblxuICAgICAgICB0aGlzLmV4dGVudCA9IEVYVEVOVDtcbiAgICAgICAgdGhpcy50eXBlID0gZmVhdHVyZS50eXBlO1xuICAgICAgICB0aGlzLnByb3BlcnRpZXMgPSBmZWF0dXJlLnRhZ3M7XG5cbiAgICAgICAgLy8gSWYgdGhlIGZlYXR1cmUgaGFzIGEgdG9wLWxldmVsIGBpZGAgcHJvcGVydHksIGNvcHkgaXQgb3ZlciwgYnV0IG9ubHlcbiAgICAgICAgLy8gaWYgaXQgY2FuIGJlIGNvZXJjZWQgdG8gYW4gaW50ZWdlciwgYmVjYXVzZSB0aGlzIHdyYXBwZXIgaXMgdXNlZCBmb3JcbiAgICAgICAgLy8gc2VyaWFsaXppbmcgZ2VvanNvbiBmZWF0dXJlIGRhdGEgaW50byB2ZWN0b3IgdGlsZSBQQkYgZGF0YSwgYW5kIHRoZVxuICAgICAgICAvLyB2ZWN0b3IgdGlsZSBzcGVjIG9ubHkgc3VwcG9ydHMgaW50ZWdlciB2YWx1ZXMgZm9yIGZlYXR1cmUgaWRzIC0tXG4gICAgICAgIC8vIGFsbG93aW5nIG5vbi1pbnRlZ2VyIHZhbHVlcyBoZXJlIHJlc3VsdHMgaW4gYSBub24tY29tcGxpYW50IFBCRlxuICAgICAgICAvLyB0aGF0IGNhdXNlcyBhbiBleGNlcHRpb24gd2hlbiBpdCBpcyBwYXJzZWQgd2l0aCB2ZWN0b3ItdGlsZS1qc1xuICAgICAgICBpZiAoJ2lkJyBpbiBmZWF0dXJlICYmICFpc05hTihmZWF0dXJlLmlkKSkge1xuICAgICAgICAgICAgdGhpcy5pZCA9IHBhcnNlSW50KGZlYXR1cmUuaWQsIDEwKTtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIGxvYWRHZW9tZXRyeSgpIHtcbiAgICAgICAgaWYgKHRoaXMuX2ZlYXR1cmUudHlwZSA9PT0gMSkge1xuICAgICAgICAgICAgY29uc3QgZ2VvbWV0cnkgPSBbXTtcbiAgICAgICAgICAgIGZvciAoY29uc3QgcG9pbnQgb2YgdGhpcy5fZmVhdHVyZS5nZW9tZXRyeSkge1xuICAgICAgICAgICAgICAgIGdlb21ldHJ5LnB1c2goW25ldyBQb2ludChwb2ludFswXSwgcG9pbnRbMV0pXSk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICByZXR1cm4gZ2VvbWV0cnk7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBjb25zdCBnZW9tZXRyeSA9IFtdO1xuICAgICAgICAgICAgZm9yIChjb25zdCByaW5nIG9mIHRoaXMuX2ZlYXR1cmUuZ2VvbWV0cnkpIHtcbiAgICAgICAgICAgICAgICBjb25zdCBuZXdSaW5nID0gW107XG4gICAgICAgICAgICAgICAgZm9yIChjb25zdCBwb2ludCBvZiByaW5nKSB7XG4gICAgICAgICAgICAgICAgICAgIG5ld1JpbmcucHVzaChuZXcgUG9pbnQocG9pbnRbMF0sIHBvaW50WzFdKSk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIGdlb21ldHJ5LnB1c2gobmV3UmluZyk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICByZXR1cm4gZ2VvbWV0cnk7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICB0b0dlb0pTT04oeDogbnVtYmVyLCB5OiBudW1iZXIsIHo6IG51bWJlcikge1xuICAgICAgICByZXR1cm4gdG9HZW9KU09OLmNhbGwodGhpcywgeCwgeSwgeik7XG4gICAgfVxufVxuXG5jbGFzcyBHZW9KU09OV3JhcHBlciBpbXBsZW1lbnRzIFZlY3RvclRpbGUsIFZlY3RvclRpbGVMYXllciB7XG4gICAgbGF5ZXJzOiB7W186IHN0cmluZ106IFZlY3RvclRpbGVMYXllcn07XG4gICAgbmFtZTogc3RyaW5nO1xuICAgIGV4dGVudDogbnVtYmVyO1xuICAgIGxlbmd0aDogbnVtYmVyO1xuICAgIF9mZWF0dXJlczogQXJyYXk8RmVhdHVyZT47XG5cbiAgICBjb25zdHJ1Y3RvcihmZWF0dXJlczogQXJyYXk8RmVhdHVyZT4pIHtcbiAgICAgICAgdGhpcy5sYXllcnMgPSB7J19nZW9qc29uVGlsZUxheWVyJzogdGhpc307XG4gICAgICAgIHRoaXMubmFtZSA9ICdfZ2VvanNvblRpbGVMYXllcic7XG4gICAgICAgIHRoaXMuZXh0ZW50ID0gRVhURU5UO1xuICAgICAgICB0aGlzLmxlbmd0aCA9IGZlYXR1cmVzLmxlbmd0aDtcbiAgICAgICAgdGhpcy5fZmVhdHVyZXMgPSBmZWF0dXJlcztcbiAgICB9XG5cbiAgICBmZWF0dXJlKGk6IG51bWJlcik6IFZlY3RvclRpbGVGZWF0dXJlIHtcbiAgICAgICAgcmV0dXJuIG5ldyBGZWF0dXJlV3JhcHBlcih0aGlzLl9mZWF0dXJlc1tpXSk7XG4gICAgfVxufVxuXG5leHBvcnQgZGVmYXVsdCBHZW9KU09OV3JhcHBlcjtcbiIsIid1c2Ugc3RyaWN0J1xuXG52YXIgUG9pbnQgPSByZXF1aXJlKCdAbWFwYm94L3BvaW50LWdlb21ldHJ5JylcbnZhciBWZWN0b3JUaWxlRmVhdHVyZSA9IHJlcXVpcmUoJ0BtYXBib3gvdmVjdG9yLXRpbGUnKS5WZWN0b3JUaWxlRmVhdHVyZVxuXG5tb2R1bGUuZXhwb3J0cyA9IEdlb0pTT05XcmFwcGVyXG5cbi8vIGNvbmZvcm0gdG8gdmVjdG9ydGlsZSBhcGlcbmZ1bmN0aW9uIEdlb0pTT05XcmFwcGVyIChmZWF0dXJlcywgb3B0aW9ucykge1xuICB0aGlzLm9wdGlvbnMgPSBvcHRpb25zIHx8IHt9XG4gIHRoaXMuZmVhdHVyZXMgPSBmZWF0dXJlc1xuICB0aGlzLmxlbmd0aCA9IGZlYXR1cmVzLmxlbmd0aFxufVxuXG5HZW9KU09OV3JhcHBlci5wcm90b3R5cGUuZmVhdHVyZSA9IGZ1bmN0aW9uIChpKSB7XG4gIHJldHVybiBuZXcgRmVhdHVyZVdyYXBwZXIodGhpcy5mZWF0dXJlc1tpXSwgdGhpcy5vcHRpb25zLmV4dGVudClcbn1cblxuZnVuY3Rpb24gRmVhdHVyZVdyYXBwZXIgKGZlYXR1cmUsIGV4dGVudCkge1xuICB0aGlzLmlkID0gdHlwZW9mIGZlYXR1cmUuaWQgPT09ICdudW1iZXInID8gZmVhdHVyZS5pZCA6IHVuZGVmaW5lZFxuICB0aGlzLnR5cGUgPSBmZWF0dXJlLnR5cGVcbiAgdGhpcy5yYXdHZW9tZXRyeSA9IGZlYXR1cmUudHlwZSA9PT0gMSA/IFtmZWF0dXJlLmdlb21ldHJ5XSA6IGZlYXR1cmUuZ2VvbWV0cnlcbiAgdGhpcy5wcm9wZXJ0aWVzID0gZmVhdHVyZS50YWdzXG4gIHRoaXMuZXh0ZW50ID0gZXh0ZW50IHx8IDQwOTZcbn1cblxuRmVhdHVyZVdyYXBwZXIucHJvdG90eXBlLmxvYWRHZW9tZXRyeSA9IGZ1bmN0aW9uICgpIHtcbiAgdmFyIHJpbmdzID0gdGhpcy5yYXdHZW9tZXRyeVxuICB0aGlzLmdlb21ldHJ5ID0gW11cblxuICBmb3IgKHZhciBpID0gMDsgaSA8IHJpbmdzLmxlbmd0aDsgaSsrKSB7XG4gICAgdmFyIHJpbmcgPSByaW5nc1tpXVxuICAgIHZhciBuZXdSaW5nID0gW11cbiAgICBmb3IgKHZhciBqID0gMDsgaiA8IHJpbmcubGVuZ3RoOyBqKyspIHtcbiAgICAgIG5ld1JpbmcucHVzaChuZXcgUG9pbnQocmluZ1tqXVswXSwgcmluZ1tqXVsxXSkpXG4gICAgfVxuICAgIHRoaXMuZ2VvbWV0cnkucHVzaChuZXdSaW5nKVxuICB9XG4gIHJldHVybiB0aGlzLmdlb21ldHJ5XG59XG5cbkZlYXR1cmVXcmFwcGVyLnByb3RvdHlwZS5iYm94ID0gZnVuY3Rpb24gKCkge1xuICBpZiAoIXRoaXMuZ2VvbWV0cnkpIHRoaXMubG9hZEdlb21ldHJ5KClcblxuICB2YXIgcmluZ3MgPSB0aGlzLmdlb21ldHJ5XG4gIHZhciB4MSA9IEluZmluaXR5XG4gIHZhciB4MiA9IC1JbmZpbml0eVxuICB2YXIgeTEgPSBJbmZpbml0eVxuICB2YXIgeTIgPSAtSW5maW5pdHlcblxuICBmb3IgKHZhciBpID0gMDsgaSA8IHJpbmdzLmxlbmd0aDsgaSsrKSB7XG4gICAgdmFyIHJpbmcgPSByaW5nc1tpXVxuXG4gICAgZm9yICh2YXIgaiA9IDA7IGogPCByaW5nLmxlbmd0aDsgaisrKSB7XG4gICAgICB2YXIgY29vcmQgPSByaW5nW2pdXG5cbiAgICAgIHgxID0gTWF0aC5taW4oeDEsIGNvb3JkLngpXG4gICAgICB4MiA9IE1hdGgubWF4KHgyLCBjb29yZC54KVxuICAgICAgeTEgPSBNYXRoLm1pbih5MSwgY29vcmQueSlcbiAgICAgIHkyID0gTWF0aC5tYXgoeTIsIGNvb3JkLnkpXG4gICAgfVxuICB9XG5cbiAgcmV0dXJuIFt4MSwgeTEsIHgyLCB5Ml1cbn1cblxuRmVhdHVyZVdyYXBwZXIucHJvdG90eXBlLnRvR2VvSlNPTiA9IFZlY3RvclRpbGVGZWF0dXJlLnByb3RvdHlwZS50b0dlb0pTT05cbiIsInZhciBQYmYgPSByZXF1aXJlKCdwYmYnKVxudmFyIEdlb0pTT05XcmFwcGVyID0gcmVxdWlyZSgnLi9saWIvZ2VvanNvbl93cmFwcGVyJylcblxubW9kdWxlLmV4cG9ydHMgPSBmcm9tVmVjdG9yVGlsZUpzXG5tb2R1bGUuZXhwb3J0cy5mcm9tVmVjdG9yVGlsZUpzID0gZnJvbVZlY3RvclRpbGVKc1xubW9kdWxlLmV4cG9ydHMuZnJvbUdlb2pzb25WdCA9IGZyb21HZW9qc29uVnRcbm1vZHVsZS5leHBvcnRzLkdlb0pTT05XcmFwcGVyID0gR2VvSlNPTldyYXBwZXJcblxuLyoqXG4gKiBTZXJpYWxpemUgYSB2ZWN0b3ItdGlsZS1qcy1jcmVhdGVkIHRpbGUgdG8gcGJmXG4gKlxuICogQHBhcmFtIHtPYmplY3R9IHRpbGVcbiAqIEByZXR1cm4ge0J1ZmZlcn0gdW5jb21wcmVzc2VkLCBwYmYtc2VyaWFsaXplZCB0aWxlIGRhdGFcbiAqL1xuZnVuY3Rpb24gZnJvbVZlY3RvclRpbGVKcyAodGlsZSkge1xuICB2YXIgb3V0ID0gbmV3IFBiZigpXG4gIHdyaXRlVGlsZSh0aWxlLCBvdXQpXG4gIHJldHVybiBvdXQuZmluaXNoKClcbn1cblxuLyoqXG4gKiBTZXJpYWxpemVkIGEgZ2VvanNvbi12dC1jcmVhdGVkIHRpbGUgdG8gcGJmLlxuICpcbiAqIEBwYXJhbSB7T2JqZWN0fSBsYXllcnMgLSBBbiBvYmplY3QgbWFwcGluZyBsYXllciBuYW1lcyB0byBnZW9qc29uLXZ0LWNyZWF0ZWQgdmVjdG9yIHRpbGUgb2JqZWN0c1xuICogQHBhcmFtIHtPYmplY3R9IFtvcHRpb25zXSAtIEFuIG9iamVjdCBzcGVjaWZ5aW5nIHRoZSB2ZWN0b3ItdGlsZSBzcGVjaWZpY2F0aW9uIHZlcnNpb24gYW5kIGV4dGVudCB0aGF0IHdlcmUgdXNlZCB0byBjcmVhdGUgYGxheWVyc2AuXG4gKiBAcGFyYW0ge051bWJlcn0gW29wdGlvbnMudmVyc2lvbj0xXSAtIFZlcnNpb24gb2YgdmVjdG9yLXRpbGUgc3BlYyB1c2VkXG4gKiBAcGFyYW0ge051bWJlcn0gW29wdGlvbnMuZXh0ZW50PTQwOTZdIC0gRXh0ZW50IG9mIHRoZSB2ZWN0b3IgdGlsZVxuICogQHJldHVybiB7QnVmZmVyfSB1bmNvbXByZXNzZWQsIHBiZi1zZXJpYWxpemVkIHRpbGUgZGF0YVxuICovXG5mdW5jdGlvbiBmcm9tR2VvanNvblZ0IChsYXllcnMsIG9wdGlvbnMpIHtcbiAgb3B0aW9ucyA9IG9wdGlvbnMgfHwge31cbiAgdmFyIGwgPSB7fVxuICBmb3IgKHZhciBrIGluIGxheWVycykge1xuICAgIGxba10gPSBuZXcgR2VvSlNPTldyYXBwZXIobGF5ZXJzW2tdLmZlYXR1cmVzLCBvcHRpb25zKVxuICAgIGxba10ubmFtZSA9IGtcbiAgICBsW2tdLnZlcnNpb24gPSBvcHRpb25zLnZlcnNpb25cbiAgICBsW2tdLmV4dGVudCA9IG9wdGlvbnMuZXh0ZW50XG4gIH1cbiAgcmV0dXJuIGZyb21WZWN0b3JUaWxlSnMoeyBsYXllcnM6IGwgfSlcbn1cblxuZnVuY3Rpb24gd3JpdGVUaWxlICh0aWxlLCBwYmYpIHtcbiAgZm9yICh2YXIga2V5IGluIHRpbGUubGF5ZXJzKSB7XG4gICAgcGJmLndyaXRlTWVzc2FnZSgzLCB3cml0ZUxheWVyLCB0aWxlLmxheWVyc1trZXldKVxuICB9XG59XG5cbmZ1bmN0aW9uIHdyaXRlTGF5ZXIgKGxheWVyLCBwYmYpIHtcbiAgcGJmLndyaXRlVmFyaW50RmllbGQoMTUsIGxheWVyLnZlcnNpb24gfHwgMSlcbiAgcGJmLndyaXRlU3RyaW5nRmllbGQoMSwgbGF5ZXIubmFtZSB8fCAnJylcbiAgcGJmLndyaXRlVmFyaW50RmllbGQoNSwgbGF5ZXIuZXh0ZW50IHx8IDQwOTYpXG5cbiAgdmFyIGlcbiAgdmFyIGNvbnRleHQgPSB7XG4gICAga2V5czogW10sXG4gICAgdmFsdWVzOiBbXSxcbiAgICBrZXljYWNoZToge30sXG4gICAgdmFsdWVjYWNoZToge31cbiAgfVxuXG4gIGZvciAoaSA9IDA7IGkgPCBsYXllci5sZW5ndGg7IGkrKykge1xuICAgIGNvbnRleHQuZmVhdHVyZSA9IGxheWVyLmZlYXR1cmUoaSlcbiAgICBwYmYud3JpdGVNZXNzYWdlKDIsIHdyaXRlRmVhdHVyZSwgY29udGV4dClcbiAgfVxuXG4gIHZhciBrZXlzID0gY29udGV4dC5rZXlzXG4gIGZvciAoaSA9IDA7IGkgPCBrZXlzLmxlbmd0aDsgaSsrKSB7XG4gICAgcGJmLndyaXRlU3RyaW5nRmllbGQoMywga2V5c1tpXSlcbiAgfVxuXG4gIHZhciB2YWx1ZXMgPSBjb250ZXh0LnZhbHVlc1xuICBmb3IgKGkgPSAwOyBpIDwgdmFsdWVzLmxlbmd0aDsgaSsrKSB7XG4gICAgcGJmLndyaXRlTWVzc2FnZSg0LCB3cml0ZVZhbHVlLCB2YWx1ZXNbaV0pXG4gIH1cbn1cblxuZnVuY3Rpb24gd3JpdGVGZWF0dXJlIChjb250ZXh0LCBwYmYpIHtcbiAgdmFyIGZlYXR1cmUgPSBjb250ZXh0LmZlYXR1cmVcblxuICBpZiAoZmVhdHVyZS5pZCAhPT0gdW5kZWZpbmVkKSB7XG4gICAgcGJmLndyaXRlVmFyaW50RmllbGQoMSwgZmVhdHVyZS5pZClcbiAgfVxuXG4gIHBiZi53cml0ZU1lc3NhZ2UoMiwgd3JpdGVQcm9wZXJ0aWVzLCBjb250ZXh0KVxuICBwYmYud3JpdGVWYXJpbnRGaWVsZCgzLCBmZWF0dXJlLnR5cGUpXG4gIHBiZi53cml0ZU1lc3NhZ2UoNCwgd3JpdGVHZW9tZXRyeSwgZmVhdHVyZSlcbn1cblxuZnVuY3Rpb24gd3JpdGVQcm9wZXJ0aWVzIChjb250ZXh0LCBwYmYpIHtcbiAgdmFyIGZlYXR1cmUgPSBjb250ZXh0LmZlYXR1cmVcbiAgdmFyIGtleXMgPSBjb250ZXh0LmtleXNcbiAgdmFyIHZhbHVlcyA9IGNvbnRleHQudmFsdWVzXG4gIHZhciBrZXljYWNoZSA9IGNvbnRleHQua2V5Y2FjaGVcbiAgdmFyIHZhbHVlY2FjaGUgPSBjb250ZXh0LnZhbHVlY2FjaGVcblxuICBmb3IgKHZhciBrZXkgaW4gZmVhdHVyZS5wcm9wZXJ0aWVzKSB7XG4gICAgdmFyIHZhbHVlID0gZmVhdHVyZS5wcm9wZXJ0aWVzW2tleV1cblxuICAgIHZhciBrZXlJbmRleCA9IGtleWNhY2hlW2tleV1cbiAgICBpZiAodmFsdWUgPT09IG51bGwpIGNvbnRpbnVlIC8vIGRvbid0IGVuY29kZSBudWxsIHZhbHVlIHByb3BlcnRpZXNcblxuICAgIGlmICh0eXBlb2Yga2V5SW5kZXggPT09ICd1bmRlZmluZWQnKSB7XG4gICAgICBrZXlzLnB1c2goa2V5KVxuICAgICAga2V5SW5kZXggPSBrZXlzLmxlbmd0aCAtIDFcbiAgICAgIGtleWNhY2hlW2tleV0gPSBrZXlJbmRleFxuICAgIH1cbiAgICBwYmYud3JpdGVWYXJpbnQoa2V5SW5kZXgpXG5cbiAgICB2YXIgdHlwZSA9IHR5cGVvZiB2YWx1ZVxuICAgIGlmICh0eXBlICE9PSAnc3RyaW5nJyAmJiB0eXBlICE9PSAnYm9vbGVhbicgJiYgdHlwZSAhPT0gJ251bWJlcicpIHtcbiAgICAgIHZhbHVlID0gSlNPTi5zdHJpbmdpZnkodmFsdWUpXG4gICAgfVxuICAgIHZhciB2YWx1ZUtleSA9IHR5cGUgKyAnOicgKyB2YWx1ZVxuICAgIHZhciB2YWx1ZUluZGV4ID0gdmFsdWVjYWNoZVt2YWx1ZUtleV1cbiAgICBpZiAodHlwZW9mIHZhbHVlSW5kZXggPT09ICd1bmRlZmluZWQnKSB7XG4gICAgICB2YWx1ZXMucHVzaCh2YWx1ZSlcbiAgICAgIHZhbHVlSW5kZXggPSB2YWx1ZXMubGVuZ3RoIC0gMVxuICAgICAgdmFsdWVjYWNoZVt2YWx1ZUtleV0gPSB2YWx1ZUluZGV4XG4gICAgfVxuICAgIHBiZi53cml0ZVZhcmludCh2YWx1ZUluZGV4KVxuICB9XG59XG5cbmZ1bmN0aW9uIGNvbW1hbmQgKGNtZCwgbGVuZ3RoKSB7XG4gIHJldHVybiAobGVuZ3RoIDw8IDMpICsgKGNtZCAmIDB4Nylcbn1cblxuZnVuY3Rpb24gemlnemFnIChudW0pIHtcbiAgcmV0dXJuIChudW0gPDwgMSkgXiAobnVtID4+IDMxKVxufVxuXG5mdW5jdGlvbiB3cml0ZUdlb21ldHJ5IChmZWF0dXJlLCBwYmYpIHtcbiAgdmFyIGdlb21ldHJ5ID0gZmVhdHVyZS5sb2FkR2VvbWV0cnkoKVxuICB2YXIgdHlwZSA9IGZlYXR1cmUudHlwZVxuICB2YXIgeCA9IDBcbiAgdmFyIHkgPSAwXG4gIHZhciByaW5ncyA9IGdlb21ldHJ5Lmxlbmd0aFxuICBmb3IgKHZhciByID0gMDsgciA8IHJpbmdzOyByKyspIHtcbiAgICB2YXIgcmluZyA9IGdlb21ldHJ5W3JdXG4gICAgdmFyIGNvdW50ID0gMVxuICAgIGlmICh0eXBlID09PSAxKSB7XG4gICAgICBjb3VudCA9IHJpbmcubGVuZ3RoXG4gICAgfVxuICAgIHBiZi53cml0ZVZhcmludChjb21tYW5kKDEsIGNvdW50KSkgLy8gbW92ZXRvXG4gICAgLy8gZG8gbm90IHdyaXRlIHBvbHlnb24gY2xvc2luZyBwYXRoIGFzIGxpbmV0b1xuICAgIHZhciBsaW5lQ291bnQgPSB0eXBlID09PSAzID8gcmluZy5sZW5ndGggLSAxIDogcmluZy5sZW5ndGhcbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IGxpbmVDb3VudDsgaSsrKSB7XG4gICAgICBpZiAoaSA9PT0gMSAmJiB0eXBlICE9PSAxKSB7XG4gICAgICAgIHBiZi53cml0ZVZhcmludChjb21tYW5kKDIsIGxpbmVDb3VudCAtIDEpKSAvLyBsaW5ldG9cbiAgICAgIH1cbiAgICAgIHZhciBkeCA9IHJpbmdbaV0ueCAtIHhcbiAgICAgIHZhciBkeSA9IHJpbmdbaV0ueSAtIHlcbiAgICAgIHBiZi53cml0ZVZhcmludCh6aWd6YWcoZHgpKVxuICAgICAgcGJmLndyaXRlVmFyaW50KHppZ3phZyhkeSkpXG4gICAgICB4ICs9IGR4XG4gICAgICB5ICs9IGR5XG4gICAgfVxuICAgIGlmICh0eXBlID09PSAzKSB7XG4gICAgICBwYmYud3JpdGVWYXJpbnQoY29tbWFuZCg3LCAxKSkgLy8gY2xvc2VwYXRoXG4gICAgfVxuICB9XG59XG5cbmZ1bmN0aW9uIHdyaXRlVmFsdWUgKHZhbHVlLCBwYmYpIHtcbiAgdmFyIHR5cGUgPSB0eXBlb2YgdmFsdWVcbiAgaWYgKHR5cGUgPT09ICdzdHJpbmcnKSB7XG4gICAgcGJmLndyaXRlU3RyaW5nRmllbGQoMSwgdmFsdWUpXG4gIH0gZWxzZSBpZiAodHlwZSA9PT0gJ2Jvb2xlYW4nKSB7XG4gICAgcGJmLndyaXRlQm9vbGVhbkZpZWxkKDcsIHZhbHVlKVxuICB9IGVsc2UgaWYgKHR5cGUgPT09ICdudW1iZXInKSB7XG4gICAgaWYgKHZhbHVlICUgMSAhPT0gMCkge1xuICAgICAgcGJmLndyaXRlRG91YmxlRmllbGQoMywgdmFsdWUpXG4gICAgfSBlbHNlIGlmICh2YWx1ZSA8IDApIHtcbiAgICAgIHBiZi53cml0ZVNWYXJpbnRGaWVsZCg2LCB2YWx1ZSlcbiAgICB9IGVsc2Uge1xuICAgICAgcGJmLndyaXRlVmFyaW50RmllbGQoNSwgdmFsdWUpXG4gICAgfVxuICB9XG59XG4iLCJcbmV4cG9ydCBkZWZhdWx0IGZ1bmN0aW9uIHNvcnRLRChpZHMsIGNvb3Jkcywgbm9kZVNpemUsIGxlZnQsIHJpZ2h0LCBkZXB0aCkge1xuICAgIGlmIChyaWdodCAtIGxlZnQgPD0gbm9kZVNpemUpIHJldHVybjtcblxuICAgIGNvbnN0IG0gPSAobGVmdCArIHJpZ2h0KSA+PiAxO1xuXG4gICAgc2VsZWN0KGlkcywgY29vcmRzLCBtLCBsZWZ0LCByaWdodCwgZGVwdGggJSAyKTtcblxuICAgIHNvcnRLRChpZHMsIGNvb3Jkcywgbm9kZVNpemUsIGxlZnQsIG0gLSAxLCBkZXB0aCArIDEpO1xuICAgIHNvcnRLRChpZHMsIGNvb3Jkcywgbm9kZVNpemUsIG0gKyAxLCByaWdodCwgZGVwdGggKyAxKTtcbn1cblxuZnVuY3Rpb24gc2VsZWN0KGlkcywgY29vcmRzLCBrLCBsZWZ0LCByaWdodCwgaW5jKSB7XG5cbiAgICB3aGlsZSAocmlnaHQgPiBsZWZ0KSB7XG4gICAgICAgIGlmIChyaWdodCAtIGxlZnQgPiA2MDApIHtcbiAgICAgICAgICAgIGNvbnN0IG4gPSByaWdodCAtIGxlZnQgKyAxO1xuICAgICAgICAgICAgY29uc3QgbSA9IGsgLSBsZWZ0ICsgMTtcbiAgICAgICAgICAgIGNvbnN0IHogPSBNYXRoLmxvZyhuKTtcbiAgICAgICAgICAgIGNvbnN0IHMgPSAwLjUgKiBNYXRoLmV4cCgyICogeiAvIDMpO1xuICAgICAgICAgICAgY29uc3Qgc2QgPSAwLjUgKiBNYXRoLnNxcnQoeiAqIHMgKiAobiAtIHMpIC8gbikgKiAobSAtIG4gLyAyIDwgMCA/IC0xIDogMSk7XG4gICAgICAgICAgICBjb25zdCBuZXdMZWZ0ID0gTWF0aC5tYXgobGVmdCwgTWF0aC5mbG9vcihrIC0gbSAqIHMgLyBuICsgc2QpKTtcbiAgICAgICAgICAgIGNvbnN0IG5ld1JpZ2h0ID0gTWF0aC5taW4ocmlnaHQsIE1hdGguZmxvb3IoayArIChuIC0gbSkgKiBzIC8gbiArIHNkKSk7XG4gICAgICAgICAgICBzZWxlY3QoaWRzLCBjb29yZHMsIGssIG5ld0xlZnQsIG5ld1JpZ2h0LCBpbmMpO1xuICAgICAgICB9XG5cbiAgICAgICAgY29uc3QgdCA9IGNvb3Jkc1syICogayArIGluY107XG4gICAgICAgIGxldCBpID0gbGVmdDtcbiAgICAgICAgbGV0IGogPSByaWdodDtcblxuICAgICAgICBzd2FwSXRlbShpZHMsIGNvb3JkcywgbGVmdCwgayk7XG4gICAgICAgIGlmIChjb29yZHNbMiAqIHJpZ2h0ICsgaW5jXSA+IHQpIHN3YXBJdGVtKGlkcywgY29vcmRzLCBsZWZ0LCByaWdodCk7XG5cbiAgICAgICAgd2hpbGUgKGkgPCBqKSB7XG4gICAgICAgICAgICBzd2FwSXRlbShpZHMsIGNvb3JkcywgaSwgaik7XG4gICAgICAgICAgICBpKys7XG4gICAgICAgICAgICBqLS07XG4gICAgICAgICAgICB3aGlsZSAoY29vcmRzWzIgKiBpICsgaW5jXSA8IHQpIGkrKztcbiAgICAgICAgICAgIHdoaWxlIChjb29yZHNbMiAqIGogKyBpbmNdID4gdCkgai0tO1xuICAgICAgICB9XG5cbiAgICAgICAgaWYgKGNvb3Jkc1syICogbGVmdCArIGluY10gPT09IHQpIHN3YXBJdGVtKGlkcywgY29vcmRzLCBsZWZ0LCBqKTtcbiAgICAgICAgZWxzZSB7XG4gICAgICAgICAgICBqKys7XG4gICAgICAgICAgICBzd2FwSXRlbShpZHMsIGNvb3JkcywgaiwgcmlnaHQpO1xuICAgICAgICB9XG5cbiAgICAgICAgaWYgKGogPD0gaykgbGVmdCA9IGogKyAxO1xuICAgICAgICBpZiAoayA8PSBqKSByaWdodCA9IGogLSAxO1xuICAgIH1cbn1cblxuZnVuY3Rpb24gc3dhcEl0ZW0oaWRzLCBjb29yZHMsIGksIGopIHtcbiAgICBzd2FwKGlkcywgaSwgaik7XG4gICAgc3dhcChjb29yZHMsIDIgKiBpLCAyICogaik7XG4gICAgc3dhcChjb29yZHMsIDIgKiBpICsgMSwgMiAqIGogKyAxKTtcbn1cblxuZnVuY3Rpb24gc3dhcChhcnIsIGksIGopIHtcbiAgICBjb25zdCB0bXAgPSBhcnJbaV07XG4gICAgYXJyW2ldID0gYXJyW2pdO1xuICAgIGFycltqXSA9IHRtcDtcbn1cbiIsIlxuZXhwb3J0IGRlZmF1bHQgZnVuY3Rpb24gcmFuZ2UoaWRzLCBjb29yZHMsIG1pblgsIG1pblksIG1heFgsIG1heFksIG5vZGVTaXplKSB7XG4gICAgY29uc3Qgc3RhY2sgPSBbMCwgaWRzLmxlbmd0aCAtIDEsIDBdO1xuICAgIGNvbnN0IHJlc3VsdCA9IFtdO1xuICAgIGxldCB4LCB5O1xuXG4gICAgd2hpbGUgKHN0YWNrLmxlbmd0aCkge1xuICAgICAgICBjb25zdCBheGlzID0gc3RhY2sucG9wKCk7XG4gICAgICAgIGNvbnN0IHJpZ2h0ID0gc3RhY2sucG9wKCk7XG4gICAgICAgIGNvbnN0IGxlZnQgPSBzdGFjay5wb3AoKTtcblxuICAgICAgICBpZiAocmlnaHQgLSBsZWZ0IDw9IG5vZGVTaXplKSB7XG4gICAgICAgICAgICBmb3IgKGxldCBpID0gbGVmdDsgaSA8PSByaWdodDsgaSsrKSB7XG4gICAgICAgICAgICAgICAgeCA9IGNvb3Jkc1syICogaV07XG4gICAgICAgICAgICAgICAgeSA9IGNvb3Jkc1syICogaSArIDFdO1xuICAgICAgICAgICAgICAgIGlmICh4ID49IG1pblggJiYgeCA8PSBtYXhYICYmIHkgPj0gbWluWSAmJiB5IDw9IG1heFkpIHJlc3VsdC5wdXNoKGlkc1tpXSk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBjb250aW51ZTtcbiAgICAgICAgfVxuXG4gICAgICAgIGNvbnN0IG0gPSBNYXRoLmZsb29yKChsZWZ0ICsgcmlnaHQpIC8gMik7XG5cbiAgICAgICAgeCA9IGNvb3Jkc1syICogbV07XG4gICAgICAgIHkgPSBjb29yZHNbMiAqIG0gKyAxXTtcblxuICAgICAgICBpZiAoeCA+PSBtaW5YICYmIHggPD0gbWF4WCAmJiB5ID49IG1pblkgJiYgeSA8PSBtYXhZKSByZXN1bHQucHVzaChpZHNbbV0pO1xuXG4gICAgICAgIGNvbnN0IG5leHRBeGlzID0gKGF4aXMgKyAxKSAlIDI7XG5cbiAgICAgICAgaWYgKGF4aXMgPT09IDAgPyBtaW5YIDw9IHggOiBtaW5ZIDw9IHkpIHtcbiAgICAgICAgICAgIHN0YWNrLnB1c2gobGVmdCk7XG4gICAgICAgICAgICBzdGFjay5wdXNoKG0gLSAxKTtcbiAgICAgICAgICAgIHN0YWNrLnB1c2gobmV4dEF4aXMpO1xuICAgICAgICB9XG4gICAgICAgIGlmIChheGlzID09PSAwID8gbWF4WCA+PSB4IDogbWF4WSA+PSB5KSB7XG4gICAgICAgICAgICBzdGFjay5wdXNoKG0gKyAxKTtcbiAgICAgICAgICAgIHN0YWNrLnB1c2gocmlnaHQpO1xuICAgICAgICAgICAgc3RhY2sucHVzaChuZXh0QXhpcyk7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICByZXR1cm4gcmVzdWx0O1xufVxuIiwiXG5leHBvcnQgZGVmYXVsdCBmdW5jdGlvbiB3aXRoaW4oaWRzLCBjb29yZHMsIHF4LCBxeSwgciwgbm9kZVNpemUpIHtcbiAgICBjb25zdCBzdGFjayA9IFswLCBpZHMubGVuZ3RoIC0gMSwgMF07XG4gICAgY29uc3QgcmVzdWx0ID0gW107XG4gICAgY29uc3QgcjIgPSByICogcjtcblxuICAgIHdoaWxlIChzdGFjay5sZW5ndGgpIHtcbiAgICAgICAgY29uc3QgYXhpcyA9IHN0YWNrLnBvcCgpO1xuICAgICAgICBjb25zdCByaWdodCA9IHN0YWNrLnBvcCgpO1xuICAgICAgICBjb25zdCBsZWZ0ID0gc3RhY2sucG9wKCk7XG5cbiAgICAgICAgaWYgKHJpZ2h0IC0gbGVmdCA8PSBub2RlU2l6ZSkge1xuICAgICAgICAgICAgZm9yIChsZXQgaSA9IGxlZnQ7IGkgPD0gcmlnaHQ7IGkrKykge1xuICAgICAgICAgICAgICAgIGlmIChzcURpc3QoY29vcmRzWzIgKiBpXSwgY29vcmRzWzIgKiBpICsgMV0sIHF4LCBxeSkgPD0gcjIpIHJlc3VsdC5wdXNoKGlkc1tpXSk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBjb250aW51ZTtcbiAgICAgICAgfVxuXG4gICAgICAgIGNvbnN0IG0gPSBNYXRoLmZsb29yKChsZWZ0ICsgcmlnaHQpIC8gMik7XG5cbiAgICAgICAgY29uc3QgeCA9IGNvb3Jkc1syICogbV07XG4gICAgICAgIGNvbnN0IHkgPSBjb29yZHNbMiAqIG0gKyAxXTtcblxuICAgICAgICBpZiAoc3FEaXN0KHgsIHksIHF4LCBxeSkgPD0gcjIpIHJlc3VsdC5wdXNoKGlkc1ttXSk7XG5cbiAgICAgICAgY29uc3QgbmV4dEF4aXMgPSAoYXhpcyArIDEpICUgMjtcblxuICAgICAgICBpZiAoYXhpcyA9PT0gMCA/IHF4IC0gciA8PSB4IDogcXkgLSByIDw9IHkpIHtcbiAgICAgICAgICAgIHN0YWNrLnB1c2gobGVmdCk7XG4gICAgICAgICAgICBzdGFjay5wdXNoKG0gLSAxKTtcbiAgICAgICAgICAgIHN0YWNrLnB1c2gobmV4dEF4aXMpO1xuICAgICAgICB9XG4gICAgICAgIGlmIChheGlzID09PSAwID8gcXggKyByID49IHggOiBxeSArIHIgPj0geSkge1xuICAgICAgICAgICAgc3RhY2sucHVzaChtICsgMSk7XG4gICAgICAgICAgICBzdGFjay5wdXNoKHJpZ2h0KTtcbiAgICAgICAgICAgIHN0YWNrLnB1c2gobmV4dEF4aXMpO1xuICAgICAgICB9XG4gICAgfVxuXG4gICAgcmV0dXJuIHJlc3VsdDtcbn1cblxuZnVuY3Rpb24gc3FEaXN0KGF4LCBheSwgYngsIGJ5KSB7XG4gICAgY29uc3QgZHggPSBheCAtIGJ4O1xuICAgIGNvbnN0IGR5ID0gYXkgLSBieTtcbiAgICByZXR1cm4gZHggKiBkeCArIGR5ICogZHk7XG59XG4iLCJcbmltcG9ydCBzb3J0IGZyb20gJy4vc29ydCc7XG5pbXBvcnQgcmFuZ2UgZnJvbSAnLi9yYW5nZSc7XG5pbXBvcnQgd2l0aGluIGZyb20gJy4vd2l0aGluJztcblxuY29uc3QgZGVmYXVsdEdldFggPSBwID0+IHBbMF07XG5jb25zdCBkZWZhdWx0R2V0WSA9IHAgPT4gcFsxXTtcblxuZXhwb3J0IGRlZmF1bHQgY2xhc3MgS0RCdXNoIHtcbiAgICBjb25zdHJ1Y3Rvcihwb2ludHMsIGdldFggPSBkZWZhdWx0R2V0WCwgZ2V0WSA9IGRlZmF1bHRHZXRZLCBub2RlU2l6ZSA9IDY0LCBBcnJheVR5cGUgPSBGbG9hdDY0QXJyYXkpIHtcbiAgICAgICAgdGhpcy5ub2RlU2l6ZSA9IG5vZGVTaXplO1xuICAgICAgICB0aGlzLnBvaW50cyA9IHBvaW50cztcblxuICAgICAgICBjb25zdCBJbmRleEFycmF5VHlwZSA9IHBvaW50cy5sZW5ndGggPCA2NTUzNiA/IFVpbnQxNkFycmF5IDogVWludDMyQXJyYXk7XG5cbiAgICAgICAgY29uc3QgaWRzID0gdGhpcy5pZHMgPSBuZXcgSW5kZXhBcnJheVR5cGUocG9pbnRzLmxlbmd0aCk7XG4gICAgICAgIGNvbnN0IGNvb3JkcyA9IHRoaXMuY29vcmRzID0gbmV3IEFycmF5VHlwZShwb2ludHMubGVuZ3RoICogMik7XG5cbiAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBwb2ludHMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgICAgIGlkc1tpXSA9IGk7XG4gICAgICAgICAgICBjb29yZHNbMiAqIGldID0gZ2V0WChwb2ludHNbaV0pO1xuICAgICAgICAgICAgY29vcmRzWzIgKiBpICsgMV0gPSBnZXRZKHBvaW50c1tpXSk7XG4gICAgICAgIH1cblxuICAgICAgICBzb3J0KGlkcywgY29vcmRzLCBub2RlU2l6ZSwgMCwgaWRzLmxlbmd0aCAtIDEsIDApO1xuICAgIH1cblxuICAgIHJhbmdlKG1pblgsIG1pblksIG1heFgsIG1heFkpIHtcbiAgICAgICAgcmV0dXJuIHJhbmdlKHRoaXMuaWRzLCB0aGlzLmNvb3JkcywgbWluWCwgbWluWSwgbWF4WCwgbWF4WSwgdGhpcy5ub2RlU2l6ZSk7XG4gICAgfVxuXG4gICAgd2l0aGluKHgsIHksIHIpIHtcbiAgICAgICAgcmV0dXJuIHdpdGhpbih0aGlzLmlkcywgdGhpcy5jb29yZHMsIHgsIHksIHIsIHRoaXMubm9kZVNpemUpO1xuICAgIH1cbn1cbiIsIlxuaW1wb3J0IEtEQnVzaCBmcm9tICdrZGJ1c2gnO1xuXG5jb25zdCBkZWZhdWx0T3B0aW9ucyA9IHtcbiAgICBtaW5ab29tOiAwLCAgIC8vIG1pbiB6b29tIHRvIGdlbmVyYXRlIGNsdXN0ZXJzIG9uXG4gICAgbWF4Wm9vbTogMTYsICAvLyBtYXggem9vbSBsZXZlbCB0byBjbHVzdGVyIHRoZSBwb2ludHMgb25cbiAgICBtaW5Qb2ludHM6IDIsIC8vIG1pbmltdW0gcG9pbnRzIHRvIGZvcm0gYSBjbHVzdGVyXG4gICAgcmFkaXVzOiA0MCwgICAvLyBjbHVzdGVyIHJhZGl1cyBpbiBwaXhlbHNcbiAgICBleHRlbnQ6IDUxMiwgIC8vIHRpbGUgZXh0ZW50IChyYWRpdXMgaXMgY2FsY3VsYXRlZCByZWxhdGl2ZSB0byBpdClcbiAgICBub2RlU2l6ZTogNjQsIC8vIHNpemUgb2YgdGhlIEtELXRyZWUgbGVhZiBub2RlLCBhZmZlY3RzIHBlcmZvcm1hbmNlXG4gICAgbG9nOiBmYWxzZSwgICAvLyB3aGV0aGVyIHRvIGxvZyB0aW1pbmcgaW5mb1xuXG4gICAgLy8gd2hldGhlciB0byBnZW5lcmF0ZSBudW1lcmljIGlkcyBmb3IgaW5wdXQgZmVhdHVyZXMgKGluIHZlY3RvciB0aWxlcylcbiAgICBnZW5lcmF0ZUlkOiBmYWxzZSxcblxuICAgIC8vIGEgcmVkdWNlIGZ1bmN0aW9uIGZvciBjYWxjdWxhdGluZyBjdXN0b20gY2x1c3RlciBwcm9wZXJ0aWVzXG4gICAgcmVkdWNlOiBudWxsLCAvLyAoYWNjdW11bGF0ZWQsIHByb3BzKSA9PiB7IGFjY3VtdWxhdGVkLnN1bSArPSBwcm9wcy5zdW07IH1cblxuICAgIC8vIHByb3BlcnRpZXMgdG8gdXNlIGZvciBpbmRpdmlkdWFsIHBvaW50cyB3aGVuIHJ1bm5pbmcgdGhlIHJlZHVjZXJcbiAgICBtYXA6IHByb3BzID0+IHByb3BzIC8vIHByb3BzID0+ICh7c3VtOiBwcm9wcy5teV92YWx1ZX0pXG59O1xuXG5jb25zdCBmcm91bmQgPSBNYXRoLmZyb3VuZCB8fCAodG1wID0+ICgoeCkgPT4geyB0bXBbMF0gPSAreDsgcmV0dXJuIHRtcFswXTsgfSkpKG5ldyBGbG9hdDMyQXJyYXkoMSkpO1xuXG5leHBvcnQgZGVmYXVsdCBjbGFzcyBTdXBlcmNsdXN0ZXIge1xuICAgIGNvbnN0cnVjdG9yKG9wdGlvbnMpIHtcbiAgICAgICAgdGhpcy5vcHRpb25zID0gZXh0ZW5kKE9iamVjdC5jcmVhdGUoZGVmYXVsdE9wdGlvbnMpLCBvcHRpb25zKTtcbiAgICAgICAgdGhpcy50cmVlcyA9IG5ldyBBcnJheSh0aGlzLm9wdGlvbnMubWF4Wm9vbSArIDEpO1xuICAgIH1cblxuICAgIGxvYWQocG9pbnRzKSB7XG4gICAgICAgIGNvbnN0IHtsb2csIG1pblpvb20sIG1heFpvb20sIG5vZGVTaXplfSA9IHRoaXMub3B0aW9ucztcblxuICAgICAgICBpZiAobG9nKSBjb25zb2xlLnRpbWUoJ3RvdGFsIHRpbWUnKTtcblxuICAgICAgICBjb25zdCB0aW1lcklkID0gYHByZXBhcmUgJHsgIHBvaW50cy5sZW5ndGggIH0gcG9pbnRzYDtcbiAgICAgICAgaWYgKGxvZykgY29uc29sZS50aW1lKHRpbWVySWQpO1xuXG4gICAgICAgIHRoaXMucG9pbnRzID0gcG9pbnRzO1xuXG4gICAgICAgIC8vIGdlbmVyYXRlIGEgY2x1c3RlciBvYmplY3QgZm9yIGVhY2ggcG9pbnQgYW5kIGluZGV4IGlucHV0IHBvaW50cyBpbnRvIGEgS0QtdHJlZVxuICAgICAgICBsZXQgY2x1c3RlcnMgPSBbXTtcbiAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBwb2ludHMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgICAgIGlmICghcG9pbnRzW2ldLmdlb21ldHJ5KSBjb250aW51ZTtcbiAgICAgICAgICAgIGNsdXN0ZXJzLnB1c2goY3JlYXRlUG9pbnRDbHVzdGVyKHBvaW50c1tpXSwgaSkpO1xuICAgICAgICB9XG4gICAgICAgIHRoaXMudHJlZXNbbWF4Wm9vbSArIDFdID0gbmV3IEtEQnVzaChjbHVzdGVycywgZ2V0WCwgZ2V0WSwgbm9kZVNpemUsIEZsb2F0MzJBcnJheSk7XG5cbiAgICAgICAgaWYgKGxvZykgY29uc29sZS50aW1lRW5kKHRpbWVySWQpO1xuXG4gICAgICAgIC8vIGNsdXN0ZXIgcG9pbnRzIG9uIG1heCB6b29tLCB0aGVuIGNsdXN0ZXIgdGhlIHJlc3VsdHMgb24gcHJldmlvdXMgem9vbSwgZXRjLjtcbiAgICAgICAgLy8gcmVzdWx0cyBpbiBhIGNsdXN0ZXIgaGllcmFyY2h5IGFjcm9zcyB6b29tIGxldmVsc1xuICAgICAgICBmb3IgKGxldCB6ID0gbWF4Wm9vbTsgeiA+PSBtaW5ab29tOyB6LS0pIHtcbiAgICAgICAgICAgIGNvbnN0IG5vdyA9ICtEYXRlLm5vdygpO1xuXG4gICAgICAgICAgICAvLyBjcmVhdGUgYSBuZXcgc2V0IG9mIGNsdXN0ZXJzIGZvciB0aGUgem9vbSBhbmQgaW5kZXggdGhlbSB3aXRoIGEgS0QtdHJlZVxuICAgICAgICAgICAgY2x1c3RlcnMgPSB0aGlzLl9jbHVzdGVyKGNsdXN0ZXJzLCB6KTtcbiAgICAgICAgICAgIHRoaXMudHJlZXNbel0gPSBuZXcgS0RCdXNoKGNsdXN0ZXJzLCBnZXRYLCBnZXRZLCBub2RlU2l6ZSwgRmxvYXQzMkFycmF5KTtcblxuICAgICAgICAgICAgaWYgKGxvZykgY29uc29sZS5sb2coJ3olZDogJWQgY2x1c3RlcnMgaW4gJWRtcycsIHosIGNsdXN0ZXJzLmxlbmd0aCwgK0RhdGUubm93KCkgLSBub3cpO1xuICAgICAgICB9XG5cbiAgICAgICAgaWYgKGxvZykgY29uc29sZS50aW1lRW5kKCd0b3RhbCB0aW1lJyk7XG5cbiAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgfVxuXG4gICAgZ2V0Q2x1c3RlcnMoYmJveCwgem9vbSkge1xuICAgICAgICBsZXQgbWluTG5nID0gKChiYm94WzBdICsgMTgwKSAlIDM2MCArIDM2MCkgJSAzNjAgLSAxODA7XG4gICAgICAgIGNvbnN0IG1pbkxhdCA9IE1hdGgubWF4KC05MCwgTWF0aC5taW4oOTAsIGJib3hbMV0pKTtcbiAgICAgICAgbGV0IG1heExuZyA9IGJib3hbMl0gPT09IDE4MCA/IDE4MCA6ICgoYmJveFsyXSArIDE4MCkgJSAzNjAgKyAzNjApICUgMzYwIC0gMTgwO1xuICAgICAgICBjb25zdCBtYXhMYXQgPSBNYXRoLm1heCgtOTAsIE1hdGgubWluKDkwLCBiYm94WzNdKSk7XG5cbiAgICAgICAgaWYgKGJib3hbMl0gLSBiYm94WzBdID49IDM2MCkge1xuICAgICAgICAgICAgbWluTG5nID0gLTE4MDtcbiAgICAgICAgICAgIG1heExuZyA9IDE4MDtcbiAgICAgICAgfSBlbHNlIGlmIChtaW5MbmcgPiBtYXhMbmcpIHtcbiAgICAgICAgICAgIGNvbnN0IGVhc3Rlcm5IZW0gPSB0aGlzLmdldENsdXN0ZXJzKFttaW5MbmcsIG1pbkxhdCwgMTgwLCBtYXhMYXRdLCB6b29tKTtcbiAgICAgICAgICAgIGNvbnN0IHdlc3Rlcm5IZW0gPSB0aGlzLmdldENsdXN0ZXJzKFstMTgwLCBtaW5MYXQsIG1heExuZywgbWF4TGF0XSwgem9vbSk7XG4gICAgICAgICAgICByZXR1cm4gZWFzdGVybkhlbS5jb25jYXQod2VzdGVybkhlbSk7XG4gICAgICAgIH1cblxuICAgICAgICBjb25zdCB0cmVlID0gdGhpcy50cmVlc1t0aGlzLl9saW1pdFpvb20oem9vbSldO1xuICAgICAgICBjb25zdCBpZHMgPSB0cmVlLnJhbmdlKGxuZ1gobWluTG5nKSwgbGF0WShtYXhMYXQpLCBsbmdYKG1heExuZyksIGxhdFkobWluTGF0KSk7XG4gICAgICAgIGNvbnN0IGNsdXN0ZXJzID0gW107XG4gICAgICAgIGZvciAoY29uc3QgaWQgb2YgaWRzKSB7XG4gICAgICAgICAgICBjb25zdCBjID0gdHJlZS5wb2ludHNbaWRdO1xuICAgICAgICAgICAgY2x1c3RlcnMucHVzaChjLm51bVBvaW50cyA/IGdldENsdXN0ZXJKU09OKGMpIDogdGhpcy5wb2ludHNbYy5pbmRleF0pO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiBjbHVzdGVycztcbiAgICB9XG5cbiAgICBnZXRDaGlsZHJlbihjbHVzdGVySWQpIHtcbiAgICAgICAgY29uc3Qgb3JpZ2luSWQgPSB0aGlzLl9nZXRPcmlnaW5JZChjbHVzdGVySWQpO1xuICAgICAgICBjb25zdCBvcmlnaW5ab29tID0gdGhpcy5fZ2V0T3JpZ2luWm9vbShjbHVzdGVySWQpO1xuICAgICAgICBjb25zdCBlcnJvck1zZyA9ICdObyBjbHVzdGVyIHdpdGggdGhlIHNwZWNpZmllZCBpZC4nO1xuXG4gICAgICAgIGNvbnN0IGluZGV4ID0gdGhpcy50cmVlc1tvcmlnaW5ab29tXTtcbiAgICAgICAgaWYgKCFpbmRleCkgdGhyb3cgbmV3IEVycm9yKGVycm9yTXNnKTtcblxuICAgICAgICBjb25zdCBvcmlnaW4gPSBpbmRleC5wb2ludHNbb3JpZ2luSWRdO1xuICAgICAgICBpZiAoIW9yaWdpbikgdGhyb3cgbmV3IEVycm9yKGVycm9yTXNnKTtcblxuICAgICAgICBjb25zdCByID0gdGhpcy5vcHRpb25zLnJhZGl1cyAvICh0aGlzLm9wdGlvbnMuZXh0ZW50ICogTWF0aC5wb3coMiwgb3JpZ2luWm9vbSAtIDEpKTtcbiAgICAgICAgY29uc3QgaWRzID0gaW5kZXgud2l0aGluKG9yaWdpbi54LCBvcmlnaW4ueSwgcik7XG4gICAgICAgIGNvbnN0IGNoaWxkcmVuID0gW107XG4gICAgICAgIGZvciAoY29uc3QgaWQgb2YgaWRzKSB7XG4gICAgICAgICAgICBjb25zdCBjID0gaW5kZXgucG9pbnRzW2lkXTtcbiAgICAgICAgICAgIGlmIChjLnBhcmVudElkID09PSBjbHVzdGVySWQpIHtcbiAgICAgICAgICAgICAgICBjaGlsZHJlbi5wdXNoKGMubnVtUG9pbnRzID8gZ2V0Q2x1c3RlckpTT04oYykgOiB0aGlzLnBvaW50c1tjLmluZGV4XSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICBpZiAoY2hpbGRyZW4ubGVuZ3RoID09PSAwKSB0aHJvdyBuZXcgRXJyb3IoZXJyb3JNc2cpO1xuXG4gICAgICAgIHJldHVybiBjaGlsZHJlbjtcbiAgICB9XG5cbiAgICBnZXRMZWF2ZXMoY2x1c3RlcklkLCBsaW1pdCwgb2Zmc2V0KSB7XG4gICAgICAgIGxpbWl0ID0gbGltaXQgfHwgMTA7XG4gICAgICAgIG9mZnNldCA9IG9mZnNldCB8fCAwO1xuXG4gICAgICAgIGNvbnN0IGxlYXZlcyA9IFtdO1xuICAgICAgICB0aGlzLl9hcHBlbmRMZWF2ZXMobGVhdmVzLCBjbHVzdGVySWQsIGxpbWl0LCBvZmZzZXQsIDApO1xuXG4gICAgICAgIHJldHVybiBsZWF2ZXM7XG4gICAgfVxuXG4gICAgZ2V0VGlsZSh6LCB4LCB5KSB7XG4gICAgICAgIGNvbnN0IHRyZWUgPSB0aGlzLnRyZWVzW3RoaXMuX2xpbWl0Wm9vbSh6KV07XG4gICAgICAgIGNvbnN0IHoyID0gTWF0aC5wb3coMiwgeik7XG4gICAgICAgIGNvbnN0IHtleHRlbnQsIHJhZGl1c30gPSB0aGlzLm9wdGlvbnM7XG4gICAgICAgIGNvbnN0IHAgPSByYWRpdXMgLyBleHRlbnQ7XG4gICAgICAgIGNvbnN0IHRvcCA9ICh5IC0gcCkgLyB6MjtcbiAgICAgICAgY29uc3QgYm90dG9tID0gKHkgKyAxICsgcCkgLyB6MjtcblxuICAgICAgICBjb25zdCB0aWxlID0ge1xuICAgICAgICAgICAgZmVhdHVyZXM6IFtdXG4gICAgICAgIH07XG5cbiAgICAgICAgdGhpcy5fYWRkVGlsZUZlYXR1cmVzKFxuICAgICAgICAgICAgdHJlZS5yYW5nZSgoeCAtIHApIC8gejIsIHRvcCwgKHggKyAxICsgcCkgLyB6MiwgYm90dG9tKSxcbiAgICAgICAgICAgIHRyZWUucG9pbnRzLCB4LCB5LCB6MiwgdGlsZSk7XG5cbiAgICAgICAgaWYgKHggPT09IDApIHtcbiAgICAgICAgICAgIHRoaXMuX2FkZFRpbGVGZWF0dXJlcyhcbiAgICAgICAgICAgICAgICB0cmVlLnJhbmdlKDEgLSBwIC8gejIsIHRvcCwgMSwgYm90dG9tKSxcbiAgICAgICAgICAgICAgICB0cmVlLnBvaW50cywgejIsIHksIHoyLCB0aWxlKTtcbiAgICAgICAgfVxuICAgICAgICBpZiAoeCA9PT0gejIgLSAxKSB7XG4gICAgICAgICAgICB0aGlzLl9hZGRUaWxlRmVhdHVyZXMoXG4gICAgICAgICAgICAgICAgdHJlZS5yYW5nZSgwLCB0b3AsIHAgLyB6MiwgYm90dG9tKSxcbiAgICAgICAgICAgICAgICB0cmVlLnBvaW50cywgLTEsIHksIHoyLCB0aWxlKTtcbiAgICAgICAgfVxuXG4gICAgICAgIHJldHVybiB0aWxlLmZlYXR1cmVzLmxlbmd0aCA/IHRpbGUgOiBudWxsO1xuICAgIH1cblxuICAgIGdldENsdXN0ZXJFeHBhbnNpb25ab29tKGNsdXN0ZXJJZCkge1xuICAgICAgICBsZXQgZXhwYW5zaW9uWm9vbSA9IHRoaXMuX2dldE9yaWdpblpvb20oY2x1c3RlcklkKSAtIDE7XG4gICAgICAgIHdoaWxlIChleHBhbnNpb25ab29tIDw9IHRoaXMub3B0aW9ucy5tYXhab29tKSB7XG4gICAgICAgICAgICBjb25zdCBjaGlsZHJlbiA9IHRoaXMuZ2V0Q2hpbGRyZW4oY2x1c3RlcklkKTtcbiAgICAgICAgICAgIGV4cGFuc2lvblpvb20rKztcbiAgICAgICAgICAgIGlmIChjaGlsZHJlbi5sZW5ndGggIT09IDEpIGJyZWFrO1xuICAgICAgICAgICAgY2x1c3RlcklkID0gY2hpbGRyZW5bMF0ucHJvcGVydGllcy5jbHVzdGVyX2lkO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiBleHBhbnNpb25ab29tO1xuICAgIH1cblxuICAgIF9hcHBlbmRMZWF2ZXMocmVzdWx0LCBjbHVzdGVySWQsIGxpbWl0LCBvZmZzZXQsIHNraXBwZWQpIHtcbiAgICAgICAgY29uc3QgY2hpbGRyZW4gPSB0aGlzLmdldENoaWxkcmVuKGNsdXN0ZXJJZCk7XG5cbiAgICAgICAgZm9yIChjb25zdCBjaGlsZCBvZiBjaGlsZHJlbikge1xuICAgICAgICAgICAgY29uc3QgcHJvcHMgPSBjaGlsZC5wcm9wZXJ0aWVzO1xuXG4gICAgICAgICAgICBpZiAocHJvcHMgJiYgcHJvcHMuY2x1c3Rlcikge1xuICAgICAgICAgICAgICAgIGlmIChza2lwcGVkICsgcHJvcHMucG9pbnRfY291bnQgPD0gb2Zmc2V0KSB7XG4gICAgICAgICAgICAgICAgICAgIC8vIHNraXAgdGhlIHdob2xlIGNsdXN0ZXJcbiAgICAgICAgICAgICAgICAgICAgc2tpcHBlZCArPSBwcm9wcy5wb2ludF9jb3VudDtcbiAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICAvLyBlbnRlciB0aGUgY2x1c3RlclxuICAgICAgICAgICAgICAgICAgICBza2lwcGVkID0gdGhpcy5fYXBwZW5kTGVhdmVzKHJlc3VsdCwgcHJvcHMuY2x1c3Rlcl9pZCwgbGltaXQsIG9mZnNldCwgc2tpcHBlZCk7XG4gICAgICAgICAgICAgICAgICAgIC8vIGV4aXQgdGhlIGNsdXN0ZXJcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9IGVsc2UgaWYgKHNraXBwZWQgPCBvZmZzZXQpIHtcbiAgICAgICAgICAgICAgICAvLyBza2lwIGEgc2luZ2xlIHBvaW50XG4gICAgICAgICAgICAgICAgc2tpcHBlZCsrO1xuICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAvLyBhZGQgYSBzaW5nbGUgcG9pbnRcbiAgICAgICAgICAgICAgICByZXN1bHQucHVzaChjaGlsZCk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBpZiAocmVzdWx0Lmxlbmd0aCA9PT0gbGltaXQpIGJyZWFrO1xuICAgICAgICB9XG5cbiAgICAgICAgcmV0dXJuIHNraXBwZWQ7XG4gICAgfVxuXG4gICAgX2FkZFRpbGVGZWF0dXJlcyhpZHMsIHBvaW50cywgeCwgeSwgejIsIHRpbGUpIHtcbiAgICAgICAgZm9yIChjb25zdCBpIG9mIGlkcykge1xuICAgICAgICAgICAgY29uc3QgYyA9IHBvaW50c1tpXTtcbiAgICAgICAgICAgIGNvbnN0IGlzQ2x1c3RlciA9IGMubnVtUG9pbnRzO1xuXG4gICAgICAgICAgICBsZXQgdGFncywgcHgsIHB5O1xuICAgICAgICAgICAgaWYgKGlzQ2x1c3Rlcikge1xuICAgICAgICAgICAgICAgIHRhZ3MgPSBnZXRDbHVzdGVyUHJvcGVydGllcyhjKTtcbiAgICAgICAgICAgICAgICBweCA9IGMueDtcbiAgICAgICAgICAgICAgICBweSA9IGMueTtcbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgY29uc3QgcCA9IHRoaXMucG9pbnRzW2MuaW5kZXhdO1xuICAgICAgICAgICAgICAgIHRhZ3MgPSBwLnByb3BlcnRpZXM7XG4gICAgICAgICAgICAgICAgcHggPSBsbmdYKHAuZ2VvbWV0cnkuY29vcmRpbmF0ZXNbMF0pO1xuICAgICAgICAgICAgICAgIHB5ID0gbGF0WShwLmdlb21ldHJ5LmNvb3JkaW5hdGVzWzFdKTtcbiAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgY29uc3QgZiA9IHtcbiAgICAgICAgICAgICAgICB0eXBlOiAxLFxuICAgICAgICAgICAgICAgIGdlb21ldHJ5OiBbW1xuICAgICAgICAgICAgICAgICAgICBNYXRoLnJvdW5kKHRoaXMub3B0aW9ucy5leHRlbnQgKiAocHggKiB6MiAtIHgpKSxcbiAgICAgICAgICAgICAgICAgICAgTWF0aC5yb3VuZCh0aGlzLm9wdGlvbnMuZXh0ZW50ICogKHB5ICogejIgLSB5KSlcbiAgICAgICAgICAgICAgICBdXSxcbiAgICAgICAgICAgICAgICB0YWdzXG4gICAgICAgICAgICB9O1xuXG4gICAgICAgICAgICAvLyBhc3NpZ24gaWRcbiAgICAgICAgICAgIGxldCBpZDtcbiAgICAgICAgICAgIGlmIChpc0NsdXN0ZXIpIHtcbiAgICAgICAgICAgICAgICBpZCA9IGMuaWQ7XG4gICAgICAgICAgICB9IGVsc2UgaWYgKHRoaXMub3B0aW9ucy5nZW5lcmF0ZUlkKSB7XG4gICAgICAgICAgICAgICAgLy8gb3B0aW9uYWxseSBnZW5lcmF0ZSBpZFxuICAgICAgICAgICAgICAgIGlkID0gYy5pbmRleDtcbiAgICAgICAgICAgIH0gZWxzZSBpZiAodGhpcy5wb2ludHNbYy5pbmRleF0uaWQpIHtcbiAgICAgICAgICAgICAgICAvLyBrZWVwIGlkIGlmIGFscmVhZHkgYXNzaWduZWRcbiAgICAgICAgICAgICAgICBpZCA9IHRoaXMucG9pbnRzW2MuaW5kZXhdLmlkO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICBpZiAoaWQgIT09IHVuZGVmaW5lZCkgZi5pZCA9IGlkO1xuXG4gICAgICAgICAgICB0aWxlLmZlYXR1cmVzLnB1c2goZik7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICBfbGltaXRab29tKHopIHtcbiAgICAgICAgcmV0dXJuIE1hdGgubWF4KHRoaXMub3B0aW9ucy5taW5ab29tLCBNYXRoLm1pbihNYXRoLmZsb29yKCt6KSwgdGhpcy5vcHRpb25zLm1heFpvb20gKyAxKSk7XG4gICAgfVxuXG4gICAgX2NsdXN0ZXIocG9pbnRzLCB6b29tKSB7XG4gICAgICAgIGNvbnN0IGNsdXN0ZXJzID0gW107XG4gICAgICAgIGNvbnN0IHtyYWRpdXMsIGV4dGVudCwgcmVkdWNlLCBtaW5Qb2ludHN9ID0gdGhpcy5vcHRpb25zO1xuICAgICAgICBjb25zdCByID0gcmFkaXVzIC8gKGV4dGVudCAqIE1hdGgucG93KDIsIHpvb20pKTtcblxuICAgICAgICAvLyBsb29wIHRocm91Z2ggZWFjaCBwb2ludFxuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IHBvaW50cy5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICAgY29uc3QgcCA9IHBvaW50c1tpXTtcbiAgICAgICAgICAgIC8vIGlmIHdlJ3ZlIGFscmVhZHkgdmlzaXRlZCB0aGUgcG9pbnQgYXQgdGhpcyB6b29tIGxldmVsLCBza2lwIGl0XG4gICAgICAgICAgICBpZiAocC56b29tIDw9IHpvb20pIGNvbnRpbnVlO1xuICAgICAgICAgICAgcC56b29tID0gem9vbTtcblxuICAgICAgICAgICAgLy8gZmluZCBhbGwgbmVhcmJ5IHBvaW50c1xuICAgICAgICAgICAgY29uc3QgdHJlZSA9IHRoaXMudHJlZXNbem9vbSArIDFdO1xuICAgICAgICAgICAgY29uc3QgbmVpZ2hib3JJZHMgPSB0cmVlLndpdGhpbihwLngsIHAueSwgcik7XG5cbiAgICAgICAgICAgIGNvbnN0IG51bVBvaW50c09yaWdpbiA9IHAubnVtUG9pbnRzIHx8IDE7XG4gICAgICAgICAgICBsZXQgbnVtUG9pbnRzID0gbnVtUG9pbnRzT3JpZ2luO1xuXG4gICAgICAgICAgICAvLyBjb3VudCB0aGUgbnVtYmVyIG9mIHBvaW50cyBpbiBhIHBvdGVudGlhbCBjbHVzdGVyXG4gICAgICAgICAgICBmb3IgKGNvbnN0IG5laWdoYm9ySWQgb2YgbmVpZ2hib3JJZHMpIHtcbiAgICAgICAgICAgICAgICBjb25zdCBiID0gdHJlZS5wb2ludHNbbmVpZ2hib3JJZF07XG4gICAgICAgICAgICAgICAgLy8gZmlsdGVyIG91dCBuZWlnaGJvcnMgdGhhdCBhcmUgYWxyZWFkeSBwcm9jZXNzZWRcbiAgICAgICAgICAgICAgICBpZiAoYi56b29tID4gem9vbSkgbnVtUG9pbnRzICs9IGIubnVtUG9pbnRzIHx8IDE7XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgICAgIC8vIGlmIHRoZXJlIHdlcmUgbmVpZ2hib3JzIHRvIG1lcmdlLCBhbmQgdGhlcmUgYXJlIGVub3VnaCBwb2ludHMgdG8gZm9ybSBhIGNsdXN0ZXJcbiAgICAgICAgICAgIGlmIChudW1Qb2ludHMgPiBudW1Qb2ludHNPcmlnaW4gJiYgbnVtUG9pbnRzID49IG1pblBvaW50cykge1xuICAgICAgICAgICAgICAgIGxldCB3eCA9IHAueCAqIG51bVBvaW50c09yaWdpbjtcbiAgICAgICAgICAgICAgICBsZXQgd3kgPSBwLnkgKiBudW1Qb2ludHNPcmlnaW47XG5cbiAgICAgICAgICAgICAgICBsZXQgY2x1c3RlclByb3BlcnRpZXMgPSByZWR1Y2UgJiYgbnVtUG9pbnRzT3JpZ2luID4gMSA/IHRoaXMuX21hcChwLCB0cnVlKSA6IG51bGw7XG5cbiAgICAgICAgICAgICAgICAvLyBlbmNvZGUgYm90aCB6b29tIGFuZCBwb2ludCBpbmRleCBvbiB3aGljaCB0aGUgY2x1c3RlciBvcmlnaW5hdGVkIC0tIG9mZnNldCBieSB0b3RhbCBsZW5ndGggb2YgZmVhdHVyZXNcbiAgICAgICAgICAgICAgICBjb25zdCBpZCA9IChpIDw8IDUpICsgKHpvb20gKyAxKSArIHRoaXMucG9pbnRzLmxlbmd0aDtcblxuICAgICAgICAgICAgICAgIGZvciAoY29uc3QgbmVpZ2hib3JJZCBvZiBuZWlnaGJvcklkcykge1xuICAgICAgICAgICAgICAgICAgICBjb25zdCBiID0gdHJlZS5wb2ludHNbbmVpZ2hib3JJZF07XG5cbiAgICAgICAgICAgICAgICAgICAgaWYgKGIuem9vbSA8PSB6b29tKSBjb250aW51ZTtcbiAgICAgICAgICAgICAgICAgICAgYi56b29tID0gem9vbTsgLy8gc2F2ZSB0aGUgem9vbSAoc28gaXQgZG9lc24ndCBnZXQgcHJvY2Vzc2VkIHR3aWNlKVxuXG4gICAgICAgICAgICAgICAgICAgIGNvbnN0IG51bVBvaW50czIgPSBiLm51bVBvaW50cyB8fCAxO1xuICAgICAgICAgICAgICAgICAgICB3eCArPSBiLnggKiBudW1Qb2ludHMyOyAvLyBhY2N1bXVsYXRlIGNvb3JkaW5hdGVzIGZvciBjYWxjdWxhdGluZyB3ZWlnaHRlZCBjZW50ZXJcbiAgICAgICAgICAgICAgICAgICAgd3kgKz0gYi55ICogbnVtUG9pbnRzMjtcblxuICAgICAgICAgICAgICAgICAgICBiLnBhcmVudElkID0gaWQ7XG5cbiAgICAgICAgICAgICAgICAgICAgaWYgKHJlZHVjZSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgaWYgKCFjbHVzdGVyUHJvcGVydGllcykgY2x1c3RlclByb3BlcnRpZXMgPSB0aGlzLl9tYXAocCwgdHJ1ZSk7XG4gICAgICAgICAgICAgICAgICAgICAgICByZWR1Y2UoY2x1c3RlclByb3BlcnRpZXMsIHRoaXMuX21hcChiKSk7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICBwLnBhcmVudElkID0gaWQ7XG4gICAgICAgICAgICAgICAgY2x1c3RlcnMucHVzaChjcmVhdGVDbHVzdGVyKHd4IC8gbnVtUG9pbnRzLCB3eSAvIG51bVBvaW50cywgaWQsIG51bVBvaW50cywgY2x1c3RlclByb3BlcnRpZXMpKTtcblxuICAgICAgICAgICAgfSBlbHNlIHsgLy8gbGVmdCBwb2ludHMgYXMgdW5jbHVzdGVyZWRcbiAgICAgICAgICAgICAgICBjbHVzdGVycy5wdXNoKHApO1xuXG4gICAgICAgICAgICAgICAgaWYgKG51bVBvaW50cyA+IDEpIHtcbiAgICAgICAgICAgICAgICAgICAgZm9yIChjb25zdCBuZWlnaGJvcklkIG9mIG5laWdoYm9ySWRzKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBjb25zdCBiID0gdHJlZS5wb2ludHNbbmVpZ2hib3JJZF07XG4gICAgICAgICAgICAgICAgICAgICAgICBpZiAoYi56b29tIDw9IHpvb20pIGNvbnRpbnVlO1xuICAgICAgICAgICAgICAgICAgICAgICAgYi56b29tID0gem9vbTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGNsdXN0ZXJzLnB1c2goYik7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICByZXR1cm4gY2x1c3RlcnM7XG4gICAgfVxuXG4gICAgLy8gZ2V0IGluZGV4IG9mIHRoZSBwb2ludCBmcm9tIHdoaWNoIHRoZSBjbHVzdGVyIG9yaWdpbmF0ZWRcbiAgICBfZ2V0T3JpZ2luSWQoY2x1c3RlcklkKSB7XG4gICAgICAgIHJldHVybiAoY2x1c3RlcklkIC0gdGhpcy5wb2ludHMubGVuZ3RoKSA+PiA1O1xuICAgIH1cblxuICAgIC8vIGdldCB6b29tIG9mIHRoZSBwb2ludCBmcm9tIHdoaWNoIHRoZSBjbHVzdGVyIG9yaWdpbmF0ZWRcbiAgICBfZ2V0T3JpZ2luWm9vbShjbHVzdGVySWQpIHtcbiAgICAgICAgcmV0dXJuIChjbHVzdGVySWQgLSB0aGlzLnBvaW50cy5sZW5ndGgpICUgMzI7XG4gICAgfVxuXG4gICAgX21hcChwb2ludCwgY2xvbmUpIHtcbiAgICAgICAgaWYgKHBvaW50Lm51bVBvaW50cykge1xuICAgICAgICAgICAgcmV0dXJuIGNsb25lID8gZXh0ZW5kKHt9LCBwb2ludC5wcm9wZXJ0aWVzKSA6IHBvaW50LnByb3BlcnRpZXM7XG4gICAgICAgIH1cbiAgICAgICAgY29uc3Qgb3JpZ2luYWwgPSB0aGlzLnBvaW50c1twb2ludC5pbmRleF0ucHJvcGVydGllcztcbiAgICAgICAgY29uc3QgcmVzdWx0ID0gdGhpcy5vcHRpb25zLm1hcChvcmlnaW5hbCk7XG4gICAgICAgIHJldHVybiBjbG9uZSAmJiByZXN1bHQgPT09IG9yaWdpbmFsID8gZXh0ZW5kKHt9LCByZXN1bHQpIDogcmVzdWx0O1xuICAgIH1cbn1cblxuZnVuY3Rpb24gY3JlYXRlQ2x1c3Rlcih4LCB5LCBpZCwgbnVtUG9pbnRzLCBwcm9wZXJ0aWVzKSB7XG4gICAgcmV0dXJuIHtcbiAgICAgICAgeDogZnJvdW5kKHgpLCAvLyB3ZWlnaHRlZCBjbHVzdGVyIGNlbnRlcjsgcm91bmQgZm9yIGNvbnNpc3RlbmN5IHdpdGggRmxvYXQzMkFycmF5IGluZGV4XG4gICAgICAgIHk6IGZyb3VuZCh5KSxcbiAgICAgICAgem9vbTogSW5maW5pdHksIC8vIHRoZSBsYXN0IHpvb20gdGhlIGNsdXN0ZXIgd2FzIHByb2Nlc3NlZCBhdFxuICAgICAgICBpZCwgLy8gZW5jb2RlcyBpbmRleCBvZiB0aGUgZmlyc3QgY2hpbGQgb2YgdGhlIGNsdXN0ZXIgYW5kIGl0cyB6b29tIGxldmVsXG4gICAgICAgIHBhcmVudElkOiAtMSwgLy8gcGFyZW50IGNsdXN0ZXIgaWRcbiAgICAgICAgbnVtUG9pbnRzLFxuICAgICAgICBwcm9wZXJ0aWVzXG4gICAgfTtcbn1cblxuZnVuY3Rpb24gY3JlYXRlUG9pbnRDbHVzdGVyKHAsIGlkKSB7XG4gICAgY29uc3QgW3gsIHldID0gcC5nZW9tZXRyeS5jb29yZGluYXRlcztcbiAgICByZXR1cm4ge1xuICAgICAgICB4OiBmcm91bmQobG5nWCh4KSksIC8vIHByb2plY3RlZCBwb2ludCBjb29yZGluYXRlc1xuICAgICAgICB5OiBmcm91bmQobGF0WSh5KSksXG4gICAgICAgIHpvb206IEluZmluaXR5LCAvLyB0aGUgbGFzdCB6b29tIHRoZSBwb2ludCB3YXMgcHJvY2Vzc2VkIGF0XG4gICAgICAgIGluZGV4OiBpZCwgLy8gaW5kZXggb2YgdGhlIHNvdXJjZSBmZWF0dXJlIGluIHRoZSBvcmlnaW5hbCBpbnB1dCBhcnJheSxcbiAgICAgICAgcGFyZW50SWQ6IC0xIC8vIHBhcmVudCBjbHVzdGVyIGlkXG4gICAgfTtcbn1cblxuZnVuY3Rpb24gZ2V0Q2x1c3RlckpTT04oY2x1c3Rlcikge1xuICAgIHJldHVybiB7XG4gICAgICAgIHR5cGU6ICdGZWF0dXJlJyxcbiAgICAgICAgaWQ6IGNsdXN0ZXIuaWQsXG4gICAgICAgIHByb3BlcnRpZXM6IGdldENsdXN0ZXJQcm9wZXJ0aWVzKGNsdXN0ZXIpLFxuICAgICAgICBnZW9tZXRyeToge1xuICAgICAgICAgICAgdHlwZTogJ1BvaW50JyxcbiAgICAgICAgICAgIGNvb3JkaW5hdGVzOiBbeExuZyhjbHVzdGVyLngpLCB5TGF0KGNsdXN0ZXIueSldXG4gICAgICAgIH1cbiAgICB9O1xufVxuXG5mdW5jdGlvbiBnZXRDbHVzdGVyUHJvcGVydGllcyhjbHVzdGVyKSB7XG4gICAgY29uc3QgY291bnQgPSBjbHVzdGVyLm51bVBvaW50cztcbiAgICBjb25zdCBhYmJyZXYgPVxuICAgICAgICBjb3VudCA+PSAxMDAwMCA/IGAke01hdGgucm91bmQoY291bnQgLyAxMDAwKSAgfWtgIDpcbiAgICAgICAgY291bnQgPj0gMTAwMCA/IGAke01hdGgucm91bmQoY291bnQgLyAxMDApIC8gMTAgIH1rYCA6IGNvdW50O1xuICAgIHJldHVybiBleHRlbmQoZXh0ZW5kKHt9LCBjbHVzdGVyLnByb3BlcnRpZXMpLCB7XG4gICAgICAgIGNsdXN0ZXI6IHRydWUsXG4gICAgICAgIGNsdXN0ZXJfaWQ6IGNsdXN0ZXIuaWQsXG4gICAgICAgIHBvaW50X2NvdW50OiBjb3VudCxcbiAgICAgICAgcG9pbnRfY291bnRfYWJicmV2aWF0ZWQ6IGFiYnJldlxuICAgIH0pO1xufVxuXG4vLyBsb25naXR1ZGUvbGF0aXR1ZGUgdG8gc3BoZXJpY2FsIG1lcmNhdG9yIGluIFswLi4xXSByYW5nZVxuZnVuY3Rpb24gbG5nWChsbmcpIHtcbiAgICByZXR1cm4gbG5nIC8gMzYwICsgMC41O1xufVxuZnVuY3Rpb24gbGF0WShsYXQpIHtcbiAgICBjb25zdCBzaW4gPSBNYXRoLnNpbihsYXQgKiBNYXRoLlBJIC8gMTgwKTtcbiAgICBjb25zdCB5ID0gKDAuNSAtIDAuMjUgKiBNYXRoLmxvZygoMSArIHNpbikgLyAoMSAtIHNpbikpIC8gTWF0aC5QSSk7XG4gICAgcmV0dXJuIHkgPCAwID8gMCA6IHkgPiAxID8gMSA6IHk7XG59XG5cbi8vIHNwaGVyaWNhbCBtZXJjYXRvciB0byBsb25naXR1ZGUvbGF0aXR1ZGVcbmZ1bmN0aW9uIHhMbmcoeCkge1xuICAgIHJldHVybiAoeCAtIDAuNSkgKiAzNjA7XG59XG5mdW5jdGlvbiB5TGF0KHkpIHtcbiAgICBjb25zdCB5MiA9ICgxODAgLSB5ICogMzYwKSAqIE1hdGguUEkgLyAxODA7XG4gICAgcmV0dXJuIDM2MCAqIE1hdGguYXRhbihNYXRoLmV4cCh5MikpIC8gTWF0aC5QSSAtIDkwO1xufVxuXG5mdW5jdGlvbiBleHRlbmQoZGVzdCwgc3JjKSB7XG4gICAgZm9yIChjb25zdCBpZCBpbiBzcmMpIGRlc3RbaWRdID0gc3JjW2lkXTtcbiAgICByZXR1cm4gZGVzdDtcbn1cblxuZnVuY3Rpb24gZ2V0WChwKSB7XG4gICAgcmV0dXJuIHAueDtcbn1cbmZ1bmN0aW9uIGdldFkocCkge1xuICAgIHJldHVybiBwLnk7XG59XG4iLCJcbi8vIGNhbGN1bGF0ZSBzaW1wbGlmaWNhdGlvbiBkYXRhIHVzaW5nIG9wdGltaXplZCBEb3VnbGFzLVBldWNrZXIgYWxnb3JpdGhtXG5cbmV4cG9ydCBkZWZhdWx0IGZ1bmN0aW9uIHNpbXBsaWZ5KGNvb3JkcywgZmlyc3QsIGxhc3QsIHNxVG9sZXJhbmNlKSB7XG4gICAgdmFyIG1heFNxRGlzdCA9IHNxVG9sZXJhbmNlO1xuICAgIHZhciBtaWQgPSAobGFzdCAtIGZpcnN0KSA+PiAxO1xuICAgIHZhciBtaW5Qb3NUb01pZCA9IGxhc3QgLSBmaXJzdDtcbiAgICB2YXIgaW5kZXg7XG5cbiAgICB2YXIgYXggPSBjb29yZHNbZmlyc3RdO1xuICAgIHZhciBheSA9IGNvb3Jkc1tmaXJzdCArIDFdO1xuICAgIHZhciBieCA9IGNvb3Jkc1tsYXN0XTtcbiAgICB2YXIgYnkgPSBjb29yZHNbbGFzdCArIDFdO1xuXG4gICAgZm9yICh2YXIgaSA9IGZpcnN0ICsgMzsgaSA8IGxhc3Q7IGkgKz0gMykge1xuICAgICAgICB2YXIgZCA9IGdldFNxU2VnRGlzdChjb29yZHNbaV0sIGNvb3Jkc1tpICsgMV0sIGF4LCBheSwgYngsIGJ5KTtcblxuICAgICAgICBpZiAoZCA+IG1heFNxRGlzdCkge1xuICAgICAgICAgICAgaW5kZXggPSBpO1xuICAgICAgICAgICAgbWF4U3FEaXN0ID0gZDtcblxuICAgICAgICB9IGVsc2UgaWYgKGQgPT09IG1heFNxRGlzdCkge1xuICAgICAgICAgICAgLy8gYSB3b3JrYXJvdW5kIHRvIGVuc3VyZSB3ZSBjaG9vc2UgYSBwaXZvdCBjbG9zZSB0byB0aGUgbWlkZGxlIG9mIHRoZSBsaXN0LFxuICAgICAgICAgICAgLy8gcmVkdWNpbmcgcmVjdXJzaW9uIGRlcHRoLCBmb3IgY2VydGFpbiBkZWdlbmVyYXRlIGlucHV0c1xuICAgICAgICAgICAgLy8gaHR0cHM6Ly9naXRodWIuY29tL21hcGJveC9nZW9qc29uLXZ0L2lzc3Vlcy8xMDRcbiAgICAgICAgICAgIHZhciBwb3NUb01pZCA9IE1hdGguYWJzKGkgLSBtaWQpO1xuICAgICAgICAgICAgaWYgKHBvc1RvTWlkIDwgbWluUG9zVG9NaWQpIHtcbiAgICAgICAgICAgICAgICBpbmRleCA9IGk7XG4gICAgICAgICAgICAgICAgbWluUG9zVG9NaWQgPSBwb3NUb01pZDtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH1cblxuICAgIGlmIChtYXhTcURpc3QgPiBzcVRvbGVyYW5jZSkge1xuICAgICAgICBpZiAoaW5kZXggLSBmaXJzdCA+IDMpIHNpbXBsaWZ5KGNvb3JkcywgZmlyc3QsIGluZGV4LCBzcVRvbGVyYW5jZSk7XG4gICAgICAgIGNvb3Jkc1tpbmRleCArIDJdID0gbWF4U3FEaXN0O1xuICAgICAgICBpZiAobGFzdCAtIGluZGV4ID4gMykgc2ltcGxpZnkoY29vcmRzLCBpbmRleCwgbGFzdCwgc3FUb2xlcmFuY2UpO1xuICAgIH1cbn1cblxuLy8gc3F1YXJlIGRpc3RhbmNlIGZyb20gYSBwb2ludCB0byBhIHNlZ21lbnRcbmZ1bmN0aW9uIGdldFNxU2VnRGlzdChweCwgcHksIHgsIHksIGJ4LCBieSkge1xuXG4gICAgdmFyIGR4ID0gYnggLSB4O1xuICAgIHZhciBkeSA9IGJ5IC0geTtcblxuICAgIGlmIChkeCAhPT0gMCB8fCBkeSAhPT0gMCkge1xuXG4gICAgICAgIHZhciB0ID0gKChweCAtIHgpICogZHggKyAocHkgLSB5KSAqIGR5KSAvIChkeCAqIGR4ICsgZHkgKiBkeSk7XG5cbiAgICAgICAgaWYgKHQgPiAxKSB7XG4gICAgICAgICAgICB4ID0gYng7XG4gICAgICAgICAgICB5ID0gYnk7XG5cbiAgICAgICAgfSBlbHNlIGlmICh0ID4gMCkge1xuICAgICAgICAgICAgeCArPSBkeCAqIHQ7XG4gICAgICAgICAgICB5ICs9IGR5ICogdDtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIGR4ID0gcHggLSB4O1xuICAgIGR5ID0gcHkgLSB5O1xuXG4gICAgcmV0dXJuIGR4ICogZHggKyBkeSAqIGR5O1xufVxuIiwiXG5leHBvcnQgZGVmYXVsdCBmdW5jdGlvbiBjcmVhdGVGZWF0dXJlKGlkLCB0eXBlLCBnZW9tLCB0YWdzKSB7XG4gICAgdmFyIGZlYXR1cmUgPSB7XG4gICAgICAgIGlkOiB0eXBlb2YgaWQgPT09ICd1bmRlZmluZWQnID8gbnVsbCA6IGlkLFxuICAgICAgICB0eXBlOiB0eXBlLFxuICAgICAgICBnZW9tZXRyeTogZ2VvbSxcbiAgICAgICAgdGFnczogdGFncyxcbiAgICAgICAgbWluWDogSW5maW5pdHksXG4gICAgICAgIG1pblk6IEluZmluaXR5LFxuICAgICAgICBtYXhYOiAtSW5maW5pdHksXG4gICAgICAgIG1heFk6IC1JbmZpbml0eVxuICAgIH07XG4gICAgY2FsY0JCb3goZmVhdHVyZSk7XG4gICAgcmV0dXJuIGZlYXR1cmU7XG59XG5cbmZ1bmN0aW9uIGNhbGNCQm94KGZlYXR1cmUpIHtcbiAgICB2YXIgZ2VvbSA9IGZlYXR1cmUuZ2VvbWV0cnk7XG4gICAgdmFyIHR5cGUgPSBmZWF0dXJlLnR5cGU7XG5cbiAgICBpZiAodHlwZSA9PT0gJ1BvaW50JyB8fCB0eXBlID09PSAnTXVsdGlQb2ludCcgfHwgdHlwZSA9PT0gJ0xpbmVTdHJpbmcnKSB7XG4gICAgICAgIGNhbGNMaW5lQkJveChmZWF0dXJlLCBnZW9tKTtcblxuICAgIH0gZWxzZSBpZiAodHlwZSA9PT0gJ1BvbHlnb24nIHx8IHR5cGUgPT09ICdNdWx0aUxpbmVTdHJpbmcnKSB7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgZ2VvbS5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICAgY2FsY0xpbmVCQm94KGZlYXR1cmUsIGdlb21baV0pO1xuICAgICAgICB9XG5cbiAgICB9IGVsc2UgaWYgKHR5cGUgPT09ICdNdWx0aVBvbHlnb24nKSB7XG4gICAgICAgIGZvciAoaSA9IDA7IGkgPCBnZW9tLmxlbmd0aDsgaSsrKSB7XG4gICAgICAgICAgICBmb3IgKHZhciBqID0gMDsgaiA8IGdlb21baV0ubGVuZ3RoOyBqKyspIHtcbiAgICAgICAgICAgICAgICBjYWxjTGluZUJCb3goZmVhdHVyZSwgZ2VvbVtpXVtqXSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICB9XG59XG5cbmZ1bmN0aW9uIGNhbGNMaW5lQkJveChmZWF0dXJlLCBnZW9tKSB7XG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBnZW9tLmxlbmd0aDsgaSArPSAzKSB7XG4gICAgICAgIGZlYXR1cmUubWluWCA9IE1hdGgubWluKGZlYXR1cmUubWluWCwgZ2VvbVtpXSk7XG4gICAgICAgIGZlYXR1cmUubWluWSA9IE1hdGgubWluKGZlYXR1cmUubWluWSwgZ2VvbVtpICsgMV0pO1xuICAgICAgICBmZWF0dXJlLm1heFggPSBNYXRoLm1heChmZWF0dXJlLm1heFgsIGdlb21baV0pO1xuICAgICAgICBmZWF0dXJlLm1heFkgPSBNYXRoLm1heChmZWF0dXJlLm1heFksIGdlb21baSArIDFdKTtcbiAgICB9XG59XG4iLCJcbmltcG9ydCBzaW1wbGlmeSBmcm9tICcuL3NpbXBsaWZ5JztcbmltcG9ydCBjcmVhdGVGZWF0dXJlIGZyb20gJy4vZmVhdHVyZSc7XG5cbi8vIGNvbnZlcnRzIEdlb0pTT04gZmVhdHVyZSBpbnRvIGFuIGludGVybWVkaWF0ZSBwcm9qZWN0ZWQgSlNPTiB2ZWN0b3IgZm9ybWF0IHdpdGggc2ltcGxpZmljYXRpb24gZGF0YVxuXG5leHBvcnQgZGVmYXVsdCBmdW5jdGlvbiBjb252ZXJ0KGRhdGEsIG9wdGlvbnMpIHtcbiAgICB2YXIgZmVhdHVyZXMgPSBbXTtcbiAgICBpZiAoZGF0YS50eXBlID09PSAnRmVhdHVyZUNvbGxlY3Rpb24nKSB7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgZGF0YS5mZWF0dXJlcy5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICAgY29udmVydEZlYXR1cmUoZmVhdHVyZXMsIGRhdGEuZmVhdHVyZXNbaV0sIG9wdGlvbnMsIGkpO1xuICAgICAgICB9XG5cbiAgICB9IGVsc2UgaWYgKGRhdGEudHlwZSA9PT0gJ0ZlYXR1cmUnKSB7XG4gICAgICAgIGNvbnZlcnRGZWF0dXJlKGZlYXR1cmVzLCBkYXRhLCBvcHRpb25zKTtcblxuICAgIH0gZWxzZSB7XG4gICAgICAgIC8vIHNpbmdsZSBnZW9tZXRyeSBvciBhIGdlb21ldHJ5IGNvbGxlY3Rpb25cbiAgICAgICAgY29udmVydEZlYXR1cmUoZmVhdHVyZXMsIHtnZW9tZXRyeTogZGF0YX0sIG9wdGlvbnMpO1xuICAgIH1cblxuICAgIHJldHVybiBmZWF0dXJlcztcbn1cblxuZnVuY3Rpb24gY29udmVydEZlYXR1cmUoZmVhdHVyZXMsIGdlb2pzb24sIG9wdGlvbnMsIGluZGV4KSB7XG4gICAgaWYgKCFnZW9qc29uLmdlb21ldHJ5KSByZXR1cm47XG5cbiAgICB2YXIgY29vcmRzID0gZ2VvanNvbi5nZW9tZXRyeS5jb29yZGluYXRlcztcbiAgICB2YXIgdHlwZSA9IGdlb2pzb24uZ2VvbWV0cnkudHlwZTtcbiAgICB2YXIgdG9sZXJhbmNlID0gTWF0aC5wb3cob3B0aW9ucy50b2xlcmFuY2UgLyAoKDEgPDwgb3B0aW9ucy5tYXhab29tKSAqIG9wdGlvbnMuZXh0ZW50KSwgMik7XG4gICAgdmFyIGdlb21ldHJ5ID0gW107XG4gICAgdmFyIGlkID0gZ2VvanNvbi5pZDtcbiAgICBpZiAob3B0aW9ucy5wcm9tb3RlSWQpIHtcbiAgICAgICAgaWQgPSBnZW9qc29uLnByb3BlcnRpZXNbb3B0aW9ucy5wcm9tb3RlSWRdO1xuICAgIH0gZWxzZSBpZiAob3B0aW9ucy5nZW5lcmF0ZUlkKSB7XG4gICAgICAgIGlkID0gaW5kZXggfHwgMDtcbiAgICB9XG4gICAgaWYgKHR5cGUgPT09ICdQb2ludCcpIHtcbiAgICAgICAgY29udmVydFBvaW50KGNvb3JkcywgZ2VvbWV0cnkpO1xuXG4gICAgfSBlbHNlIGlmICh0eXBlID09PSAnTXVsdGlQb2ludCcpIHtcbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBjb29yZHMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgICAgIGNvbnZlcnRQb2ludChjb29yZHNbaV0sIGdlb21ldHJ5KTtcbiAgICAgICAgfVxuXG4gICAgfSBlbHNlIGlmICh0eXBlID09PSAnTGluZVN0cmluZycpIHtcbiAgICAgICAgY29udmVydExpbmUoY29vcmRzLCBnZW9tZXRyeSwgdG9sZXJhbmNlLCBmYWxzZSk7XG5cbiAgICB9IGVsc2UgaWYgKHR5cGUgPT09ICdNdWx0aUxpbmVTdHJpbmcnKSB7XG4gICAgICAgIGlmIChvcHRpb25zLmxpbmVNZXRyaWNzKSB7XG4gICAgICAgICAgICAvLyBleHBsb2RlIGludG8gbGluZXN0cmluZ3MgdG8gYmUgYWJsZSB0byB0cmFjayBtZXRyaWNzXG4gICAgICAgICAgICBmb3IgKGkgPSAwOyBpIDwgY29vcmRzLmxlbmd0aDsgaSsrKSB7XG4gICAgICAgICAgICAgICAgZ2VvbWV0cnkgPSBbXTtcbiAgICAgICAgICAgICAgICBjb252ZXJ0TGluZShjb29yZHNbaV0sIGdlb21ldHJ5LCB0b2xlcmFuY2UsIGZhbHNlKTtcbiAgICAgICAgICAgICAgICBmZWF0dXJlcy5wdXNoKGNyZWF0ZUZlYXR1cmUoaWQsICdMaW5lU3RyaW5nJywgZ2VvbWV0cnksIGdlb2pzb24ucHJvcGVydGllcykpO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcmV0dXJuO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgY29udmVydExpbmVzKGNvb3JkcywgZ2VvbWV0cnksIHRvbGVyYW5jZSwgZmFsc2UpO1xuICAgICAgICB9XG5cbiAgICB9IGVsc2UgaWYgKHR5cGUgPT09ICdQb2x5Z29uJykge1xuICAgICAgICBjb252ZXJ0TGluZXMoY29vcmRzLCBnZW9tZXRyeSwgdG9sZXJhbmNlLCB0cnVlKTtcblxuICAgIH0gZWxzZSBpZiAodHlwZSA9PT0gJ011bHRpUG9seWdvbicpIHtcbiAgICAgICAgZm9yIChpID0gMDsgaSA8IGNvb3Jkcy5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICAgdmFyIHBvbHlnb24gPSBbXTtcbiAgICAgICAgICAgIGNvbnZlcnRMaW5lcyhjb29yZHNbaV0sIHBvbHlnb24sIHRvbGVyYW5jZSwgdHJ1ZSk7XG4gICAgICAgICAgICBnZW9tZXRyeS5wdXNoKHBvbHlnb24pO1xuICAgICAgICB9XG4gICAgfSBlbHNlIGlmICh0eXBlID09PSAnR2VvbWV0cnlDb2xsZWN0aW9uJykge1xuICAgICAgICBmb3IgKGkgPSAwOyBpIDwgZ2VvanNvbi5nZW9tZXRyeS5nZW9tZXRyaWVzLmxlbmd0aDsgaSsrKSB7XG4gICAgICAgICAgICBjb252ZXJ0RmVhdHVyZShmZWF0dXJlcywge1xuICAgICAgICAgICAgICAgIGlkOiBpZCxcbiAgICAgICAgICAgICAgICBnZW9tZXRyeTogZ2VvanNvbi5nZW9tZXRyeS5nZW9tZXRyaWVzW2ldLFxuICAgICAgICAgICAgICAgIHByb3BlcnRpZXM6IGdlb2pzb24ucHJvcGVydGllc1xuICAgICAgICAgICAgfSwgb3B0aW9ucywgaW5kZXgpO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybjtcbiAgICB9IGVsc2Uge1xuICAgICAgICB0aHJvdyBuZXcgRXJyb3IoJ0lucHV0IGRhdGEgaXMgbm90IGEgdmFsaWQgR2VvSlNPTiBvYmplY3QuJyk7XG4gICAgfVxuXG4gICAgZmVhdHVyZXMucHVzaChjcmVhdGVGZWF0dXJlKGlkLCB0eXBlLCBnZW9tZXRyeSwgZ2VvanNvbi5wcm9wZXJ0aWVzKSk7XG59XG5cbmZ1bmN0aW9uIGNvbnZlcnRQb2ludChjb29yZHMsIG91dCkge1xuICAgIG91dC5wdXNoKHByb2plY3RYKGNvb3Jkc1swXSkpO1xuICAgIG91dC5wdXNoKHByb2plY3RZKGNvb3Jkc1sxXSkpO1xuICAgIG91dC5wdXNoKDApO1xufVxuXG5mdW5jdGlvbiBjb252ZXJ0TGluZShyaW5nLCBvdXQsIHRvbGVyYW5jZSwgaXNQb2x5Z29uKSB7XG4gICAgdmFyIHgwLCB5MDtcbiAgICB2YXIgc2l6ZSA9IDA7XG5cbiAgICBmb3IgKHZhciBqID0gMDsgaiA8IHJpbmcubGVuZ3RoOyBqKyspIHtcbiAgICAgICAgdmFyIHggPSBwcm9qZWN0WChyaW5nW2pdWzBdKTtcbiAgICAgICAgdmFyIHkgPSBwcm9qZWN0WShyaW5nW2pdWzFdKTtcblxuICAgICAgICBvdXQucHVzaCh4KTtcbiAgICAgICAgb3V0LnB1c2goeSk7XG4gICAgICAgIG91dC5wdXNoKDApO1xuXG4gICAgICAgIGlmIChqID4gMCkge1xuICAgICAgICAgICAgaWYgKGlzUG9seWdvbikge1xuICAgICAgICAgICAgICAgIHNpemUgKz0gKHgwICogeSAtIHggKiB5MCkgLyAyOyAvLyBhcmVhXG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIHNpemUgKz0gTWF0aC5zcXJ0KE1hdGgucG93KHggLSB4MCwgMikgKyBNYXRoLnBvdyh5IC0geTAsIDIpKTsgLy8gbGVuZ3RoXG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgeDAgPSB4O1xuICAgICAgICB5MCA9IHk7XG4gICAgfVxuXG4gICAgdmFyIGxhc3QgPSBvdXQubGVuZ3RoIC0gMztcbiAgICBvdXRbMl0gPSAxO1xuICAgIHNpbXBsaWZ5KG91dCwgMCwgbGFzdCwgdG9sZXJhbmNlKTtcbiAgICBvdXRbbGFzdCArIDJdID0gMTtcblxuICAgIG91dC5zaXplID0gTWF0aC5hYnMoc2l6ZSk7XG4gICAgb3V0LnN0YXJ0ID0gMDtcbiAgICBvdXQuZW5kID0gb3V0LnNpemU7XG59XG5cbmZ1bmN0aW9uIGNvbnZlcnRMaW5lcyhyaW5ncywgb3V0LCB0b2xlcmFuY2UsIGlzUG9seWdvbikge1xuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgcmluZ3MubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgdmFyIGdlb20gPSBbXTtcbiAgICAgICAgY29udmVydExpbmUocmluZ3NbaV0sIGdlb20sIHRvbGVyYW5jZSwgaXNQb2x5Z29uKTtcbiAgICAgICAgb3V0LnB1c2goZ2VvbSk7XG4gICAgfVxufVxuXG5mdW5jdGlvbiBwcm9qZWN0WCh4KSB7XG4gICAgcmV0dXJuIHggLyAzNjAgKyAwLjU7XG59XG5cbmZ1bmN0aW9uIHByb2plY3RZKHkpIHtcbiAgICB2YXIgc2luID0gTWF0aC5zaW4oeSAqIE1hdGguUEkgLyAxODApO1xuICAgIHZhciB5MiA9IDAuNSAtIDAuMjUgKiBNYXRoLmxvZygoMSArIHNpbikgLyAoMSAtIHNpbikpIC8gTWF0aC5QSTtcbiAgICByZXR1cm4geTIgPCAwID8gMCA6IHkyID4gMSA/IDEgOiB5Mjtcbn1cbiIsIlxuaW1wb3J0IGNyZWF0ZUZlYXR1cmUgZnJvbSAnLi9mZWF0dXJlJztcblxuLyogY2xpcCBmZWF0dXJlcyBiZXR3ZWVuIHR3byBheGlzLXBhcmFsbGVsIGxpbmVzOlxuICogICAgIHwgICAgICAgIHxcbiAqICBfX198X19fICAgICB8ICAgICAvXG4gKiAvICAgfCAgIFxcX19fX3xfX19fL1xuICogICAgIHwgICAgICAgIHxcbiAqL1xuXG5leHBvcnQgZGVmYXVsdCBmdW5jdGlvbiBjbGlwKGZlYXR1cmVzLCBzY2FsZSwgazEsIGsyLCBheGlzLCBtaW5BbGwsIG1heEFsbCwgb3B0aW9ucykge1xuXG4gICAgazEgLz0gc2NhbGU7XG4gICAgazIgLz0gc2NhbGU7XG5cbiAgICBpZiAobWluQWxsID49IGsxICYmIG1heEFsbCA8IGsyKSByZXR1cm4gZmVhdHVyZXM7IC8vIHRyaXZpYWwgYWNjZXB0XG4gICAgZWxzZSBpZiAobWF4QWxsIDwgazEgfHwgbWluQWxsID49IGsyKSByZXR1cm4gbnVsbDsgLy8gdHJpdmlhbCByZWplY3RcblxuICAgIHZhciBjbGlwcGVkID0gW107XG5cbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IGZlYXR1cmVzLmxlbmd0aDsgaSsrKSB7XG5cbiAgICAgICAgdmFyIGZlYXR1cmUgPSBmZWF0dXJlc1tpXTtcbiAgICAgICAgdmFyIGdlb21ldHJ5ID0gZmVhdHVyZS5nZW9tZXRyeTtcbiAgICAgICAgdmFyIHR5cGUgPSBmZWF0dXJlLnR5cGU7XG5cbiAgICAgICAgdmFyIG1pbiA9IGF4aXMgPT09IDAgPyBmZWF0dXJlLm1pblggOiBmZWF0dXJlLm1pblk7XG4gICAgICAgIHZhciBtYXggPSBheGlzID09PSAwID8gZmVhdHVyZS5tYXhYIDogZmVhdHVyZS5tYXhZO1xuXG4gICAgICAgIGlmIChtaW4gPj0gazEgJiYgbWF4IDwgazIpIHsgLy8gdHJpdmlhbCBhY2NlcHRcbiAgICAgICAgICAgIGNsaXBwZWQucHVzaChmZWF0dXJlKTtcbiAgICAgICAgICAgIGNvbnRpbnVlO1xuICAgICAgICB9IGVsc2UgaWYgKG1heCA8IGsxIHx8IG1pbiA+PSBrMikgeyAvLyB0cml2aWFsIHJlamVjdFxuICAgICAgICAgICAgY29udGludWU7XG4gICAgICAgIH1cblxuICAgICAgICB2YXIgbmV3R2VvbWV0cnkgPSBbXTtcblxuICAgICAgICBpZiAodHlwZSA9PT0gJ1BvaW50JyB8fCB0eXBlID09PSAnTXVsdGlQb2ludCcpIHtcbiAgICAgICAgICAgIGNsaXBQb2ludHMoZ2VvbWV0cnksIG5ld0dlb21ldHJ5LCBrMSwgazIsIGF4aXMpO1xuXG4gICAgICAgIH0gZWxzZSBpZiAodHlwZSA9PT0gJ0xpbmVTdHJpbmcnKSB7XG4gICAgICAgICAgICBjbGlwTGluZShnZW9tZXRyeSwgbmV3R2VvbWV0cnksIGsxLCBrMiwgYXhpcywgZmFsc2UsIG9wdGlvbnMubGluZU1ldHJpY3MpO1xuXG4gICAgICAgIH0gZWxzZSBpZiAodHlwZSA9PT0gJ011bHRpTGluZVN0cmluZycpIHtcbiAgICAgICAgICAgIGNsaXBMaW5lcyhnZW9tZXRyeSwgbmV3R2VvbWV0cnksIGsxLCBrMiwgYXhpcywgZmFsc2UpO1xuXG4gICAgICAgIH0gZWxzZSBpZiAodHlwZSA9PT0gJ1BvbHlnb24nKSB7XG4gICAgICAgICAgICBjbGlwTGluZXMoZ2VvbWV0cnksIG5ld0dlb21ldHJ5LCBrMSwgazIsIGF4aXMsIHRydWUpO1xuXG4gICAgICAgIH0gZWxzZSBpZiAodHlwZSA9PT0gJ011bHRpUG9seWdvbicpIHtcbiAgICAgICAgICAgIGZvciAodmFyIGogPSAwOyBqIDwgZ2VvbWV0cnkubGVuZ3RoOyBqKyspIHtcbiAgICAgICAgICAgICAgICB2YXIgcG9seWdvbiA9IFtdO1xuICAgICAgICAgICAgICAgIGNsaXBMaW5lcyhnZW9tZXRyeVtqXSwgcG9seWdvbiwgazEsIGsyLCBheGlzLCB0cnVlKTtcbiAgICAgICAgICAgICAgICBpZiAocG9seWdvbi5sZW5ndGgpIHtcbiAgICAgICAgICAgICAgICAgICAgbmV3R2VvbWV0cnkucHVzaChwb2x5Z29uKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICBpZiAobmV3R2VvbWV0cnkubGVuZ3RoKSB7XG4gICAgICAgICAgICBpZiAob3B0aW9ucy5saW5lTWV0cmljcyAmJiB0eXBlID09PSAnTGluZVN0cmluZycpIHtcbiAgICAgICAgICAgICAgICBmb3IgKGogPSAwOyBqIDwgbmV3R2VvbWV0cnkubGVuZ3RoOyBqKyspIHtcbiAgICAgICAgICAgICAgICAgICAgY2xpcHBlZC5wdXNoKGNyZWF0ZUZlYXR1cmUoZmVhdHVyZS5pZCwgdHlwZSwgbmV3R2VvbWV0cnlbal0sIGZlYXR1cmUudGFncykpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBjb250aW51ZTtcbiAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgaWYgKHR5cGUgPT09ICdMaW5lU3RyaW5nJyB8fCB0eXBlID09PSAnTXVsdGlMaW5lU3RyaW5nJykge1xuICAgICAgICAgICAgICAgIGlmIChuZXdHZW9tZXRyeS5sZW5ndGggPT09IDEpIHtcbiAgICAgICAgICAgICAgICAgICAgdHlwZSA9ICdMaW5lU3RyaW5nJztcbiAgICAgICAgICAgICAgICAgICAgbmV3R2VvbWV0cnkgPSBuZXdHZW9tZXRyeVswXTtcbiAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICB0eXBlID0gJ011bHRpTGluZVN0cmluZyc7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICAgICAgaWYgKHR5cGUgPT09ICdQb2ludCcgfHwgdHlwZSA9PT0gJ011bHRpUG9pbnQnKSB7XG4gICAgICAgICAgICAgICAgdHlwZSA9IG5ld0dlb21ldHJ5Lmxlbmd0aCA9PT0gMyA/ICdQb2ludCcgOiAnTXVsdGlQb2ludCc7XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgICAgIGNsaXBwZWQucHVzaChjcmVhdGVGZWF0dXJlKGZlYXR1cmUuaWQsIHR5cGUsIG5ld0dlb21ldHJ5LCBmZWF0dXJlLnRhZ3MpKTtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIHJldHVybiBjbGlwcGVkLmxlbmd0aCA/IGNsaXBwZWQgOiBudWxsO1xufVxuXG5mdW5jdGlvbiBjbGlwUG9pbnRzKGdlb20sIG5ld0dlb20sIGsxLCBrMiwgYXhpcykge1xuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgZ2VvbS5sZW5ndGg7IGkgKz0gMykge1xuICAgICAgICB2YXIgYSA9IGdlb21baSArIGF4aXNdO1xuXG4gICAgICAgIGlmIChhID49IGsxICYmIGEgPD0gazIpIHtcbiAgICAgICAgICAgIG5ld0dlb20ucHVzaChnZW9tW2ldKTtcbiAgICAgICAgICAgIG5ld0dlb20ucHVzaChnZW9tW2kgKyAxXSk7XG4gICAgICAgICAgICBuZXdHZW9tLnB1c2goZ2VvbVtpICsgMl0pO1xuICAgICAgICB9XG4gICAgfVxufVxuXG5mdW5jdGlvbiBjbGlwTGluZShnZW9tLCBuZXdHZW9tLCBrMSwgazIsIGF4aXMsIGlzUG9seWdvbiwgdHJhY2tNZXRyaWNzKSB7XG5cbiAgICB2YXIgc2xpY2UgPSBuZXdTbGljZShnZW9tKTtcbiAgICB2YXIgaW50ZXJzZWN0ID0gYXhpcyA9PT0gMCA/IGludGVyc2VjdFggOiBpbnRlcnNlY3RZO1xuICAgIHZhciBsZW4gPSBnZW9tLnN0YXJ0O1xuICAgIHZhciBzZWdMZW4sIHQ7XG5cbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IGdlb20ubGVuZ3RoIC0gMzsgaSArPSAzKSB7XG4gICAgICAgIHZhciBheCA9IGdlb21baV07XG4gICAgICAgIHZhciBheSA9IGdlb21baSArIDFdO1xuICAgICAgICB2YXIgYXogPSBnZW9tW2kgKyAyXTtcbiAgICAgICAgdmFyIGJ4ID0gZ2VvbVtpICsgM107XG4gICAgICAgIHZhciBieSA9IGdlb21baSArIDRdO1xuICAgICAgICB2YXIgYSA9IGF4aXMgPT09IDAgPyBheCA6IGF5O1xuICAgICAgICB2YXIgYiA9IGF4aXMgPT09IDAgPyBieCA6IGJ5O1xuICAgICAgICB2YXIgZXhpdGVkID0gZmFsc2U7XG5cbiAgICAgICAgaWYgKHRyYWNrTWV0cmljcykgc2VnTGVuID0gTWF0aC5zcXJ0KE1hdGgucG93KGF4IC0gYngsIDIpICsgTWF0aC5wb3coYXkgLSBieSwgMikpO1xuXG4gICAgICAgIGlmIChhIDwgazEpIHtcbiAgICAgICAgICAgIC8vIC0tLXwtLT4gIHwgKGxpbmUgZW50ZXJzIHRoZSBjbGlwIHJlZ2lvbiBmcm9tIHRoZSBsZWZ0KVxuICAgICAgICAgICAgaWYgKGIgPiBrMSkge1xuICAgICAgICAgICAgICAgIHQgPSBpbnRlcnNlY3Qoc2xpY2UsIGF4LCBheSwgYngsIGJ5LCBrMSk7XG4gICAgICAgICAgICAgICAgaWYgKHRyYWNrTWV0cmljcykgc2xpY2Uuc3RhcnQgPSBsZW4gKyBzZWdMZW4gKiB0O1xuICAgICAgICAgICAgfVxuICAgICAgICB9IGVsc2UgaWYgKGEgPiBrMikge1xuICAgICAgICAgICAgLy8gfCAgPC0tfC0tLSAobGluZSBlbnRlcnMgdGhlIGNsaXAgcmVnaW9uIGZyb20gdGhlIHJpZ2h0KVxuICAgICAgICAgICAgaWYgKGIgPCBrMikge1xuICAgICAgICAgICAgICAgIHQgPSBpbnRlcnNlY3Qoc2xpY2UsIGF4LCBheSwgYngsIGJ5LCBrMik7XG4gICAgICAgICAgICAgICAgaWYgKHRyYWNrTWV0cmljcykgc2xpY2Uuc3RhcnQgPSBsZW4gKyBzZWdMZW4gKiB0O1xuICAgICAgICAgICAgfVxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgYWRkUG9pbnQoc2xpY2UsIGF4LCBheSwgYXopO1xuICAgICAgICB9XG4gICAgICAgIGlmIChiIDwgazEgJiYgYSA+PSBrMSkge1xuICAgICAgICAgICAgLy8gPC0tfC0tLSAgfCBvciA8LS18LS0tLS18LS0tIChsaW5lIGV4aXRzIHRoZSBjbGlwIHJlZ2lvbiBvbiB0aGUgbGVmdClcbiAgICAgICAgICAgIHQgPSBpbnRlcnNlY3Qoc2xpY2UsIGF4LCBheSwgYngsIGJ5LCBrMSk7XG4gICAgICAgICAgICBleGl0ZWQgPSB0cnVlO1xuICAgICAgICB9XG4gICAgICAgIGlmIChiID4gazIgJiYgYSA8PSBrMikge1xuICAgICAgICAgICAgLy8gfCAgLS0tfC0tPiBvciAtLS18LS0tLS18LS0+IChsaW5lIGV4aXRzIHRoZSBjbGlwIHJlZ2lvbiBvbiB0aGUgcmlnaHQpXG4gICAgICAgICAgICB0ID0gaW50ZXJzZWN0KHNsaWNlLCBheCwgYXksIGJ4LCBieSwgazIpO1xuICAgICAgICAgICAgZXhpdGVkID0gdHJ1ZTtcbiAgICAgICAgfVxuXG4gICAgICAgIGlmICghaXNQb2x5Z29uICYmIGV4aXRlZCkge1xuICAgICAgICAgICAgaWYgKHRyYWNrTWV0cmljcykgc2xpY2UuZW5kID0gbGVuICsgc2VnTGVuICogdDtcbiAgICAgICAgICAgIG5ld0dlb20ucHVzaChzbGljZSk7XG4gICAgICAgICAgICBzbGljZSA9IG5ld1NsaWNlKGdlb20pO1xuICAgICAgICB9XG5cbiAgICAgICAgaWYgKHRyYWNrTWV0cmljcykgbGVuICs9IHNlZ0xlbjtcbiAgICB9XG5cbiAgICAvLyBhZGQgdGhlIGxhc3QgcG9pbnRcbiAgICB2YXIgbGFzdCA9IGdlb20ubGVuZ3RoIC0gMztcbiAgICBheCA9IGdlb21bbGFzdF07XG4gICAgYXkgPSBnZW9tW2xhc3QgKyAxXTtcbiAgICBheiA9IGdlb21bbGFzdCArIDJdO1xuICAgIGEgPSBheGlzID09PSAwID8gYXggOiBheTtcbiAgICBpZiAoYSA+PSBrMSAmJiBhIDw9IGsyKSBhZGRQb2ludChzbGljZSwgYXgsIGF5LCBheik7XG5cbiAgICAvLyBjbG9zZSB0aGUgcG9seWdvbiBpZiBpdHMgZW5kcG9pbnRzIGFyZSBub3QgdGhlIHNhbWUgYWZ0ZXIgY2xpcHBpbmdcbiAgICBsYXN0ID0gc2xpY2UubGVuZ3RoIC0gMztcbiAgICBpZiAoaXNQb2x5Z29uICYmIGxhc3QgPj0gMyAmJiAoc2xpY2VbbGFzdF0gIT09IHNsaWNlWzBdIHx8IHNsaWNlW2xhc3QgKyAxXSAhPT0gc2xpY2VbMV0pKSB7XG4gICAgICAgIGFkZFBvaW50KHNsaWNlLCBzbGljZVswXSwgc2xpY2VbMV0sIHNsaWNlWzJdKTtcbiAgICB9XG5cbiAgICAvLyBhZGQgdGhlIGZpbmFsIHNsaWNlXG4gICAgaWYgKHNsaWNlLmxlbmd0aCkge1xuICAgICAgICBuZXdHZW9tLnB1c2goc2xpY2UpO1xuICAgIH1cbn1cblxuZnVuY3Rpb24gbmV3U2xpY2UobGluZSkge1xuICAgIHZhciBzbGljZSA9IFtdO1xuICAgIHNsaWNlLnNpemUgPSBsaW5lLnNpemU7XG4gICAgc2xpY2Uuc3RhcnQgPSBsaW5lLnN0YXJ0O1xuICAgIHNsaWNlLmVuZCA9IGxpbmUuZW5kO1xuICAgIHJldHVybiBzbGljZTtcbn1cblxuZnVuY3Rpb24gY2xpcExpbmVzKGdlb20sIG5ld0dlb20sIGsxLCBrMiwgYXhpcywgaXNQb2x5Z29uKSB7XG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBnZW9tLmxlbmd0aDsgaSsrKSB7XG4gICAgICAgIGNsaXBMaW5lKGdlb21baV0sIG5ld0dlb20sIGsxLCBrMiwgYXhpcywgaXNQb2x5Z29uLCBmYWxzZSk7XG4gICAgfVxufVxuXG5mdW5jdGlvbiBhZGRQb2ludChvdXQsIHgsIHksIHopIHtcbiAgICBvdXQucHVzaCh4KTtcbiAgICBvdXQucHVzaCh5KTtcbiAgICBvdXQucHVzaCh6KTtcbn1cblxuZnVuY3Rpb24gaW50ZXJzZWN0WChvdXQsIGF4LCBheSwgYngsIGJ5LCB4KSB7XG4gICAgdmFyIHQgPSAoeCAtIGF4KSAvIChieCAtIGF4KTtcbiAgICBvdXQucHVzaCh4KTtcbiAgICBvdXQucHVzaChheSArIChieSAtIGF5KSAqIHQpO1xuICAgIG91dC5wdXNoKDEpO1xuICAgIHJldHVybiB0O1xufVxuXG5mdW5jdGlvbiBpbnRlcnNlY3RZKG91dCwgYXgsIGF5LCBieCwgYnksIHkpIHtcbiAgICB2YXIgdCA9ICh5IC0gYXkpIC8gKGJ5IC0gYXkpO1xuICAgIG91dC5wdXNoKGF4ICsgKGJ4IC0gYXgpICogdCk7XG4gICAgb3V0LnB1c2goeSk7XG4gICAgb3V0LnB1c2goMSk7XG4gICAgcmV0dXJuIHQ7XG59XG4iLCJcbmltcG9ydCBjbGlwIGZyb20gJy4vY2xpcCc7XG5pbXBvcnQgY3JlYXRlRmVhdHVyZSBmcm9tICcuL2ZlYXR1cmUnO1xuXG5leHBvcnQgZGVmYXVsdCBmdW5jdGlvbiB3cmFwKGZlYXR1cmVzLCBvcHRpb25zKSB7XG4gICAgdmFyIGJ1ZmZlciA9IG9wdGlvbnMuYnVmZmVyIC8gb3B0aW9ucy5leHRlbnQ7XG4gICAgdmFyIG1lcmdlZCA9IGZlYXR1cmVzO1xuICAgIHZhciBsZWZ0ICA9IGNsaXAoZmVhdHVyZXMsIDEsIC0xIC0gYnVmZmVyLCBidWZmZXIsICAgICAwLCAtMSwgMiwgb3B0aW9ucyk7IC8vIGxlZnQgd29ybGQgY29weVxuICAgIHZhciByaWdodCA9IGNsaXAoZmVhdHVyZXMsIDEsICAxIC0gYnVmZmVyLCAyICsgYnVmZmVyLCAwLCAtMSwgMiwgb3B0aW9ucyk7IC8vIHJpZ2h0IHdvcmxkIGNvcHlcblxuICAgIGlmIChsZWZ0IHx8IHJpZ2h0KSB7XG4gICAgICAgIG1lcmdlZCA9IGNsaXAoZmVhdHVyZXMsIDEsIC1idWZmZXIsIDEgKyBidWZmZXIsIDAsIC0xLCAyLCBvcHRpb25zKSB8fCBbXTsgLy8gY2VudGVyIHdvcmxkIGNvcHlcblxuICAgICAgICBpZiAobGVmdCkgbWVyZ2VkID0gc2hpZnRGZWF0dXJlQ29vcmRzKGxlZnQsIDEpLmNvbmNhdChtZXJnZWQpOyAvLyBtZXJnZSBsZWZ0IGludG8gY2VudGVyXG4gICAgICAgIGlmIChyaWdodCkgbWVyZ2VkID0gbWVyZ2VkLmNvbmNhdChzaGlmdEZlYXR1cmVDb29yZHMocmlnaHQsIC0xKSk7IC8vIG1lcmdlIHJpZ2h0IGludG8gY2VudGVyXG4gICAgfVxuXG4gICAgcmV0dXJuIG1lcmdlZDtcbn1cblxuZnVuY3Rpb24gc2hpZnRGZWF0dXJlQ29vcmRzKGZlYXR1cmVzLCBvZmZzZXQpIHtcbiAgICB2YXIgbmV3RmVhdHVyZXMgPSBbXTtcblxuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgZmVhdHVyZXMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgdmFyIGZlYXR1cmUgPSBmZWF0dXJlc1tpXSxcbiAgICAgICAgICAgIHR5cGUgPSBmZWF0dXJlLnR5cGU7XG5cbiAgICAgICAgdmFyIG5ld0dlb21ldHJ5O1xuXG4gICAgICAgIGlmICh0eXBlID09PSAnUG9pbnQnIHx8IHR5cGUgPT09ICdNdWx0aVBvaW50JyB8fCB0eXBlID09PSAnTGluZVN0cmluZycpIHtcbiAgICAgICAgICAgIG5ld0dlb21ldHJ5ID0gc2hpZnRDb29yZHMoZmVhdHVyZS5nZW9tZXRyeSwgb2Zmc2V0KTtcblxuICAgICAgICB9IGVsc2UgaWYgKHR5cGUgPT09ICdNdWx0aUxpbmVTdHJpbmcnIHx8IHR5cGUgPT09ICdQb2x5Z29uJykge1xuICAgICAgICAgICAgbmV3R2VvbWV0cnkgPSBbXTtcbiAgICAgICAgICAgIGZvciAodmFyIGogPSAwOyBqIDwgZmVhdHVyZS5nZW9tZXRyeS5sZW5ndGg7IGorKykge1xuICAgICAgICAgICAgICAgIG5ld0dlb21ldHJ5LnB1c2goc2hpZnRDb29yZHMoZmVhdHVyZS5nZW9tZXRyeVtqXSwgb2Zmc2V0KSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH0gZWxzZSBpZiAodHlwZSA9PT0gJ011bHRpUG9seWdvbicpIHtcbiAgICAgICAgICAgIG5ld0dlb21ldHJ5ID0gW107XG4gICAgICAgICAgICBmb3IgKGogPSAwOyBqIDwgZmVhdHVyZS5nZW9tZXRyeS5sZW5ndGg7IGorKykge1xuICAgICAgICAgICAgICAgIHZhciBuZXdQb2x5Z29uID0gW107XG4gICAgICAgICAgICAgICAgZm9yICh2YXIgayA9IDA7IGsgPCBmZWF0dXJlLmdlb21ldHJ5W2pdLmxlbmd0aDsgaysrKSB7XG4gICAgICAgICAgICAgICAgICAgIG5ld1BvbHlnb24ucHVzaChzaGlmdENvb3JkcyhmZWF0dXJlLmdlb21ldHJ5W2pdW2tdLCBvZmZzZXQpKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgbmV3R2VvbWV0cnkucHVzaChuZXdQb2x5Z29uKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuXG4gICAgICAgIG5ld0ZlYXR1cmVzLnB1c2goY3JlYXRlRmVhdHVyZShmZWF0dXJlLmlkLCB0eXBlLCBuZXdHZW9tZXRyeSwgZmVhdHVyZS50YWdzKSk7XG4gICAgfVxuXG4gICAgcmV0dXJuIG5ld0ZlYXR1cmVzO1xufVxuXG5mdW5jdGlvbiBzaGlmdENvb3Jkcyhwb2ludHMsIG9mZnNldCkge1xuICAgIHZhciBuZXdQb2ludHMgPSBbXTtcbiAgICBuZXdQb2ludHMuc2l6ZSA9IHBvaW50cy5zaXplO1xuXG4gICAgaWYgKHBvaW50cy5zdGFydCAhPT0gdW5kZWZpbmVkKSB7XG4gICAgICAgIG5ld1BvaW50cy5zdGFydCA9IHBvaW50cy5zdGFydDtcbiAgICAgICAgbmV3UG9pbnRzLmVuZCA9IHBvaW50cy5lbmQ7XG4gICAgfVxuXG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBwb2ludHMubGVuZ3RoOyBpICs9IDMpIHtcbiAgICAgICAgbmV3UG9pbnRzLnB1c2gocG9pbnRzW2ldICsgb2Zmc2V0LCBwb2ludHNbaSArIDFdLCBwb2ludHNbaSArIDJdKTtcbiAgICB9XG4gICAgcmV0dXJuIG5ld1BvaW50cztcbn1cbiIsIlxuLy8gVHJhbnNmb3JtcyB0aGUgY29vcmRpbmF0ZXMgb2YgZWFjaCBmZWF0dXJlIGluIHRoZSBnaXZlbiB0aWxlIGZyb21cbi8vIG1lcmNhdG9yLXByb2plY3RlZCBzcGFjZSBpbnRvIChleHRlbnQgeCBleHRlbnQpIHRpbGUgc3BhY2UuXG5leHBvcnQgZGVmYXVsdCBmdW5jdGlvbiB0cmFuc2Zvcm1UaWxlKHRpbGUsIGV4dGVudCkge1xuICAgIGlmICh0aWxlLnRyYW5zZm9ybWVkKSByZXR1cm4gdGlsZTtcblxuICAgIHZhciB6MiA9IDEgPDwgdGlsZS56LFxuICAgICAgICB0eCA9IHRpbGUueCxcbiAgICAgICAgdHkgPSB0aWxlLnksXG4gICAgICAgIGksIGosIGs7XG5cbiAgICBmb3IgKGkgPSAwOyBpIDwgdGlsZS5mZWF0dXJlcy5sZW5ndGg7IGkrKykge1xuICAgICAgICB2YXIgZmVhdHVyZSA9IHRpbGUuZmVhdHVyZXNbaV0sXG4gICAgICAgICAgICBnZW9tID0gZmVhdHVyZS5nZW9tZXRyeSxcbiAgICAgICAgICAgIHR5cGUgPSBmZWF0dXJlLnR5cGU7XG5cbiAgICAgICAgZmVhdHVyZS5nZW9tZXRyeSA9IFtdO1xuXG4gICAgICAgIGlmICh0eXBlID09PSAxKSB7XG4gICAgICAgICAgICBmb3IgKGogPSAwOyBqIDwgZ2VvbS5sZW5ndGg7IGogKz0gMikge1xuICAgICAgICAgICAgICAgIGZlYXR1cmUuZ2VvbWV0cnkucHVzaCh0cmFuc2Zvcm1Qb2ludChnZW9tW2pdLCBnZW9tW2ogKyAxXSwgZXh0ZW50LCB6MiwgdHgsIHR5KSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBmb3IgKGogPSAwOyBqIDwgZ2VvbS5sZW5ndGg7IGorKykge1xuICAgICAgICAgICAgICAgIHZhciByaW5nID0gW107XG4gICAgICAgICAgICAgICAgZm9yIChrID0gMDsgayA8IGdlb21bal0ubGVuZ3RoOyBrICs9IDIpIHtcbiAgICAgICAgICAgICAgICAgICAgcmluZy5wdXNoKHRyYW5zZm9ybVBvaW50KGdlb21bal1ba10sIGdlb21bal1bayArIDFdLCBleHRlbnQsIHoyLCB0eCwgdHkpKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgZmVhdHVyZS5nZW9tZXRyeS5wdXNoKHJpbmcpO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgfVxuXG4gICAgdGlsZS50cmFuc2Zvcm1lZCA9IHRydWU7XG5cbiAgICByZXR1cm4gdGlsZTtcbn1cblxuZnVuY3Rpb24gdHJhbnNmb3JtUG9pbnQoeCwgeSwgZXh0ZW50LCB6MiwgdHgsIHR5KSB7XG4gICAgcmV0dXJuIFtcbiAgICAgICAgTWF0aC5yb3VuZChleHRlbnQgKiAoeCAqIHoyIC0gdHgpKSxcbiAgICAgICAgTWF0aC5yb3VuZChleHRlbnQgKiAoeSAqIHoyIC0gdHkpKV07XG59XG4iLCJcbmV4cG9ydCBkZWZhdWx0IGZ1bmN0aW9uIGNyZWF0ZVRpbGUoZmVhdHVyZXMsIHosIHR4LCB0eSwgb3B0aW9ucykge1xuICAgIHZhciB0b2xlcmFuY2UgPSB6ID09PSBvcHRpb25zLm1heFpvb20gPyAwIDogb3B0aW9ucy50b2xlcmFuY2UgLyAoKDEgPDwgeikgKiBvcHRpb25zLmV4dGVudCk7XG4gICAgdmFyIHRpbGUgPSB7XG4gICAgICAgIGZlYXR1cmVzOiBbXSxcbiAgICAgICAgbnVtUG9pbnRzOiAwLFxuICAgICAgICBudW1TaW1wbGlmaWVkOiAwLFxuICAgICAgICBudW1GZWF0dXJlczogMCxcbiAgICAgICAgc291cmNlOiBudWxsLFxuICAgICAgICB4OiB0eCxcbiAgICAgICAgeTogdHksXG4gICAgICAgIHo6IHosXG4gICAgICAgIHRyYW5zZm9ybWVkOiBmYWxzZSxcbiAgICAgICAgbWluWDogMixcbiAgICAgICAgbWluWTogMSxcbiAgICAgICAgbWF4WDogLTEsXG4gICAgICAgIG1heFk6IDBcbiAgICB9O1xuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgZmVhdHVyZXMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgdGlsZS5udW1GZWF0dXJlcysrO1xuICAgICAgICBhZGRGZWF0dXJlKHRpbGUsIGZlYXR1cmVzW2ldLCB0b2xlcmFuY2UsIG9wdGlvbnMpO1xuXG4gICAgICAgIHZhciBtaW5YID0gZmVhdHVyZXNbaV0ubWluWDtcbiAgICAgICAgdmFyIG1pblkgPSBmZWF0dXJlc1tpXS5taW5ZO1xuICAgICAgICB2YXIgbWF4WCA9IGZlYXR1cmVzW2ldLm1heFg7XG4gICAgICAgIHZhciBtYXhZID0gZmVhdHVyZXNbaV0ubWF4WTtcblxuICAgICAgICBpZiAobWluWCA8IHRpbGUubWluWCkgdGlsZS5taW5YID0gbWluWDtcbiAgICAgICAgaWYgKG1pblkgPCB0aWxlLm1pblkpIHRpbGUubWluWSA9IG1pblk7XG4gICAgICAgIGlmIChtYXhYID4gdGlsZS5tYXhYKSB0aWxlLm1heFggPSBtYXhYO1xuICAgICAgICBpZiAobWF4WSA+IHRpbGUubWF4WSkgdGlsZS5tYXhZID0gbWF4WTtcbiAgICB9XG4gICAgcmV0dXJuIHRpbGU7XG59XG5cbmZ1bmN0aW9uIGFkZEZlYXR1cmUodGlsZSwgZmVhdHVyZSwgdG9sZXJhbmNlLCBvcHRpb25zKSB7XG5cbiAgICB2YXIgZ2VvbSA9IGZlYXR1cmUuZ2VvbWV0cnksXG4gICAgICAgIHR5cGUgPSBmZWF0dXJlLnR5cGUsXG4gICAgICAgIHNpbXBsaWZpZWQgPSBbXTtcblxuICAgIGlmICh0eXBlID09PSAnUG9pbnQnIHx8IHR5cGUgPT09ICdNdWx0aVBvaW50Jykge1xuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IGdlb20ubGVuZ3RoOyBpICs9IDMpIHtcbiAgICAgICAgICAgIHNpbXBsaWZpZWQucHVzaChnZW9tW2ldKTtcbiAgICAgICAgICAgIHNpbXBsaWZpZWQucHVzaChnZW9tW2kgKyAxXSk7XG4gICAgICAgICAgICB0aWxlLm51bVBvaW50cysrO1xuICAgICAgICAgICAgdGlsZS5udW1TaW1wbGlmaWVkKys7XG4gICAgICAgIH1cblxuICAgIH0gZWxzZSBpZiAodHlwZSA9PT0gJ0xpbmVTdHJpbmcnKSB7XG4gICAgICAgIGFkZExpbmUoc2ltcGxpZmllZCwgZ2VvbSwgdGlsZSwgdG9sZXJhbmNlLCBmYWxzZSwgZmFsc2UpO1xuXG4gICAgfSBlbHNlIGlmICh0eXBlID09PSAnTXVsdGlMaW5lU3RyaW5nJyB8fCB0eXBlID09PSAnUG9seWdvbicpIHtcbiAgICAgICAgZm9yIChpID0gMDsgaSA8IGdlb20ubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgICAgIGFkZExpbmUoc2ltcGxpZmllZCwgZ2VvbVtpXSwgdGlsZSwgdG9sZXJhbmNlLCB0eXBlID09PSAnUG9seWdvbicsIGkgPT09IDApO1xuICAgICAgICB9XG5cbiAgICB9IGVsc2UgaWYgKHR5cGUgPT09ICdNdWx0aVBvbHlnb24nKSB7XG5cbiAgICAgICAgZm9yICh2YXIgayA9IDA7IGsgPCBnZW9tLmxlbmd0aDsgaysrKSB7XG4gICAgICAgICAgICB2YXIgcG9seWdvbiA9IGdlb21ba107XG4gICAgICAgICAgICBmb3IgKGkgPSAwOyBpIDwgcG9seWdvbi5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICAgICAgIGFkZExpbmUoc2ltcGxpZmllZCwgcG9seWdvbltpXSwgdGlsZSwgdG9sZXJhbmNlLCB0cnVlLCBpID09PSAwKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH1cblxuICAgIGlmIChzaW1wbGlmaWVkLmxlbmd0aCkge1xuICAgICAgICB2YXIgdGFncyA9IGZlYXR1cmUudGFncyB8fCBudWxsO1xuICAgICAgICBpZiAodHlwZSA9PT0gJ0xpbmVTdHJpbmcnICYmIG9wdGlvbnMubGluZU1ldHJpY3MpIHtcbiAgICAgICAgICAgIHRhZ3MgPSB7fTtcbiAgICAgICAgICAgIGZvciAodmFyIGtleSBpbiBmZWF0dXJlLnRhZ3MpIHRhZ3Nba2V5XSA9IGZlYXR1cmUudGFnc1trZXldO1xuICAgICAgICAgICAgdGFnc1snbWFwYm94X2NsaXBfc3RhcnQnXSA9IGdlb20uc3RhcnQgLyBnZW9tLnNpemU7XG4gICAgICAgICAgICB0YWdzWydtYXBib3hfY2xpcF9lbmQnXSA9IGdlb20uZW5kIC8gZ2VvbS5zaXplO1xuICAgICAgICB9XG4gICAgICAgIHZhciB0aWxlRmVhdHVyZSA9IHtcbiAgICAgICAgICAgIGdlb21ldHJ5OiBzaW1wbGlmaWVkLFxuICAgICAgICAgICAgdHlwZTogdHlwZSA9PT0gJ1BvbHlnb24nIHx8IHR5cGUgPT09ICdNdWx0aVBvbHlnb24nID8gMyA6XG4gICAgICAgICAgICAgICAgdHlwZSA9PT0gJ0xpbmVTdHJpbmcnIHx8IHR5cGUgPT09ICdNdWx0aUxpbmVTdHJpbmcnID8gMiA6IDEsXG4gICAgICAgICAgICB0YWdzOiB0YWdzXG4gICAgICAgIH07XG4gICAgICAgIGlmIChmZWF0dXJlLmlkICE9PSBudWxsKSB7XG4gICAgICAgICAgICB0aWxlRmVhdHVyZS5pZCA9IGZlYXR1cmUuaWQ7XG4gICAgICAgIH1cbiAgICAgICAgdGlsZS5mZWF0dXJlcy5wdXNoKHRpbGVGZWF0dXJlKTtcbiAgICB9XG59XG5cbmZ1bmN0aW9uIGFkZExpbmUocmVzdWx0LCBnZW9tLCB0aWxlLCB0b2xlcmFuY2UsIGlzUG9seWdvbiwgaXNPdXRlcikge1xuICAgIHZhciBzcVRvbGVyYW5jZSA9IHRvbGVyYW5jZSAqIHRvbGVyYW5jZTtcblxuICAgIGlmICh0b2xlcmFuY2UgPiAwICYmIChnZW9tLnNpemUgPCAoaXNQb2x5Z29uID8gc3FUb2xlcmFuY2UgOiB0b2xlcmFuY2UpKSkge1xuICAgICAgICB0aWxlLm51bVBvaW50cyArPSBnZW9tLmxlbmd0aCAvIDM7XG4gICAgICAgIHJldHVybjtcbiAgICB9XG5cbiAgICB2YXIgcmluZyA9IFtdO1xuXG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBnZW9tLmxlbmd0aDsgaSArPSAzKSB7XG4gICAgICAgIGlmICh0b2xlcmFuY2UgPT09IDAgfHwgZ2VvbVtpICsgMl0gPiBzcVRvbGVyYW5jZSkge1xuICAgICAgICAgICAgdGlsZS5udW1TaW1wbGlmaWVkKys7XG4gICAgICAgICAgICByaW5nLnB1c2goZ2VvbVtpXSk7XG4gICAgICAgICAgICByaW5nLnB1c2goZ2VvbVtpICsgMV0pO1xuICAgICAgICB9XG4gICAgICAgIHRpbGUubnVtUG9pbnRzKys7XG4gICAgfVxuXG4gICAgaWYgKGlzUG9seWdvbikgcmV3aW5kKHJpbmcsIGlzT3V0ZXIpO1xuXG4gICAgcmVzdWx0LnB1c2gocmluZyk7XG59XG5cbmZ1bmN0aW9uIHJld2luZChyaW5nLCBjbG9ja3dpc2UpIHtcbiAgICB2YXIgYXJlYSA9IDA7XG4gICAgZm9yICh2YXIgaSA9IDAsIGxlbiA9IHJpbmcubGVuZ3RoLCBqID0gbGVuIC0gMjsgaSA8IGxlbjsgaiA9IGksIGkgKz0gMikge1xuICAgICAgICBhcmVhICs9IChyaW5nW2ldIC0gcmluZ1tqXSkgKiAocmluZ1tpICsgMV0gKyByaW5nW2ogKyAxXSk7XG4gICAgfVxuICAgIGlmIChhcmVhID4gMCA9PT0gY2xvY2t3aXNlKSB7XG4gICAgICAgIGZvciAoaSA9IDAsIGxlbiA9IHJpbmcubGVuZ3RoOyBpIDwgbGVuIC8gMjsgaSArPSAyKSB7XG4gICAgICAgICAgICB2YXIgeCA9IHJpbmdbaV07XG4gICAgICAgICAgICB2YXIgeSA9IHJpbmdbaSArIDFdO1xuICAgICAgICAgICAgcmluZ1tpXSA9IHJpbmdbbGVuIC0gMiAtIGldO1xuICAgICAgICAgICAgcmluZ1tpICsgMV0gPSByaW5nW2xlbiAtIDEgLSBpXTtcbiAgICAgICAgICAgIHJpbmdbbGVuIC0gMiAtIGldID0geDtcbiAgICAgICAgICAgIHJpbmdbbGVuIC0gMSAtIGldID0geTtcbiAgICAgICAgfVxuICAgIH1cbn1cbiIsIlxuaW1wb3J0IGNvbnZlcnQgZnJvbSAnLi9jb252ZXJ0JzsgICAgIC8vIEdlb0pTT04gY29udmVyc2lvbiBhbmQgcHJlcHJvY2Vzc2luZ1xuaW1wb3J0IGNsaXAgZnJvbSAnLi9jbGlwJzsgICAgICAgICAgIC8vIHN0cmlwZSBjbGlwcGluZyBhbGdvcml0aG1cbmltcG9ydCB3cmFwIGZyb20gJy4vd3JhcCc7ICAgICAgICAgICAvLyBkYXRlIGxpbmUgcHJvY2Vzc2luZ1xuaW1wb3J0IHRyYW5zZm9ybSBmcm9tICcuL3RyYW5zZm9ybSc7IC8vIGNvb3JkaW5hdGUgdHJhbnNmb3JtYXRpb25cbmltcG9ydCBjcmVhdGVUaWxlIGZyb20gJy4vdGlsZSc7ICAgICAvLyBmaW5hbCBzaW1wbGlmaWVkIHRpbGUgZ2VuZXJhdGlvblxuXG5leHBvcnQgZGVmYXVsdCBmdW5jdGlvbiBnZW9qc29udnQoZGF0YSwgb3B0aW9ucykge1xuICAgIHJldHVybiBuZXcgR2VvSlNPTlZUKGRhdGEsIG9wdGlvbnMpO1xufVxuXG5mdW5jdGlvbiBHZW9KU09OVlQoZGF0YSwgb3B0aW9ucykge1xuICAgIG9wdGlvbnMgPSB0aGlzLm9wdGlvbnMgPSBleHRlbmQoT2JqZWN0LmNyZWF0ZSh0aGlzLm9wdGlvbnMpLCBvcHRpb25zKTtcblxuICAgIHZhciBkZWJ1ZyA9IG9wdGlvbnMuZGVidWc7XG5cbiAgICBpZiAoZGVidWcpIGNvbnNvbGUudGltZSgncHJlcHJvY2VzcyBkYXRhJyk7XG5cbiAgICBpZiAob3B0aW9ucy5tYXhab29tIDwgMCB8fCBvcHRpb25zLm1heFpvb20gPiAyNCkgdGhyb3cgbmV3IEVycm9yKCdtYXhab29tIHNob3VsZCBiZSBpbiB0aGUgMC0yNCByYW5nZScpO1xuICAgIGlmIChvcHRpb25zLnByb21vdGVJZCAmJiBvcHRpb25zLmdlbmVyYXRlSWQpIHRocm93IG5ldyBFcnJvcigncHJvbW90ZUlkIGFuZCBnZW5lcmF0ZUlkIGNhbm5vdCBiZSB1c2VkIHRvZ2V0aGVyLicpO1xuXG4gICAgdmFyIGZlYXR1cmVzID0gY29udmVydChkYXRhLCBvcHRpb25zKTtcblxuICAgIHRoaXMudGlsZXMgPSB7fTtcbiAgICB0aGlzLnRpbGVDb29yZHMgPSBbXTtcblxuICAgIGlmIChkZWJ1Zykge1xuICAgICAgICBjb25zb2xlLnRpbWVFbmQoJ3ByZXByb2Nlc3MgZGF0YScpO1xuICAgICAgICBjb25zb2xlLmxvZygnaW5kZXg6IG1heFpvb206ICVkLCBtYXhQb2ludHM6ICVkJywgb3B0aW9ucy5pbmRleE1heFpvb20sIG9wdGlvbnMuaW5kZXhNYXhQb2ludHMpO1xuICAgICAgICBjb25zb2xlLnRpbWUoJ2dlbmVyYXRlIHRpbGVzJyk7XG4gICAgICAgIHRoaXMuc3RhdHMgPSB7fTtcbiAgICAgICAgdGhpcy50b3RhbCA9IDA7XG4gICAgfVxuXG4gICAgZmVhdHVyZXMgPSB3cmFwKGZlYXR1cmVzLCBvcHRpb25zKTtcblxuICAgIC8vIHN0YXJ0IHNsaWNpbmcgZnJvbSB0aGUgdG9wIHRpbGUgZG93blxuICAgIGlmIChmZWF0dXJlcy5sZW5ndGgpIHRoaXMuc3BsaXRUaWxlKGZlYXR1cmVzLCAwLCAwLCAwKTtcblxuICAgIGlmIChkZWJ1Zykge1xuICAgICAgICBpZiAoZmVhdHVyZXMubGVuZ3RoKSBjb25zb2xlLmxvZygnZmVhdHVyZXM6ICVkLCBwb2ludHM6ICVkJywgdGhpcy50aWxlc1swXS5udW1GZWF0dXJlcywgdGhpcy50aWxlc1swXS5udW1Qb2ludHMpO1xuICAgICAgICBjb25zb2xlLnRpbWVFbmQoJ2dlbmVyYXRlIHRpbGVzJyk7XG4gICAgICAgIGNvbnNvbGUubG9nKCd0aWxlcyBnZW5lcmF0ZWQ6JywgdGhpcy50b3RhbCwgSlNPTi5zdHJpbmdpZnkodGhpcy5zdGF0cykpO1xuICAgIH1cbn1cblxuR2VvSlNPTlZULnByb3RvdHlwZS5vcHRpb25zID0ge1xuICAgIG1heFpvb206IDE0LCAgICAgICAgICAgIC8vIG1heCB6b29tIHRvIHByZXNlcnZlIGRldGFpbCBvblxuICAgIGluZGV4TWF4Wm9vbTogNSwgICAgICAgIC8vIG1heCB6b29tIGluIHRoZSB0aWxlIGluZGV4XG4gICAgaW5kZXhNYXhQb2ludHM6IDEwMDAwMCwgLy8gbWF4IG51bWJlciBvZiBwb2ludHMgcGVyIHRpbGUgaW4gdGhlIHRpbGUgaW5kZXhcbiAgICB0b2xlcmFuY2U6IDMsICAgICAgICAgICAvLyBzaW1wbGlmaWNhdGlvbiB0b2xlcmFuY2UgKGhpZ2hlciBtZWFucyBzaW1wbGVyKVxuICAgIGV4dGVudDogNDA5NiwgICAgICAgICAgIC8vIHRpbGUgZXh0ZW50XG4gICAgYnVmZmVyOiA2NCwgICAgICAgICAgICAgLy8gdGlsZSBidWZmZXIgb24gZWFjaCBzaWRlXG4gICAgbGluZU1ldHJpY3M6IGZhbHNlLCAgICAgLy8gd2hldGhlciB0byBjYWxjdWxhdGUgbGluZSBtZXRyaWNzXG4gICAgcHJvbW90ZUlkOiBudWxsLCAgICAgICAgLy8gbmFtZSBvZiBhIGZlYXR1cmUgcHJvcGVydHkgdG8gYmUgcHJvbW90ZWQgdG8gZmVhdHVyZS5pZFxuICAgIGdlbmVyYXRlSWQ6IGZhbHNlLCAgICAgIC8vIHdoZXRoZXIgdG8gZ2VuZXJhdGUgZmVhdHVyZSBpZHMuIENhbm5vdCBiZSB1c2VkIHdpdGggcHJvbW90ZUlkXG4gICAgZGVidWc6IDAgICAgICAgICAgICAgICAgLy8gbG9nZ2luZyBsZXZlbCAoMCwgMSBvciAyKVxufTtcblxuR2VvSlNPTlZULnByb3RvdHlwZS5zcGxpdFRpbGUgPSBmdW5jdGlvbiAoZmVhdHVyZXMsIHosIHgsIHksIGN6LCBjeCwgY3kpIHtcblxuICAgIHZhciBzdGFjayA9IFtmZWF0dXJlcywgeiwgeCwgeV0sXG4gICAgICAgIG9wdGlvbnMgPSB0aGlzLm9wdGlvbnMsXG4gICAgICAgIGRlYnVnID0gb3B0aW9ucy5kZWJ1ZztcblxuICAgIC8vIGF2b2lkIHJlY3Vyc2lvbiBieSB1c2luZyBhIHByb2Nlc3NpbmcgcXVldWVcbiAgICB3aGlsZSAoc3RhY2subGVuZ3RoKSB7XG4gICAgICAgIHkgPSBzdGFjay5wb3AoKTtcbiAgICAgICAgeCA9IHN0YWNrLnBvcCgpO1xuICAgICAgICB6ID0gc3RhY2sucG9wKCk7XG4gICAgICAgIGZlYXR1cmVzID0gc3RhY2sucG9wKCk7XG5cbiAgICAgICAgdmFyIHoyID0gMSA8PCB6LFxuICAgICAgICAgICAgaWQgPSB0b0lEKHosIHgsIHkpLFxuICAgICAgICAgICAgdGlsZSA9IHRoaXMudGlsZXNbaWRdO1xuXG4gICAgICAgIGlmICghdGlsZSkge1xuICAgICAgICAgICAgaWYgKGRlYnVnID4gMSkgY29uc29sZS50aW1lKCdjcmVhdGlvbicpO1xuXG4gICAgICAgICAgICB0aWxlID0gdGhpcy50aWxlc1tpZF0gPSBjcmVhdGVUaWxlKGZlYXR1cmVzLCB6LCB4LCB5LCBvcHRpb25zKTtcbiAgICAgICAgICAgIHRoaXMudGlsZUNvb3Jkcy5wdXNoKHt6OiB6LCB4OiB4LCB5OiB5fSk7XG5cbiAgICAgICAgICAgIGlmIChkZWJ1Zykge1xuICAgICAgICAgICAgICAgIGlmIChkZWJ1ZyA+IDEpIHtcbiAgICAgICAgICAgICAgICAgICAgY29uc29sZS5sb2coJ3RpbGUgeiVkLSVkLSVkIChmZWF0dXJlczogJWQsIHBvaW50czogJWQsIHNpbXBsaWZpZWQ6ICVkKScsXG4gICAgICAgICAgICAgICAgICAgICAgICB6LCB4LCB5LCB0aWxlLm51bUZlYXR1cmVzLCB0aWxlLm51bVBvaW50cywgdGlsZS5udW1TaW1wbGlmaWVkKTtcbiAgICAgICAgICAgICAgICAgICAgY29uc29sZS50aW1lRW5kKCdjcmVhdGlvbicpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB2YXIga2V5ID0gJ3onICsgejtcbiAgICAgICAgICAgICAgICB0aGlzLnN0YXRzW2tleV0gPSAodGhpcy5zdGF0c1trZXldIHx8IDApICsgMTtcbiAgICAgICAgICAgICAgICB0aGlzLnRvdGFsKys7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICAvLyBzYXZlIHJlZmVyZW5jZSB0byBvcmlnaW5hbCBnZW9tZXRyeSBpbiB0aWxlIHNvIHRoYXQgd2UgY2FuIGRyaWxsIGRvd24gbGF0ZXIgaWYgd2Ugc3RvcCBub3dcbiAgICAgICAgdGlsZS5zb3VyY2UgPSBmZWF0dXJlcztcblxuICAgICAgICAvLyBpZiBpdCdzIHRoZSBmaXJzdC1wYXNzIHRpbGluZ1xuICAgICAgICBpZiAoIWN6KSB7XG4gICAgICAgICAgICAvLyBzdG9wIHRpbGluZyBpZiB3ZSByZWFjaGVkIG1heCB6b29tLCBvciBpZiB0aGUgdGlsZSBpcyB0b28gc2ltcGxlXG4gICAgICAgICAgICBpZiAoeiA9PT0gb3B0aW9ucy5pbmRleE1heFpvb20gfHwgdGlsZS5udW1Qb2ludHMgPD0gb3B0aW9ucy5pbmRleE1heFBvaW50cykgY29udGludWU7XG5cbiAgICAgICAgLy8gaWYgYSBkcmlsbGRvd24gdG8gYSBzcGVjaWZpYyB0aWxlXG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAvLyBzdG9wIHRpbGluZyBpZiB3ZSByZWFjaGVkIGJhc2Ugem9vbSBvciBvdXIgdGFyZ2V0IHRpbGUgem9vbVxuICAgICAgICAgICAgaWYgKHogPT09IG9wdGlvbnMubWF4Wm9vbSB8fCB6ID09PSBjeikgY29udGludWU7XG5cbiAgICAgICAgICAgIC8vIHN0b3AgdGlsaW5nIGlmIGl0J3Mgbm90IGFuIGFuY2VzdG9yIG9mIHRoZSB0YXJnZXQgdGlsZVxuICAgICAgICAgICAgdmFyIG0gPSAxIDw8IChjeiAtIHopO1xuICAgICAgICAgICAgaWYgKHggIT09IE1hdGguZmxvb3IoY3ggLyBtKSB8fCB5ICE9PSBNYXRoLmZsb29yKGN5IC8gbSkpIGNvbnRpbnVlO1xuICAgICAgICB9XG5cbiAgICAgICAgLy8gaWYgd2Ugc2xpY2UgZnVydGhlciBkb3duLCBubyBuZWVkIHRvIGtlZXAgc291cmNlIGdlb21ldHJ5XG4gICAgICAgIHRpbGUuc291cmNlID0gbnVsbDtcblxuICAgICAgICBpZiAoZmVhdHVyZXMubGVuZ3RoID09PSAwKSBjb250aW51ZTtcblxuICAgICAgICBpZiAoZGVidWcgPiAxKSBjb25zb2xlLnRpbWUoJ2NsaXBwaW5nJyk7XG5cbiAgICAgICAgLy8gdmFsdWVzIHdlJ2xsIHVzZSBmb3IgY2xpcHBpbmdcbiAgICAgICAgdmFyIGsxID0gMC41ICogb3B0aW9ucy5idWZmZXIgLyBvcHRpb25zLmV4dGVudCxcbiAgICAgICAgICAgIGsyID0gMC41IC0gazEsXG4gICAgICAgICAgICBrMyA9IDAuNSArIGsxLFxuICAgICAgICAgICAgazQgPSAxICsgazEsXG4gICAgICAgICAgICB0bCwgYmwsIHRyLCBiciwgbGVmdCwgcmlnaHQ7XG5cbiAgICAgICAgdGwgPSBibCA9IHRyID0gYnIgPSBudWxsO1xuXG4gICAgICAgIGxlZnQgID0gY2xpcChmZWF0dXJlcywgejIsIHggLSBrMSwgeCArIGszLCAwLCB0aWxlLm1pblgsIHRpbGUubWF4WCwgb3B0aW9ucyk7XG4gICAgICAgIHJpZ2h0ID0gY2xpcChmZWF0dXJlcywgejIsIHggKyBrMiwgeCArIGs0LCAwLCB0aWxlLm1pblgsIHRpbGUubWF4WCwgb3B0aW9ucyk7XG4gICAgICAgIGZlYXR1cmVzID0gbnVsbDtcblxuICAgICAgICBpZiAobGVmdCkge1xuICAgICAgICAgICAgdGwgPSBjbGlwKGxlZnQsIHoyLCB5IC0gazEsIHkgKyBrMywgMSwgdGlsZS5taW5ZLCB0aWxlLm1heFksIG9wdGlvbnMpO1xuICAgICAgICAgICAgYmwgPSBjbGlwKGxlZnQsIHoyLCB5ICsgazIsIHkgKyBrNCwgMSwgdGlsZS5taW5ZLCB0aWxlLm1heFksIG9wdGlvbnMpO1xuICAgICAgICAgICAgbGVmdCA9IG51bGw7XG4gICAgICAgIH1cblxuICAgICAgICBpZiAocmlnaHQpIHtcbiAgICAgICAgICAgIHRyID0gY2xpcChyaWdodCwgejIsIHkgLSBrMSwgeSArIGszLCAxLCB0aWxlLm1pblksIHRpbGUubWF4WSwgb3B0aW9ucyk7XG4gICAgICAgICAgICBiciA9IGNsaXAocmlnaHQsIHoyLCB5ICsgazIsIHkgKyBrNCwgMSwgdGlsZS5taW5ZLCB0aWxlLm1heFksIG9wdGlvbnMpO1xuICAgICAgICAgICAgcmlnaHQgPSBudWxsO1xuICAgICAgICB9XG5cbiAgICAgICAgaWYgKGRlYnVnID4gMSkgY29uc29sZS50aW1lRW5kKCdjbGlwcGluZycpO1xuXG4gICAgICAgIHN0YWNrLnB1c2godGwgfHwgW10sIHogKyAxLCB4ICogMiwgICAgIHkgKiAyKTtcbiAgICAgICAgc3RhY2sucHVzaChibCB8fCBbXSwgeiArIDEsIHggKiAyLCAgICAgeSAqIDIgKyAxKTtcbiAgICAgICAgc3RhY2sucHVzaCh0ciB8fCBbXSwgeiArIDEsIHggKiAyICsgMSwgeSAqIDIpO1xuICAgICAgICBzdGFjay5wdXNoKGJyIHx8IFtdLCB6ICsgMSwgeCAqIDIgKyAxLCB5ICogMiArIDEpO1xuICAgIH1cbn07XG5cbkdlb0pTT05WVC5wcm90b3R5cGUuZ2V0VGlsZSA9IGZ1bmN0aW9uICh6LCB4LCB5KSB7XG4gICAgdmFyIG9wdGlvbnMgPSB0aGlzLm9wdGlvbnMsXG4gICAgICAgIGV4dGVudCA9IG9wdGlvbnMuZXh0ZW50LFxuICAgICAgICBkZWJ1ZyA9IG9wdGlvbnMuZGVidWc7XG5cbiAgICBpZiAoeiA8IDAgfHwgeiA+IDI0KSByZXR1cm4gbnVsbDtcblxuICAgIHZhciB6MiA9IDEgPDwgejtcbiAgICB4ID0gKCh4ICUgejIpICsgejIpICUgejI7IC8vIHdyYXAgdGlsZSB4IGNvb3JkaW5hdGVcblxuICAgIHZhciBpZCA9IHRvSUQoeiwgeCwgeSk7XG4gICAgaWYgKHRoaXMudGlsZXNbaWRdKSByZXR1cm4gdHJhbnNmb3JtKHRoaXMudGlsZXNbaWRdLCBleHRlbnQpO1xuXG4gICAgaWYgKGRlYnVnID4gMSkgY29uc29sZS5sb2coJ2RyaWxsaW5nIGRvd24gdG8geiVkLSVkLSVkJywgeiwgeCwgeSk7XG5cbiAgICB2YXIgejAgPSB6LFxuICAgICAgICB4MCA9IHgsXG4gICAgICAgIHkwID0geSxcbiAgICAgICAgcGFyZW50O1xuXG4gICAgd2hpbGUgKCFwYXJlbnQgJiYgejAgPiAwKSB7XG4gICAgICAgIHowLS07XG4gICAgICAgIHgwID0gTWF0aC5mbG9vcih4MCAvIDIpO1xuICAgICAgICB5MCA9IE1hdGguZmxvb3IoeTAgLyAyKTtcbiAgICAgICAgcGFyZW50ID0gdGhpcy50aWxlc1t0b0lEKHowLCB4MCwgeTApXTtcbiAgICB9XG5cbiAgICBpZiAoIXBhcmVudCB8fCAhcGFyZW50LnNvdXJjZSkgcmV0dXJuIG51bGw7XG5cbiAgICAvLyBpZiB3ZSBmb3VuZCBhIHBhcmVudCB0aWxlIGNvbnRhaW5pbmcgdGhlIG9yaWdpbmFsIGdlb21ldHJ5LCB3ZSBjYW4gZHJpbGwgZG93biBmcm9tIGl0XG4gICAgaWYgKGRlYnVnID4gMSkgY29uc29sZS5sb2coJ2ZvdW5kIHBhcmVudCB0aWxlIHolZC0lZC0lZCcsIHowLCB4MCwgeTApO1xuXG4gICAgaWYgKGRlYnVnID4gMSkgY29uc29sZS50aW1lKCdkcmlsbGluZyBkb3duJyk7XG4gICAgdGhpcy5zcGxpdFRpbGUocGFyZW50LnNvdXJjZSwgejAsIHgwLCB5MCwgeiwgeCwgeSk7XG4gICAgaWYgKGRlYnVnID4gMSkgY29uc29sZS50aW1lRW5kKCdkcmlsbGluZyBkb3duJyk7XG5cbiAgICByZXR1cm4gdGhpcy50aWxlc1tpZF0gPyB0cmFuc2Zvcm0odGhpcy50aWxlc1tpZF0sIGV4dGVudCkgOiBudWxsO1xufTtcblxuZnVuY3Rpb24gdG9JRCh6LCB4LCB5KSB7XG4gICAgcmV0dXJuICgoKDEgPDwgeikgKiB5ICsgeCkgKiAzMikgKyB6O1xufVxuXG5mdW5jdGlvbiBleHRlbmQoZGVzdCwgc3JjKSB7XG4gICAgZm9yICh2YXIgaSBpbiBzcmMpIGRlc3RbaV0gPSBzcmNbaV07XG4gICAgcmV0dXJuIGRlc3Q7XG59XG4iLCIvLyBAZmxvd1xuXG5pbXBvcnQge2dldEpTT059IGZyb20gJy4uL3V0aWwvYWpheC5qcyc7XG5cbmltcG9ydCB7Z2V0UGVyZm9ybWFuY2VNZWFzdXJlbWVudH0gZnJvbSAnLi4vdXRpbC9wZXJmb3JtYW5jZS5qcyc7XG5pbXBvcnQgcmV3aW5kIGZyb20gJ0BtYXBib3gvZ2VvanNvbi1yZXdpbmQnO1xuaW1wb3J0IEdlb0pTT05XcmFwcGVyIGZyb20gJy4vZ2VvanNvbl93cmFwcGVyLmpzJztcbmltcG9ydCB2dHBiZiBmcm9tICd2dC1wYmYnO1xuaW1wb3J0IFN1cGVyY2x1c3RlciBmcm9tICdzdXBlcmNsdXN0ZXInO1xuaW1wb3J0IGdlb2pzb252dCBmcm9tICdnZW9qc29uLXZ0JztcbmltcG9ydCBhc3NlcnQgZnJvbSAnYXNzZXJ0JztcbmltcG9ydCBWZWN0b3JUaWxlV29ya2VyU291cmNlIGZyb20gJy4vdmVjdG9yX3RpbGVfd29ya2VyX3NvdXJjZS5qcyc7XG5pbXBvcnQge2NyZWF0ZUV4cHJlc3Npb259IGZyb20gJy4uL3N0eWxlLXNwZWMvZXhwcmVzc2lvbi9pbmRleC5qcyc7XG5cbmltcG9ydCB0eXBlIHtcbiAgICBSZXF1ZXN0ZWRUaWxlUGFyYW1ldGVycyxcbiAgICBXb3JrZXJUaWxlUGFyYW1ldGVycyxcbiAgICBXb3JrZXJUaWxlQ2FsbGJhY2ssXG59IGZyb20gJy4uL3NvdXJjZS93b3JrZXJfc291cmNlLmpzJztcblxuaW1wb3J0IHR5cGUgQWN0b3IgZnJvbSAnLi4vdXRpbC9hY3Rvci5qcyc7XG5pbXBvcnQgdHlwZSBTdHlsZUxheWVySW5kZXggZnJvbSAnLi4vc3R5bGUvc3R5bGVfbGF5ZXJfaW5kZXguanMnO1xuXG5pbXBvcnQgdHlwZSB7TG9hZFZlY3RvckRhdGFDYWxsYmFja30gZnJvbSAnLi92ZWN0b3JfdGlsZV93b3JrZXJfc291cmNlLmpzJztcbmltcG9ydCB0eXBlIHtSZXF1ZXN0UGFyYW1ldGVycywgUmVzcG9uc2VDYWxsYmFja30gZnJvbSAnLi4vdXRpbC9hamF4LmpzJztcbmltcG9ydCB0eXBlIHtDYWxsYmFja30gZnJvbSAnLi4vdHlwZXMvY2FsbGJhY2suanMnO1xuaW1wb3J0IHR5cGUge0dlb0pTT05GZWF0dXJlfSBmcm9tICdAbWFwYm94L2dlb2pzb24tdHlwZXMnO1xuXG5leHBvcnQgdHlwZSBHZW9KU09OV29ya2VyT3B0aW9ucyA9IHtcbiAgICBzb3VyY2U6IHN0cmluZyxcbiAgICBjbHVzdGVyOiBib29sZWFuLFxuICAgIHN1cGVyY2x1c3Rlck9wdGlvbnM/OiBPYmplY3QsXG4gICAgZ2VvanNvblZ0T3B0aW9ucz86IE9iamVjdCxcbiAgICBjbHVzdGVyUHJvcGVydGllcz86IE9iamVjdCxcbiAgICBmaWx0ZXI/OiBBcnJheTxtaXhlZD5cbn07XG5cbmV4cG9ydCB0eXBlIExvYWRHZW9KU09OUGFyYW1ldGVycyA9IEdlb0pTT05Xb3JrZXJPcHRpb25zICYge1xuICAgIHJlcXVlc3Q/OiBSZXF1ZXN0UGFyYW1ldGVycyxcbiAgICBkYXRhPzogc3RyaW5nXG59O1xuXG5leHBvcnQgdHlwZSBMb2FkR2VvSlNPTiA9IChwYXJhbXM6IExvYWRHZW9KU09OUGFyYW1ldGVycywgY2FsbGJhY2s6IFJlc3BvbnNlQ2FsbGJhY2s8T2JqZWN0PikgPT4gdm9pZDtcblxuZXhwb3J0IGludGVyZmFjZSBHZW9KU09OSW5kZXgge1xuICAgIGdldFRpbGUoejogbnVtYmVyLCB4OiBudW1iZXIsIHk6IG51bWJlcik6IE9iamVjdDtcblxuICAgIC8vIHN1cGVyY2x1c3RlciBtZXRob2RzXG4gICAgZ2V0Q2x1c3RlckV4cGFuc2lvblpvb20oY2x1c3RlcklkOiBudW1iZXIpOiBudW1iZXI7XG4gICAgZ2V0Q2hpbGRyZW4oY2x1c3RlcklkOiBudW1iZXIpOiBBcnJheTxHZW9KU09ORmVhdHVyZT47XG4gICAgZ2V0TGVhdmVzKGNsdXN0ZXJJZDogbnVtYmVyLCBsaW1pdDogbnVtYmVyLCBvZmZzZXQ6IG51bWJlcik6IEFycmF5PEdlb0pTT05GZWF0dXJlPjtcbn1cblxuZnVuY3Rpb24gbG9hZEdlb0pTT05UaWxlKHBhcmFtczogUmVxdWVzdGVkVGlsZVBhcmFtZXRlcnMsIGNhbGxiYWNrOiBMb2FkVmVjdG9yRGF0YUNhbGxiYWNrKSB7XG4gICAgY29uc3QgY2Fub25pY2FsID0gcGFyYW1zLnRpbGVJRC5jYW5vbmljYWw7XG5cbiAgICBpZiAoIXRoaXMuX2dlb0pTT05JbmRleCkge1xuICAgICAgICByZXR1cm4gY2FsbGJhY2sobnVsbCwgbnVsbCk7ICAvLyB3ZSBjb3VsZG4ndCBsb2FkIHRoZSBmaWxlXG4gICAgfVxuXG4gICAgY29uc3QgZ2VvSlNPTlRpbGUgPSB0aGlzLl9nZW9KU09OSW5kZXguZ2V0VGlsZShjYW5vbmljYWwueiwgY2Fub25pY2FsLngsIGNhbm9uaWNhbC55KTtcbiAgICBpZiAoIWdlb0pTT05UaWxlKSB7XG4gICAgICAgIHJldHVybiBjYWxsYmFjayhudWxsLCBudWxsKTsgLy8gbm90aGluZyBpbiB0aGUgZ2l2ZW4gdGlsZVxuICAgIH1cblxuICAgIGNvbnN0IGdlb2pzb25XcmFwcGVyID0gbmV3IEdlb0pTT05XcmFwcGVyKGdlb0pTT05UaWxlLmZlYXR1cmVzKTtcblxuICAgIC8vIEVuY29kZSB0aGUgZ2VvanNvbi12dCB0aWxlIGludG8gYmluYXJ5IHZlY3RvciB0aWxlIGZvcm0uICBUaGlzXG4gICAgLy8gaXMgYSBjb252ZW5pZW5jZSB0aGF0IGFsbG93cyBgRmVhdHVyZUluZGV4YCB0byBvcGVyYXRlIHRoZSBzYW1lIHdheVxuICAgIC8vIGFjcm9zcyBgVmVjdG9yVGlsZVNvdXJjZWAgYW5kIGBHZW9KU09OU291cmNlYCBkYXRhLlxuICAgIGxldCBwYmYgPSB2dHBiZihnZW9qc29uV3JhcHBlcik7XG4gICAgaWYgKHBiZi5ieXRlT2Zmc2V0ICE9PSAwIHx8IHBiZi5ieXRlTGVuZ3RoICE9PSBwYmYuYnVmZmVyLmJ5dGVMZW5ndGgpIHtcbiAgICAgICAgLy8gQ29tcGF0aWJpbGl0eSB3aXRoIG5vZGUgQnVmZmVyIChodHRwczovL2dpdGh1Yi5jb20vbWFwYm94L3BiZi9pc3N1ZXMvMzUpXG4gICAgICAgIHBiZiA9IG5ldyBVaW50OEFycmF5KHBiZik7XG4gICAgfVxuXG4gICAgY2FsbGJhY2sobnVsbCwge1xuICAgICAgICB2ZWN0b3JUaWxlOiBnZW9qc29uV3JhcHBlcixcbiAgICAgICAgcmF3RGF0YTogcGJmLmJ1ZmZlclxuICAgIH0pO1xufVxuXG4vKipcbiAqIFRoZSB7QGxpbmsgV29ya2VyU291cmNlfSBpbXBsZW1lbnRhdGlvbiB0aGF0IHN1cHBvcnRzIHtAbGluayBHZW9KU09OU291cmNlfS5cbiAqIFRoaXMgY2xhc3MgaXMgZGVzaWduZWQgdG8gYmUgZWFzaWx5IHJldXNlZCB0byBzdXBwb3J0IGN1c3RvbSBzb3VyY2UgdHlwZXNcbiAqIGZvciBkYXRhIGZvcm1hdHMgdGhhdCBjYW4gYmUgcGFyc2VkL2NvbnZlcnRlZCBpbnRvIGFuIGluLW1lbW9yeSBHZW9KU09OXG4gKiByZXByZXNlbnRhdGlvbi4gIFRvIGRvIHNvLCBjcmVhdGUgaXQgd2l0aFxuICogYG5ldyBHZW9KU09OV29ya2VyU291cmNlKGFjdG9yLCBsYXllckluZGV4LCBjdXN0b21Mb2FkR2VvSlNPTkZ1bmN0aW9uKWAuXG4gKiBGb3IgYSBmdWxsIGV4YW1wbGUsIHNlZSBbbWFwYm94LWdsLXRvcG9qc29uXShodHRwczovL2dpdGh1Yi5jb20vZGV2ZWxvcG1lbnRzZWVkL21hcGJveC1nbC10b3BvanNvbikuXG4gKlxuICogQHByaXZhdGVcbiAqL1xuY2xhc3MgR2VvSlNPTldvcmtlclNvdXJjZSBleHRlbmRzIFZlY3RvclRpbGVXb3JrZXJTb3VyY2Uge1xuICAgIGxvYWRHZW9KU09OOiBMb2FkR2VvSlNPTjtcbiAgICBfZ2VvSlNPTkluZGV4OiBHZW9KU09OSW5kZXhcblxuICAgIC8qKlxuICAgICAqIEBwYXJhbSBbbG9hZEdlb0pTT05dIE9wdGlvbmFsIG1ldGhvZCBmb3IgY3VzdG9tIGxvYWRpbmcvcGFyc2luZyBvZlxuICAgICAqIEdlb0pTT04gYmFzZWQgb24gcGFyYW1ldGVycyBwYXNzZWQgZnJvbSB0aGUgbWFpbi10aHJlYWQgU291cmNlLlxuICAgICAqIFNlZSB7QGxpbmsgR2VvSlNPTldvcmtlclNvdXJjZSNsb2FkR2VvSlNPTn0uXG4gICAgICogQHByaXZhdGVcbiAgICAgKi9cbiAgICBjb25zdHJ1Y3RvcihhY3RvcjogQWN0b3IsIGxheWVySW5kZXg6IFN0eWxlTGF5ZXJJbmRleCwgYXZhaWxhYmxlSW1hZ2VzOiBBcnJheTxzdHJpbmc+LCBpc1Nwcml0ZUxvYWRlZDogYm9vbGVhbiwgbG9hZEdlb0pTT046ID9Mb2FkR2VvSlNPTikge1xuICAgICAgICBzdXBlcihhY3RvciwgbGF5ZXJJbmRleCwgYXZhaWxhYmxlSW1hZ2VzLCBpc1Nwcml0ZUxvYWRlZCwgbG9hZEdlb0pTT05UaWxlKTtcbiAgICAgICAgaWYgKGxvYWRHZW9KU09OKSB7XG4gICAgICAgICAgICB0aGlzLmxvYWRHZW9KU09OID0gbG9hZEdlb0pTT047XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBGZXRjaGVzIChpZiBhcHByb3ByaWF0ZSksIHBhcnNlcywgYW5kIGluZGV4IGdlb2pzb24gZGF0YSBpbnRvIHRpbGVzLiBUaGlzXG4gICAgICogcHJlcGFyYXRvcnkgbWV0aG9kIG11c3QgYmUgY2FsbGVkIGJlZm9yZSB7QGxpbmsgR2VvSlNPTldvcmtlclNvdXJjZSNsb2FkVGlsZX1cbiAgICAgKiBjYW4gY29ycmVjdGx5IHNlcnZlIHVwIHRpbGVzLlxuICAgICAqXG4gICAgICogRGVmZXJzIHRvIHtAbGluayBHZW9KU09OV29ya2VyU291cmNlI2xvYWRHZW9KU09OfSBmb3IgdGhlIGZldGNoaW5nL3BhcnNpbmcsXG4gICAgICogZXhwZWN0aW5nIGBjYWxsYmFjayhlcnJvciwgZGF0YSlgIHRvIGJlIGNhbGxlZCB3aXRoIGVpdGhlciBhbiBlcnJvciBvciBhXG4gICAgICogcGFyc2VkIEdlb0pTT04gb2JqZWN0LlxuICAgICAqXG4gICAgICogV2hlbiBgbG9hZERhdGFgIHJlcXVlc3RzIGNvbWUgaW4gZmFzdGVyIHRoYW4gdGhleSBjYW4gYmUgcHJvY2Vzc2VkLFxuICAgICAqIHRoZXkgYXJlIGNvYWxlc2NlZCBpbnRvIGEgc2luZ2xlIHJlcXVlc3QgdXNpbmcgdGhlIGxhdGVzdCBkYXRhLlxuICAgICAqIFNlZSB7QGxpbmsgR2VvSlNPTldvcmtlclNvdXJjZSNjb2FsZXNjZX1cbiAgICAgKlxuICAgICAqIEBwYXJhbSBwYXJhbXNcbiAgICAgKiBAcGFyYW0gY2FsbGJhY2tcbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqL1xuICAgIGxvYWREYXRhKHBhcmFtczogTG9hZEdlb0pTT05QYXJhbWV0ZXJzLCBjYWxsYmFjazogQ2FsbGJhY2s8e3Jlc291cmNlVGltaW5nPzoge1tfOiBzdHJpbmddOiBBcnJheTxQZXJmb3JtYW5jZVJlc291cmNlVGltaW5nPn19Pikge1xuICAgICAgICBjb25zdCByZXF1ZXN0UGFyYW0gPSBwYXJhbXMgJiYgcGFyYW1zLnJlcXVlc3Q7XG4gICAgICAgIGNvbnN0IHBlcmYgPSByZXF1ZXN0UGFyYW0gJiYgcmVxdWVzdFBhcmFtLmNvbGxlY3RSZXNvdXJjZVRpbWluZztcblxuICAgICAgICB0aGlzLmxvYWRHZW9KU09OKHBhcmFtcywgKGVycjogP0Vycm9yLCBkYXRhOiA/T2JqZWN0KSA9PiB7XG4gICAgICAgICAgICBpZiAoZXJyIHx8ICFkYXRhKSB7XG4gICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKGVycik7XG4gICAgICAgICAgICB9IGVsc2UgaWYgKHR5cGVvZiBkYXRhICE9PSAnb2JqZWN0Jykge1xuICAgICAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhuZXcgRXJyb3IoYElucHV0IGRhdGEgZ2l2ZW4gdG8gJyR7cGFyYW1zLnNvdXJjZX0nIGlzIG5vdCBhIHZhbGlkIEdlb0pTT04gb2JqZWN0LmApKTtcbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgcmV3aW5kKGRhdGEsIHRydWUpO1xuXG4gICAgICAgICAgICAgICAgdHJ5IHtcbiAgICAgICAgICAgICAgICAgICAgaWYgKHBhcmFtcy5maWx0ZXIpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIGNvbnN0IGNvbXBpbGVkID0gY3JlYXRlRXhwcmVzc2lvbihwYXJhbXMuZmlsdGVyLCB7dHlwZTogJ2Jvb2xlYW4nLCAncHJvcGVydHktdHlwZSc6ICdkYXRhLWRyaXZlbicsIG92ZXJyaWRhYmxlOiBmYWxzZSwgdHJhbnNpdGlvbjogZmFsc2V9KTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGlmIChjb21waWxlZC5yZXN1bHQgPT09ICdlcnJvcicpXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdGhyb3cgbmV3IEVycm9yKGNvbXBpbGVkLnZhbHVlLm1hcChlcnIgPT4gYCR7ZXJyLmtleX06ICR7ZXJyLm1lc3NhZ2V9YCkuam9pbignLCAnKSk7XG5cbiAgICAgICAgICAgICAgICAgICAgICAgIGNvbnN0IGZlYXR1cmVzID0gZGF0YS5mZWF0dXJlcy5maWx0ZXIoZmVhdHVyZSA9PiBjb21waWxlZC52YWx1ZS5ldmFsdWF0ZSh7em9vbTogMH0sIGZlYXR1cmUpKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGRhdGEgPSB7dHlwZTogJ0ZlYXR1cmVDb2xsZWN0aW9uJywgZmVhdHVyZXN9O1xuICAgICAgICAgICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICAgICAgdGhpcy5fZ2VvSlNPTkluZGV4ID0gcGFyYW1zLmNsdXN0ZXIgP1xuICAgICAgICAgICAgICAgICAgICAgICAgbmV3IFN1cGVyY2x1c3RlcihnZXRTdXBlcmNsdXN0ZXJPcHRpb25zKHBhcmFtcykpLmxvYWQoZGF0YS5mZWF0dXJlcykgOlxuICAgICAgICAgICAgICAgICAgICAgICAgZ2VvanNvbnZ0KGRhdGEsIHBhcmFtcy5nZW9qc29uVnRPcHRpb25zKTtcbiAgICAgICAgICAgICAgICB9IGNhdGNoIChlcnIpIHtcbiAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKGVycik7XG4gICAgICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICAgICAgdGhpcy5sb2FkZWQgPSB7fTtcblxuICAgICAgICAgICAgICAgIGNvbnN0IHJlc3VsdCA9IHt9O1xuICAgICAgICAgICAgICAgIGlmIChwZXJmKSB7XG4gICAgICAgICAgICAgICAgICAgIGNvbnN0IHJlc291cmNlVGltaW5nRGF0YSA9IGdldFBlcmZvcm1hbmNlTWVhc3VyZW1lbnQocmVxdWVzdFBhcmFtKTtcbiAgICAgICAgICAgICAgICAgICAgLy8gaXQncyBuZWNlc3NhcnkgdG8gZXZhbCB0aGUgcmVzdWx0IG9mIGdldEVudHJpZXNCeU5hbWUoKSBoZXJlIHZpYSBwYXJzZS9zdHJpbmdpZnlcbiAgICAgICAgICAgICAgICAgICAgLy8gbGF0ZSBldmFsdWF0aW9uIGluIHRoZSBtYWluIHRocmVhZCBjYXVzZXMgVHlwZUVycm9yOiBpbGxlZ2FsIGludm9jYXRpb25cbiAgICAgICAgICAgICAgICAgICAgaWYgKHJlc291cmNlVGltaW5nRGF0YSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgcmVzdWx0LnJlc291cmNlVGltaW5nID0ge307XG4gICAgICAgICAgICAgICAgICAgICAgICByZXN1bHQucmVzb3VyY2VUaW1pbmdbcGFyYW1zLnNvdXJjZV0gPSBKU09OLnBhcnNlKEpTT04uc3RyaW5naWZ5KHJlc291cmNlVGltaW5nRGF0YSkpO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIGNhbGxiYWNrKG51bGwsIHJlc3VsdCk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH0pO1xuICAgIH1cblxuICAgIC8qKlxuICAgICogSW1wbGVtZW50cyB7QGxpbmsgV29ya2VyU291cmNlI3JlbG9hZFRpbGV9LlxuICAgICpcbiAgICAqIElmIHRoZSB0aWxlIGlzIGxvYWRlZCwgdXNlcyB0aGUgaW1wbGVtZW50YXRpb24gaW4gVmVjdG9yVGlsZVdvcmtlclNvdXJjZS5cbiAgICAqIE90aGVyd2lzZSwgc3VjaCBhcyBhZnRlciBhIHNldERhdGEoKSBjYWxsLCB3ZSBsb2FkIHRoZSB0aWxlIGZyZXNoLlxuICAgICpcbiAgICAqIEBwYXJhbSBwYXJhbXNcbiAgICAqIEBwYXJhbSBwYXJhbXMudWlkIFRoZSBVSUQgZm9yIHRoaXMgdGlsZS5cbiAgICAqIEBwcml2YXRlXG4gICAgKi9cbiAgICByZWxvYWRUaWxlKHBhcmFtczogV29ya2VyVGlsZVBhcmFtZXRlcnMsIGNhbGxiYWNrOiBXb3JrZXJUaWxlQ2FsbGJhY2spOiB2b2lkIHtcbiAgICAgICAgY29uc3QgbG9hZGVkID0gdGhpcy5sb2FkZWQsXG4gICAgICAgICAgICB1aWQgPSBwYXJhbXMudWlkO1xuXG4gICAgICAgIGlmIChsb2FkZWQgJiYgbG9hZGVkW3VpZF0pIHtcbiAgICAgICAgICAgIHJldHVybiBzdXBlci5yZWxvYWRUaWxlKHBhcmFtcywgY2FsbGJhY2spO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgcmV0dXJuIHRoaXMubG9hZFRpbGUocGFyYW1zLCBjYWxsYmFjayk7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBGZXRjaCBhbmQgcGFyc2UgR2VvSlNPTiBhY2NvcmRpbmcgdG8gdGhlIGdpdmVuIHBhcmFtcy4gIENhbGxzIGBjYWxsYmFja2BcbiAgICAgKiB3aXRoIGAoZXJyLCBkYXRhKWAsIHdoZXJlIGBkYXRhYCBpcyBhIHBhcnNlZCBHZW9KU09OIG9iamVjdC5cbiAgICAgKlxuICAgICAqIEdlb0pTT04gaXMgbG9hZGVkIGFuZCBwYXJzZWQgZnJvbSBgcGFyYW1zLnVybGAgaWYgaXQgZXhpc3RzLCBvciBlbHNlXG4gICAgICogZXhwZWN0ZWQgYXMgYSBsaXRlcmFsIChzdHJpbmcgb3Igb2JqZWN0KSBgcGFyYW1zLmRhdGFgLlxuICAgICAqXG4gICAgICogQHBhcmFtIHBhcmFtc1xuICAgICAqIEBwYXJhbSBbcGFyYW1zLnVybF0gQSBVUkwgdG8gdGhlIHJlbW90ZSBHZW9KU09OIGRhdGEuXG4gICAgICogQHBhcmFtIFtwYXJhbXMuZGF0YV0gTGl0ZXJhbCBHZW9KU09OIGRhdGEuIE11c3QgYmUgcHJvdmlkZWQgaWYgYHBhcmFtcy51cmxgIGlzIG5vdC5cbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqL1xuICAgIGxvYWRHZW9KU09OKHBhcmFtczogTG9hZEdlb0pTT05QYXJhbWV0ZXJzLCBjYWxsYmFjazogUmVzcG9uc2VDYWxsYmFjazxPYmplY3Q+KTogdm9pZCB7XG4gICAgICAgIC8vIEJlY2F1c2Ugb2Ygc2FtZSBvcmlnaW4gaXNzdWVzLCB1cmxzIG11c3QgZWl0aGVyIGluY2x1ZGUgYW4gZXhwbGljaXRcbiAgICAgICAgLy8gb3JpZ2luIG9yIGFic29sdXRlIHBhdGguXG4gICAgICAgIC8vIGllOiAvZm9vL2Jhci5qc29uIG9yIGh0dHA6Ly9leGFtcGxlLmNvbS9iYXIuanNvblxuICAgICAgICAvLyBidXQgbm90IC4uL2Zvby9iYXIuanNvblxuICAgICAgICBpZiAocGFyYW1zLnJlcXVlc3QpIHtcbiAgICAgICAgICAgIGdldEpTT04ocGFyYW1zLnJlcXVlc3QsIGNhbGxiYWNrKTtcbiAgICAgICAgfSBlbHNlIGlmICh0eXBlb2YgcGFyYW1zLmRhdGEgPT09ICdzdHJpbmcnKSB7XG4gICAgICAgICAgICB0cnkge1xuICAgICAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhudWxsLCBKU09OLnBhcnNlKHBhcmFtcy5kYXRhKSk7XG4gICAgICAgICAgICB9IGNhdGNoIChlKSB7XG4gICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKG5ldyBFcnJvcihgSW5wdXQgZGF0YSBnaXZlbiB0byAnJHtwYXJhbXMuc291cmNlfScgaXMgbm90IGEgdmFsaWQgR2VvSlNPTiBvYmplY3QuYCkpO1xuICAgICAgICAgICAgfVxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKG5ldyBFcnJvcihgSW5wdXQgZGF0YSBnaXZlbiB0byAnJHtwYXJhbXMuc291cmNlfScgaXMgbm90IGEgdmFsaWQgR2VvSlNPTiBvYmplY3QuYCkpO1xuICAgICAgICB9XG4gICAgfVxuXG4gICAgZ2V0Q2x1c3RlckV4cGFuc2lvblpvb20ocGFyYW1zOiB7Y2x1c3RlcklkOiBudW1iZXJ9LCBjYWxsYmFjazogQ2FsbGJhY2s8bnVtYmVyPikge1xuICAgICAgICB0cnkge1xuICAgICAgICAgICAgY2FsbGJhY2sobnVsbCwgdGhpcy5fZ2VvSlNPTkluZGV4LmdldENsdXN0ZXJFeHBhbnNpb25ab29tKHBhcmFtcy5jbHVzdGVySWQpKTtcbiAgICAgICAgfSBjYXRjaCAoZSkge1xuICAgICAgICAgICAgY2FsbGJhY2soZSk7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICBnZXRDbHVzdGVyQ2hpbGRyZW4ocGFyYW1zOiB7Y2x1c3RlcklkOiBudW1iZXJ9LCBjYWxsYmFjazogQ2FsbGJhY2s8QXJyYXk8R2VvSlNPTkZlYXR1cmU+Pikge1xuICAgICAgICB0cnkge1xuICAgICAgICAgICAgY2FsbGJhY2sobnVsbCwgdGhpcy5fZ2VvSlNPTkluZGV4LmdldENoaWxkcmVuKHBhcmFtcy5jbHVzdGVySWQpKTtcbiAgICAgICAgfSBjYXRjaCAoZSkge1xuICAgICAgICAgICAgY2FsbGJhY2soZSk7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICBnZXRDbHVzdGVyTGVhdmVzKHBhcmFtczoge2NsdXN0ZXJJZDogbnVtYmVyLCBsaW1pdDogbnVtYmVyLCBvZmZzZXQ6IG51bWJlcn0sIGNhbGxiYWNrOiBDYWxsYmFjazxBcnJheTxHZW9KU09ORmVhdHVyZT4+KSB7XG4gICAgICAgIHRyeSB7XG4gICAgICAgICAgICBjYWxsYmFjayhudWxsLCB0aGlzLl9nZW9KU09OSW5kZXguZ2V0TGVhdmVzKHBhcmFtcy5jbHVzdGVySWQsIHBhcmFtcy5saW1pdCwgcGFyYW1zLm9mZnNldCkpO1xuICAgICAgICB9IGNhdGNoIChlKSB7XG4gICAgICAgICAgICBjYWxsYmFjayhlKTtcbiAgICAgICAgfVxuICAgIH1cbn1cblxuZnVuY3Rpb24gZ2V0U3VwZXJjbHVzdGVyT3B0aW9ucyh7c3VwZXJjbHVzdGVyT3B0aW9ucywgY2x1c3RlclByb3BlcnRpZXN9KSB7XG4gICAgaWYgKCFjbHVzdGVyUHJvcGVydGllcyB8fCAhc3VwZXJjbHVzdGVyT3B0aW9ucykgcmV0dXJuIHN1cGVyY2x1c3Rlck9wdGlvbnM7XG5cbiAgICBjb25zdCBtYXBFeHByZXNzaW9ucyA9IHt9O1xuICAgIGNvbnN0IHJlZHVjZUV4cHJlc3Npb25zID0ge307XG4gICAgY29uc3QgZ2xvYmFscyA9IHthY2N1bXVsYXRlZDogbnVsbCwgem9vbTogMH07XG4gICAgY29uc3QgZmVhdHVyZSA9IHtwcm9wZXJ0aWVzOiBudWxsfTtcbiAgICBjb25zdCBwcm9wZXJ0eU5hbWVzID0gT2JqZWN0LmtleXMoY2x1c3RlclByb3BlcnRpZXMpO1xuXG4gICAgZm9yIChjb25zdCBrZXkgb2YgcHJvcGVydHlOYW1lcykge1xuICAgICAgICBjb25zdCBbb3BlcmF0b3IsIG1hcEV4cHJlc3Npb25dID0gY2x1c3RlclByb3BlcnRpZXNba2V5XTtcblxuICAgICAgICBjb25zdCBtYXBFeHByZXNzaW9uUGFyc2VkID0gY3JlYXRlRXhwcmVzc2lvbihtYXBFeHByZXNzaW9uKTtcbiAgICAgICAgY29uc3QgcmVkdWNlRXhwcmVzc2lvblBhcnNlZCA9IGNyZWF0ZUV4cHJlc3Npb24oXG4gICAgICAgICAgICB0eXBlb2Ygb3BlcmF0b3IgPT09ICdzdHJpbmcnID8gW29wZXJhdG9yLCBbJ2FjY3VtdWxhdGVkJ10sIFsnZ2V0Jywga2V5XV0gOiBvcGVyYXRvcik7XG5cbiAgICAgICAgYXNzZXJ0KG1hcEV4cHJlc3Npb25QYXJzZWQucmVzdWx0ID09PSAnc3VjY2VzcycpO1xuICAgICAgICBhc3NlcnQocmVkdWNlRXhwcmVzc2lvblBhcnNlZC5yZXN1bHQgPT09ICdzdWNjZXNzJyk7XG5cbiAgICAgICAgbWFwRXhwcmVzc2lvbnNba2V5XSA9IG1hcEV4cHJlc3Npb25QYXJzZWQudmFsdWU7XG4gICAgICAgIHJlZHVjZUV4cHJlc3Npb25zW2tleV0gPSByZWR1Y2VFeHByZXNzaW9uUGFyc2VkLnZhbHVlO1xuICAgIH1cblxuICAgIHN1cGVyY2x1c3Rlck9wdGlvbnMubWFwID0gKHBvaW50UHJvcGVydGllcykgPT4ge1xuICAgICAgICBmZWF0dXJlLnByb3BlcnRpZXMgPSBwb2ludFByb3BlcnRpZXM7XG4gICAgICAgIGNvbnN0IHByb3BlcnRpZXMgPSB7fTtcbiAgICAgICAgZm9yIChjb25zdCBrZXkgb2YgcHJvcGVydHlOYW1lcykge1xuICAgICAgICAgICAgcHJvcGVydGllc1trZXldID0gbWFwRXhwcmVzc2lvbnNba2V5XS5ldmFsdWF0ZShnbG9iYWxzLCBmZWF0dXJlKTtcbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gcHJvcGVydGllcztcbiAgICB9O1xuICAgIHN1cGVyY2x1c3Rlck9wdGlvbnMucmVkdWNlID0gKGFjY3VtdWxhdGVkLCBjbHVzdGVyUHJvcGVydGllcykgPT4ge1xuICAgICAgICBmZWF0dXJlLnByb3BlcnRpZXMgPSBjbHVzdGVyUHJvcGVydGllcztcbiAgICAgICAgZm9yIChjb25zdCBrZXkgb2YgcHJvcGVydHlOYW1lcykge1xuICAgICAgICAgICAgZ2xvYmFscy5hY2N1bXVsYXRlZCA9IGFjY3VtdWxhdGVkW2tleV07XG4gICAgICAgICAgICBhY2N1bXVsYXRlZFtrZXldID0gcmVkdWNlRXhwcmVzc2lvbnNba2V5XS5ldmFsdWF0ZShnbG9iYWxzLCBmZWF0dXJlKTtcbiAgICAgICAgfVxuICAgIH07XG5cbiAgICByZXR1cm4gc3VwZXJjbHVzdGVyT3B0aW9ucztcbn1cblxuZXhwb3J0IGRlZmF1bHQgR2VvSlNPTldvcmtlclNvdXJjZTtcbiIsIi8vIEBmbG93XG5cbmltcG9ydCBBY3RvciBmcm9tICcuLi91dGlsL2FjdG9yLmpzJztcblxuaW1wb3J0IFN0eWxlTGF5ZXJJbmRleCBmcm9tICcuLi9zdHlsZS9zdHlsZV9sYXllcl9pbmRleC5qcyc7XG5pbXBvcnQgVmVjdG9yVGlsZVdvcmtlclNvdXJjZSBmcm9tICcuL3ZlY3Rvcl90aWxlX3dvcmtlcl9zb3VyY2UuanMnO1xuaW1wb3J0IFJhc3RlckRFTVRpbGVXb3JrZXJTb3VyY2UgZnJvbSAnLi9yYXN0ZXJfZGVtX3RpbGVfd29ya2VyX3NvdXJjZS5qcyc7XG5pbXBvcnQgR2VvSlNPTldvcmtlclNvdXJjZSBmcm9tICcuL2dlb2pzb25fd29ya2VyX3NvdXJjZS5qcyc7XG5pbXBvcnQgYXNzZXJ0IGZyb20gJ2Fzc2VydCc7XG5pbXBvcnQge3BsdWdpbiBhcyBnbG9iYWxSVExUZXh0UGx1Z2lufSBmcm9tICcuL3J0bF90ZXh0X3BsdWdpbi5qcyc7XG5pbXBvcnQge2VuZm9yY2VDYWNoZVNpemVMaW1pdH0gZnJvbSAnLi4vdXRpbC90aWxlX3JlcXVlc3RfY2FjaGUuanMnO1xuaW1wb3J0IHtleHRlbmR9IGZyb20gJy4uL3V0aWwvdXRpbC5qcyc7XG5pbXBvcnQge1BlcmZvcm1hbmNlVXRpbHN9IGZyb20gJy4uL3V0aWwvcGVyZm9ybWFuY2UuanMnO1xuaW1wb3J0IHtFdmVudH0gZnJvbSAnLi4vdXRpbC9ldmVudGVkLmpzJztcbmltcG9ydCB7Z2V0UHJvamVjdGlvbn0gZnJvbSAnLi4vZ2VvL3Byb2plY3Rpb24vaW5kZXguanMnO1xuXG5pbXBvcnQgdHlwZSB7XG4gICAgV29ya2VyU291cmNlLFxuICAgIFdvcmtlclRpbGVQYXJhbWV0ZXJzLFxuICAgIFdvcmtlckRFTVRpbGVQYXJhbWV0ZXJzLFxuICAgIFdvcmtlclRpbGVDYWxsYmFjayxcbiAgICBXb3JrZXJERU1UaWxlQ2FsbGJhY2ssXG4gICAgVGlsZVBhcmFtZXRlcnNcbn0gZnJvbSAnLi4vc291cmNlL3dvcmtlcl9zb3VyY2UuanMnO1xuXG5pbXBvcnQgdHlwZSB7V29ya2VyR2xvYmFsU2NvcGVJbnRlcmZhY2V9IGZyb20gJy4uL3V0aWwvd2ViX3dvcmtlci5qcyc7XG5pbXBvcnQgdHlwZSB7Q2FsbGJhY2t9IGZyb20gJy4uL3R5cGVzL2NhbGxiYWNrLmpzJztcbmltcG9ydCB0eXBlIHtMYXllclNwZWNpZmljYXRpb24sIFByb2plY3Rpb25TcGVjaWZpY2F0aW9ufSBmcm9tICcuLi9zdHlsZS1zcGVjL3R5cGVzLmpzJztcbmltcG9ydCB0eXBlIHtQbHVnaW5TdGF0ZX0gZnJvbSAnLi9ydGxfdGV4dF9wbHVnaW4uanMnO1xuaW1wb3J0IHR5cGUgUHJvamVjdGlvbiBmcm9tICcuLi9nZW8vcHJvamVjdGlvbi9wcm9qZWN0aW9uLmpzJztcblxuLyoqXG4gKiBAcHJpdmF0ZVxuICovXG5leHBvcnQgZGVmYXVsdCBjbGFzcyBXb3JrZXIge1xuICAgIHNlbGY6IFdvcmtlckdsb2JhbFNjb3BlSW50ZXJmYWNlO1xuICAgIGFjdG9yOiBBY3RvcjtcbiAgICBsYXllckluZGV4ZXM6IHtbXzogc3RyaW5nXTogU3R5bGVMYXllckluZGV4IH07XG4gICAgYXZhaWxhYmxlSW1hZ2VzOiB7W186IHN0cmluZ106IEFycmF5PHN0cmluZz4gfTtcbiAgICB3b3JrZXJTb3VyY2VUeXBlczoge1tfOiBzdHJpbmddOiBDbGFzczxXb3JrZXJTb3VyY2U+IH07XG4gICAgd29ya2VyU291cmNlczoge1tfOiBzdHJpbmddOiB7W186IHN0cmluZ106IHtbXzogc3RyaW5nXTogV29ya2VyU291cmNlIH0gfSB9O1xuICAgIGRlbVdvcmtlclNvdXJjZXM6IHtbXzogc3RyaW5nXToge1tfOiBzdHJpbmddOiBSYXN0ZXJERU1UaWxlV29ya2VyU291cmNlIH0gfTtcbiAgICBwcm9qZWN0aW9uczoge1tfOiBzdHJpbmddOiBQcm9qZWN0aW9uIH07XG4gICAgZGVmYXVsdFByb2plY3Rpb246IFByb2plY3Rpb247XG4gICAgaXNTcHJpdGVMb2FkZWQ6IHtbXzogc3RyaW5nXTogYm9vbGVhbiB9O1xuICAgIHJlZmVycmVyOiA/c3RyaW5nO1xuICAgIHRlcnJhaW46ID9ib29sZWFuO1xuXG4gICAgY29uc3RydWN0b3Ioc2VsZjogV29ya2VyR2xvYmFsU2NvcGVJbnRlcmZhY2UpIHtcbiAgICAgICAgUGVyZm9ybWFuY2VVdGlscy5tZWFzdXJlKCd3b3JrZXJFdmFsdWF0ZVNjcmlwdCcpO1xuICAgICAgICB0aGlzLnNlbGYgPSBzZWxmO1xuICAgICAgICB0aGlzLmFjdG9yID0gbmV3IEFjdG9yKHNlbGYsIHRoaXMpO1xuXG4gICAgICAgIHRoaXMubGF5ZXJJbmRleGVzID0ge307XG4gICAgICAgIHRoaXMuYXZhaWxhYmxlSW1hZ2VzID0ge307XG4gICAgICAgIHRoaXMuaXNTcHJpdGVMb2FkZWQgPSB7fTtcblxuICAgICAgICB0aGlzLnByb2plY3Rpb25zID0ge307XG4gICAgICAgIHRoaXMuZGVmYXVsdFByb2plY3Rpb24gPSBnZXRQcm9qZWN0aW9uKHtuYW1lOiAnbWVyY2F0b3InfSk7XG5cbiAgICAgICAgdGhpcy53b3JrZXJTb3VyY2VUeXBlcyA9IHtcbiAgICAgICAgICAgIHZlY3RvcjogVmVjdG9yVGlsZVdvcmtlclNvdXJjZSxcbiAgICAgICAgICAgIGdlb2pzb246IEdlb0pTT05Xb3JrZXJTb3VyY2VcbiAgICAgICAgfTtcblxuICAgICAgICAvLyBbbWFwSWRdW3NvdXJjZVR5cGVdW3NvdXJjZU5hbWVdID0+IHdvcmtlciBzb3VyY2UgaW5zdGFuY2VcbiAgICAgICAgdGhpcy53b3JrZXJTb3VyY2VzID0ge307XG4gICAgICAgIHRoaXMuZGVtV29ya2VyU291cmNlcyA9IHt9O1xuXG4gICAgICAgIHRoaXMuc2VsZi5yZWdpc3RlcldvcmtlclNvdXJjZSA9IChuYW1lOiBzdHJpbmcsIFdvcmtlclNvdXJjZTogQ2xhc3M8V29ya2VyU291cmNlPikgPT4ge1xuICAgICAgICAgICAgaWYgKHRoaXMud29ya2VyU291cmNlVHlwZXNbbmFtZV0pIHtcbiAgICAgICAgICAgICAgICB0aHJvdyBuZXcgRXJyb3IoYFdvcmtlciBzb3VyY2Ugd2l0aCBuYW1lIFwiJHtuYW1lfVwiIGFscmVhZHkgcmVnaXN0ZXJlZC5gKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHRoaXMud29ya2VyU291cmNlVHlwZXNbbmFtZV0gPSBXb3JrZXJTb3VyY2U7XG4gICAgICAgIH07XG5cbiAgICAgICAgLy8gVGhpcyBpcyBpbnZva2VkIGJ5IHRoZSBSVEwgdGV4dCBwbHVnaW4gd2hlbiB0aGUgZG93bmxvYWQgdmlhIHRoZSBgaW1wb3J0U2NyaXB0c2AgY2FsbCBoYXMgZmluaXNoZWQsIGFuZCB0aGUgY29kZSBoYXMgYmVlbiBwYXJzZWQuXG4gICAgICAgIHRoaXMuc2VsZi5yZWdpc3RlclJUTFRleHRQbHVnaW4gPSAocnRsVGV4dFBsdWdpbjoge2FwcGx5QXJhYmljU2hhcGluZzogRnVuY3Rpb24sIHByb2Nlc3NCaWRpcmVjdGlvbmFsVGV4dDogRnVuY3Rpb24sIHByb2Nlc3NTdHlsZWRCaWRpcmVjdGlvbmFsVGV4dD86IEZ1bmN0aW9ufSkgPT4ge1xuICAgICAgICAgICAgaWYgKGdsb2JhbFJUTFRleHRQbHVnaW4uaXNQYXJzZWQoKSkge1xuICAgICAgICAgICAgICAgIHRocm93IG5ldyBFcnJvcignUlRMIHRleHQgcGx1Z2luIGFscmVhZHkgcmVnaXN0ZXJlZC4nKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGdsb2JhbFJUTFRleHRQbHVnaW5bJ2FwcGx5QXJhYmljU2hhcGluZyddID0gcnRsVGV4dFBsdWdpbi5hcHBseUFyYWJpY1NoYXBpbmc7XG4gICAgICAgICAgICBnbG9iYWxSVExUZXh0UGx1Z2luWydwcm9jZXNzQmlkaXJlY3Rpb25hbFRleHQnXSA9IHJ0bFRleHRQbHVnaW4ucHJvY2Vzc0JpZGlyZWN0aW9uYWxUZXh0O1xuICAgICAgICAgICAgZ2xvYmFsUlRMVGV4dFBsdWdpblsncHJvY2Vzc1N0eWxlZEJpZGlyZWN0aW9uYWxUZXh0J10gPSBydGxUZXh0UGx1Z2luLnByb2Nlc3NTdHlsZWRCaWRpcmVjdGlvbmFsVGV4dDtcbiAgICAgICAgfTtcbiAgICB9XG5cbiAgICBjbGVhckNhY2hlcyhtYXBJZDogc3RyaW5nLCB1bnVzZWQ6IG1peGVkLCBjYWxsYmFjazogV29ya2VyVGlsZUNhbGxiYWNrKSB7XG4gICAgICAgIGRlbGV0ZSB0aGlzLmxheWVySW5kZXhlc1ttYXBJZF07XG4gICAgICAgIGRlbGV0ZSB0aGlzLmF2YWlsYWJsZUltYWdlc1ttYXBJZF07XG4gICAgICAgIGRlbGV0ZSB0aGlzLndvcmtlclNvdXJjZXNbbWFwSWRdO1xuICAgICAgICBkZWxldGUgdGhpcy5kZW1Xb3JrZXJTb3VyY2VzW21hcElkXTtcbiAgICAgICAgY2FsbGJhY2soKTtcbiAgICB9XG5cbiAgICBjaGVja0lmUmVhZHkobWFwSUQ6IHN0cmluZywgdW51c2VkOiBtaXhlZCwgY2FsbGJhY2s6IFdvcmtlclRpbGVDYWxsYmFjaykge1xuICAgICAgICAvLyBub29wLCB1c2VkIHRvIGNoZWNrIGlmIGEgd29ya2VyIGlzIGZ1bGx5IHNldCB1cCBhbmQgcmVhZHkgdG8gcmVjZWl2ZSBtZXNzYWdlc1xuICAgICAgICBjYWxsYmFjaygpO1xuICAgIH1cblxuICAgIHNldFJlZmVycmVyKG1hcElEOiBzdHJpbmcsIHJlZmVycmVyOiBzdHJpbmcpIHtcbiAgICAgICAgdGhpcy5yZWZlcnJlciA9IHJlZmVycmVyO1xuICAgIH1cblxuICAgIHNwcml0ZUxvYWRlZChtYXBJZDogc3RyaW5nLCBib29sOiBib29sZWFuKSB7XG4gICAgICAgIHRoaXMuaXNTcHJpdGVMb2FkZWRbbWFwSWRdID0gYm9vbDtcbiAgICAgICAgZm9yIChjb25zdCB3b3JrZXJTb3VyY2UgaW4gdGhpcy53b3JrZXJTb3VyY2VzW21hcElkXSkge1xuICAgICAgICAgICAgY29uc3Qgd3MgPSB0aGlzLndvcmtlclNvdXJjZXNbbWFwSWRdW3dvcmtlclNvdXJjZV07XG4gICAgICAgICAgICBmb3IgKGNvbnN0IHNvdXJjZSBpbiB3cykge1xuICAgICAgICAgICAgICAgIGlmICh3c1tzb3VyY2VdIGluc3RhbmNlb2YgVmVjdG9yVGlsZVdvcmtlclNvdXJjZSkge1xuICAgICAgICAgICAgICAgICAgICB3c1tzb3VyY2VdLmlzU3ByaXRlTG9hZGVkID0gYm9vbDtcbiAgICAgICAgICAgICAgICAgICAgd3Nbc291cmNlXS5maXJlKG5ldyBFdmVudCgnaXNTcHJpdGVMb2FkZWQnKSk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgfVxuXG4gICAgc2V0SW1hZ2VzKG1hcElkOiBzdHJpbmcsIGltYWdlczogQXJyYXk8c3RyaW5nPiwgY2FsbGJhY2s6IFdvcmtlclRpbGVDYWxsYmFjaykge1xuICAgICAgICB0aGlzLmF2YWlsYWJsZUltYWdlc1ttYXBJZF0gPSBpbWFnZXM7XG4gICAgICAgIGZvciAoY29uc3Qgd29ya2VyU291cmNlIGluIHRoaXMud29ya2VyU291cmNlc1ttYXBJZF0pIHtcbiAgICAgICAgICAgIGNvbnN0IHdzID0gdGhpcy53b3JrZXJTb3VyY2VzW21hcElkXVt3b3JrZXJTb3VyY2VdO1xuICAgICAgICAgICAgZm9yIChjb25zdCBzb3VyY2UgaW4gd3MpIHtcbiAgICAgICAgICAgICAgICB3c1tzb3VyY2VdLmF2YWlsYWJsZUltYWdlcyA9IGltYWdlcztcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBjYWxsYmFjaygpO1xuICAgIH1cblxuICAgIGVuYWJsZVRlcnJhaW4obWFwSWQ6IHN0cmluZywgZW5hYmxlOiBib29sZWFuLCBjYWxsYmFjazogV29ya2VyVGlsZUNhbGxiYWNrKSB7XG4gICAgICAgIHRoaXMudGVycmFpbiA9IGVuYWJsZTtcbiAgICAgICAgY2FsbGJhY2soKTtcbiAgICB9XG5cbiAgICBzZXRQcm9qZWN0aW9uKG1hcElkOiBzdHJpbmcsIGNvbmZpZzogUHJvamVjdGlvblNwZWNpZmljYXRpb24pIHtcbiAgICAgICAgdGhpcy5wcm9qZWN0aW9uc1ttYXBJZF0gPSBnZXRQcm9qZWN0aW9uKGNvbmZpZyk7XG4gICAgfVxuXG4gICAgc2V0TGF5ZXJzKG1hcElkOiBzdHJpbmcsIGxheWVyczogQXJyYXk8TGF5ZXJTcGVjaWZpY2F0aW9uPiwgY2FsbGJhY2s6IFdvcmtlclRpbGVDYWxsYmFjaykge1xuICAgICAgICB0aGlzLmdldExheWVySW5kZXgobWFwSWQpLnJlcGxhY2UobGF5ZXJzKTtcbiAgICAgICAgY2FsbGJhY2soKTtcbiAgICB9XG5cbiAgICB1cGRhdGVMYXllcnMobWFwSWQ6IHN0cmluZywgcGFyYW1zOiB7bGF5ZXJzOiBBcnJheTxMYXllclNwZWNpZmljYXRpb24+LCByZW1vdmVkSWRzOiBBcnJheTxzdHJpbmc+fSwgY2FsbGJhY2s6IFdvcmtlclRpbGVDYWxsYmFjaykge1xuICAgICAgICB0aGlzLmdldExheWVySW5kZXgobWFwSWQpLnVwZGF0ZShwYXJhbXMubGF5ZXJzLCBwYXJhbXMucmVtb3ZlZElkcyk7XG4gICAgICAgIGNhbGxiYWNrKCk7XG4gICAgfVxuXG4gICAgbG9hZFRpbGUobWFwSWQ6IHN0cmluZywgcGFyYW1zOiBXb3JrZXJUaWxlUGFyYW1ldGVycyAmIHt0eXBlOiBzdHJpbmd9LCBjYWxsYmFjazogV29ya2VyVGlsZUNhbGxiYWNrKSB7XG4gICAgICAgIGFzc2VydChwYXJhbXMudHlwZSk7XG4gICAgICAgIGNvbnN0IHAgPSB0aGlzLmVuYWJsZVRlcnJhaW4gPyBleHRlbmQoe2VuYWJsZVRlcnJhaW46IHRoaXMudGVycmFpbn0sIHBhcmFtcykgOiBwYXJhbXM7XG4gICAgICAgIHAucHJvamVjdGlvbiA9IHRoaXMucHJvamVjdGlvbnNbbWFwSWRdIHx8IHRoaXMuZGVmYXVsdFByb2plY3Rpb247XG4gICAgICAgIHRoaXMuZ2V0V29ya2VyU291cmNlKG1hcElkLCBwYXJhbXMudHlwZSwgcGFyYW1zLnNvdXJjZSkubG9hZFRpbGUocCwgY2FsbGJhY2spO1xuICAgIH1cblxuICAgIGxvYWRERU1UaWxlKG1hcElkOiBzdHJpbmcsIHBhcmFtczogV29ya2VyREVNVGlsZVBhcmFtZXRlcnMsIGNhbGxiYWNrOiBXb3JrZXJERU1UaWxlQ2FsbGJhY2spIHtcbiAgICAgICAgY29uc3QgcCA9IHRoaXMuZW5hYmxlVGVycmFpbiA/IGV4dGVuZCh7YnVpbGRRdWFkVHJlZTogdGhpcy50ZXJyYWlufSwgcGFyYW1zKSA6IHBhcmFtcztcbiAgICAgICAgdGhpcy5nZXRERU1Xb3JrZXJTb3VyY2UobWFwSWQsIHBhcmFtcy5zb3VyY2UpLmxvYWRUaWxlKHAsIGNhbGxiYWNrKTtcbiAgICB9XG5cbiAgICByZWxvYWRUaWxlKG1hcElkOiBzdHJpbmcsIHBhcmFtczogV29ya2VyVGlsZVBhcmFtZXRlcnMgJiB7dHlwZTogc3RyaW5nfSwgY2FsbGJhY2s6IFdvcmtlclRpbGVDYWxsYmFjaykge1xuICAgICAgICBhc3NlcnQocGFyYW1zLnR5cGUpO1xuICAgICAgICBjb25zdCBwID0gdGhpcy5lbmFibGVUZXJyYWluID8gZXh0ZW5kKHtlbmFibGVUZXJyYWluOiB0aGlzLnRlcnJhaW59LCBwYXJhbXMpIDogcGFyYW1zO1xuICAgICAgICBwLnByb2plY3Rpb24gPSB0aGlzLnByb2plY3Rpb25zW21hcElkXSB8fCB0aGlzLmRlZmF1bHRQcm9qZWN0aW9uO1xuICAgICAgICB0aGlzLmdldFdvcmtlclNvdXJjZShtYXBJZCwgcGFyYW1zLnR5cGUsIHBhcmFtcy5zb3VyY2UpLnJlbG9hZFRpbGUocCwgY2FsbGJhY2spO1xuICAgIH1cblxuICAgIGFib3J0VGlsZShtYXBJZDogc3RyaW5nLCBwYXJhbXM6IFRpbGVQYXJhbWV0ZXJzICYge3R5cGU6IHN0cmluZ30sIGNhbGxiYWNrOiBXb3JrZXJUaWxlQ2FsbGJhY2spIHtcbiAgICAgICAgYXNzZXJ0KHBhcmFtcy50eXBlKTtcbiAgICAgICAgdGhpcy5nZXRXb3JrZXJTb3VyY2UobWFwSWQsIHBhcmFtcy50eXBlLCBwYXJhbXMuc291cmNlKS5hYm9ydFRpbGUocGFyYW1zLCBjYWxsYmFjayk7XG4gICAgfVxuXG4gICAgcmVtb3ZlVGlsZShtYXBJZDogc3RyaW5nLCBwYXJhbXM6IFRpbGVQYXJhbWV0ZXJzICYge3R5cGU6IHN0cmluZ30sIGNhbGxiYWNrOiBXb3JrZXJUaWxlQ2FsbGJhY2spIHtcbiAgICAgICAgYXNzZXJ0KHBhcmFtcy50eXBlKTtcbiAgICAgICAgdGhpcy5nZXRXb3JrZXJTb3VyY2UobWFwSWQsIHBhcmFtcy50eXBlLCBwYXJhbXMuc291cmNlKS5yZW1vdmVUaWxlKHBhcmFtcywgY2FsbGJhY2spO1xuICAgIH1cblxuICAgIHJlbW92ZVNvdXJjZShtYXBJZDogc3RyaW5nLCBwYXJhbXM6IHtzb3VyY2U6IHN0cmluZ30gJiB7dHlwZTogc3RyaW5nfSwgY2FsbGJhY2s6IFdvcmtlclRpbGVDYWxsYmFjaykge1xuICAgICAgICBhc3NlcnQocGFyYW1zLnR5cGUpO1xuICAgICAgICBhc3NlcnQocGFyYW1zLnNvdXJjZSk7XG5cbiAgICAgICAgaWYgKCF0aGlzLndvcmtlclNvdXJjZXNbbWFwSWRdIHx8XG4gICAgICAgICAgICAhdGhpcy53b3JrZXJTb3VyY2VzW21hcElkXVtwYXJhbXMudHlwZV0gfHxcbiAgICAgICAgICAgICF0aGlzLndvcmtlclNvdXJjZXNbbWFwSWRdW3BhcmFtcy50eXBlXVtwYXJhbXMuc291cmNlXSkge1xuICAgICAgICAgICAgcmV0dXJuO1xuICAgICAgICB9XG5cbiAgICAgICAgY29uc3Qgd29ya2VyID0gdGhpcy53b3JrZXJTb3VyY2VzW21hcElkXVtwYXJhbXMudHlwZV1bcGFyYW1zLnNvdXJjZV07XG4gICAgICAgIGRlbGV0ZSB0aGlzLndvcmtlclNvdXJjZXNbbWFwSWRdW3BhcmFtcy50eXBlXVtwYXJhbXMuc291cmNlXTtcblxuICAgICAgICBpZiAod29ya2VyLnJlbW92ZVNvdXJjZSAhPT0gdW5kZWZpbmVkKSB7XG4gICAgICAgICAgICB3b3JrZXIucmVtb3ZlU291cmNlKHBhcmFtcywgY2FsbGJhY2spO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgY2FsbGJhY2soKTtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIC8qKlxuICAgICAqIExvYWQgYSB7QGxpbmsgV29ya2VyU291cmNlfSBzY3JpcHQgYXQgcGFyYW1zLnVybC4gIFRoZSBzY3JpcHQgaXMgcnVuXG4gICAgICogKHVzaW5nIGltcG9ydFNjcmlwdHMpIHdpdGggYHJlZ2lzdGVyV29ya2VyU291cmNlYCBpbiBzY29wZSwgd2hpY2ggaXMgYVxuICAgICAqIGZ1bmN0aW9uIHRha2luZyBgKG5hbWUsIHdvcmtlclNvdXJjZU9iamVjdClgLlxuICAgICAqICBAcHJpdmF0ZVxuICAgICAqL1xuICAgIGxvYWRXb3JrZXJTb3VyY2UobWFwOiBzdHJpbmcsIHBhcmFtczogeyB1cmw6IHN0cmluZyB9LCBjYWxsYmFjazogQ2FsbGJhY2s8dm9pZD4pIHtcbiAgICAgICAgdHJ5IHtcbiAgICAgICAgICAgIHRoaXMuc2VsZi5pbXBvcnRTY3JpcHRzKHBhcmFtcy51cmwpO1xuICAgICAgICAgICAgY2FsbGJhY2soKTtcbiAgICAgICAgfSBjYXRjaCAoZSkge1xuICAgICAgICAgICAgY2FsbGJhY2soZS50b1N0cmluZygpKTtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIHN5bmNSVExQbHVnaW5TdGF0ZShtYXA6IHN0cmluZywgc3RhdGU6IFBsdWdpblN0YXRlLCBjYWxsYmFjazogQ2FsbGJhY2s8Ym9vbGVhbj4pIHtcbiAgICAgICAgdHJ5IHtcbiAgICAgICAgICAgIGdsb2JhbFJUTFRleHRQbHVnaW4uc2V0U3RhdGUoc3RhdGUpO1xuICAgICAgICAgICAgY29uc3QgcGx1Z2luVVJMID0gZ2xvYmFsUlRMVGV4dFBsdWdpbi5nZXRQbHVnaW5VUkwoKTtcbiAgICAgICAgICAgIGlmIChcbiAgICAgICAgICAgICAgICBnbG9iYWxSVExUZXh0UGx1Z2luLmlzTG9hZGVkKCkgJiZcbiAgICAgICAgICAgICAgICAhZ2xvYmFsUlRMVGV4dFBsdWdpbi5pc1BhcnNlZCgpICYmXG4gICAgICAgICAgICAgICAgcGx1Z2luVVJMICE9IG51bGwgLy8gTm90IHBvc3NpYmxlIHdoZW4gYGlzTG9hZGVkYCBpcyB0cnVlLCBidXQga2VlcHMgZmxvdyBoYXBweVxuICAgICAgICAgICAgKSB7XG4gICAgICAgICAgICAgICAgdGhpcy5zZWxmLmltcG9ydFNjcmlwdHMocGx1Z2luVVJMKTtcbiAgICAgICAgICAgICAgICBjb25zdCBjb21wbGV0ZSA9IGdsb2JhbFJUTFRleHRQbHVnaW4uaXNQYXJzZWQoKTtcbiAgICAgICAgICAgICAgICBjb25zdCBlcnJvciA9IGNvbXBsZXRlID8gdW5kZWZpbmVkIDogbmV3IEVycm9yKGBSVEwgVGV4dCBQbHVnaW4gZmFpbGVkIHRvIGltcG9ydCBzY3JpcHRzIGZyb20gJHtwbHVnaW5VUkx9YCk7XG4gICAgICAgICAgICAgICAgY2FsbGJhY2soZXJyb3IsIGNvbXBsZXRlKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSBjYXRjaCAoZSkge1xuICAgICAgICAgICAgY2FsbGJhY2soZS50b1N0cmluZygpKTtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIGdldEF2YWlsYWJsZUltYWdlcyhtYXBJZDogc3RyaW5nKTogQXJyYXk8c3RyaW5nPiB7XG4gICAgICAgIGxldCBhdmFpbGFibGVJbWFnZXMgPSB0aGlzLmF2YWlsYWJsZUltYWdlc1ttYXBJZF07XG5cbiAgICAgICAgaWYgKCFhdmFpbGFibGVJbWFnZXMpIHtcbiAgICAgICAgICAgIGF2YWlsYWJsZUltYWdlcyA9IFtdO1xuICAgICAgICB9XG5cbiAgICAgICAgcmV0dXJuIGF2YWlsYWJsZUltYWdlcztcbiAgICB9XG5cbiAgICBnZXRMYXllckluZGV4KG1hcElkOiBzdHJpbmcpOiBTdHlsZUxheWVySW5kZXgge1xuICAgICAgICBsZXQgbGF5ZXJJbmRleGVzID0gdGhpcy5sYXllckluZGV4ZXNbbWFwSWRdO1xuICAgICAgICBpZiAoIWxheWVySW5kZXhlcykge1xuICAgICAgICAgICAgbGF5ZXJJbmRleGVzID0gdGhpcy5sYXllckluZGV4ZXNbbWFwSWRdID0gbmV3IFN0eWxlTGF5ZXJJbmRleCgpO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiBsYXllckluZGV4ZXM7XG4gICAgfVxuXG4gICAgZ2V0V29ya2VyU291cmNlKG1hcElkOiBzdHJpbmcsIHR5cGU6IHN0cmluZywgc291cmNlOiBzdHJpbmcpOiBXb3JrZXJTb3VyY2Uge1xuICAgICAgICBpZiAoIXRoaXMud29ya2VyU291cmNlc1ttYXBJZF0pXG4gICAgICAgICAgICB0aGlzLndvcmtlclNvdXJjZXNbbWFwSWRdID0ge307XG4gICAgICAgIGlmICghdGhpcy53b3JrZXJTb3VyY2VzW21hcElkXVt0eXBlXSlcbiAgICAgICAgICAgIHRoaXMud29ya2VyU291cmNlc1ttYXBJZF1bdHlwZV0gPSB7fTtcblxuICAgICAgICBpZiAoIXRoaXMud29ya2VyU291cmNlc1ttYXBJZF1bdHlwZV1bc291cmNlXSkge1xuICAgICAgICAgICAgLy8gdXNlIGEgd3JhcHBlZCBhY3RvciBzbyB0aGF0IHdlIGNhbiBhdHRhY2ggYSB0YXJnZXQgbWFwSWQgcGFyYW1cbiAgICAgICAgICAgIC8vIHRvIGFueSBtZXNzYWdlcyBpbnZva2VkIGJ5IHRoZSBXb3JrZXJTb3VyY2VcbiAgICAgICAgICAgIGNvbnN0IGFjdG9yID0ge1xuICAgICAgICAgICAgICAgIHNlbmQ6ICh0eXBlLCBkYXRhLCBjYWxsYmFjaywgXywgbXVzdFF1ZXVlLCBtZXRhZGF0YSkgPT4ge1xuICAgICAgICAgICAgICAgICAgICB0aGlzLmFjdG9yLnNlbmQodHlwZSwgZGF0YSwgY2FsbGJhY2ssIG1hcElkLCBtdXN0UXVldWUsIG1ldGFkYXRhKTtcbiAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgIHNjaGVkdWxlcjogdGhpcy5hY3Rvci5zY2hlZHVsZXJcbiAgICAgICAgICAgIH07XG4gICAgICAgICAgICB0aGlzLndvcmtlclNvdXJjZXNbbWFwSWRdW3R5cGVdW3NvdXJjZV0gPSBuZXcgKHRoaXMud29ya2VyU291cmNlVHlwZXNbdHlwZV06IGFueSkoKGFjdG9yOiBhbnkpLCB0aGlzLmdldExheWVySW5kZXgobWFwSWQpLCB0aGlzLmdldEF2YWlsYWJsZUltYWdlcyhtYXBJZCksIHRoaXMuaXNTcHJpdGVMb2FkZWRbbWFwSWRdKTtcbiAgICAgICAgfVxuXG4gICAgICAgIHJldHVybiB0aGlzLndvcmtlclNvdXJjZXNbbWFwSWRdW3R5cGVdW3NvdXJjZV07XG4gICAgfVxuXG4gICAgZ2V0REVNV29ya2VyU291cmNlKG1hcElkOiBzdHJpbmcsIHNvdXJjZTogc3RyaW5nKTogUmFzdGVyREVNVGlsZVdvcmtlclNvdXJjZSB7XG4gICAgICAgIGlmICghdGhpcy5kZW1Xb3JrZXJTb3VyY2VzW21hcElkXSlcbiAgICAgICAgICAgIHRoaXMuZGVtV29ya2VyU291cmNlc1ttYXBJZF0gPSB7fTtcblxuICAgICAgICBpZiAoIXRoaXMuZGVtV29ya2VyU291cmNlc1ttYXBJZF1bc291cmNlXSkge1xuICAgICAgICAgICAgdGhpcy5kZW1Xb3JrZXJTb3VyY2VzW21hcElkXVtzb3VyY2VdID0gbmV3IFJhc3RlckRFTVRpbGVXb3JrZXJTb3VyY2UoKTtcbiAgICAgICAgfVxuXG4gICAgICAgIHJldHVybiB0aGlzLmRlbVdvcmtlclNvdXJjZXNbbWFwSWRdW3NvdXJjZV07XG4gICAgfVxuXG4gICAgZW5mb3JjZUNhY2hlU2l6ZUxpbWl0KG1hcElkOiBzdHJpbmcsIGxpbWl0OiBudW1iZXIpIHtcbiAgICAgICAgZW5mb3JjZUNhY2hlU2l6ZUxpbWl0KGxpbWl0KTtcbiAgICB9XG5cbiAgICBnZXRXb3JrZXJQZXJmb3JtYW5jZU1ldHJpY3MobWFwSWQ6IHN0cmluZywgcGFyYW1zOiBhbnksIGNhbGxiYWNrOiAoZXJyb3I6ID9FcnJvciwgcmVzdWx0OiA/T2JqZWN0KSA9PiB2b2lkKSB7XG4gICAgICAgIGNhbGxiYWNrKHVuZGVmaW5lZCwgUGVyZm9ybWFuY2VVdGlscy5nZXRXb3JrZXJQZXJmb3JtYW5jZU1ldHJpY3MoKSk7XG4gICAgfVxufVxuXG4vKiBnbG9iYWwgc2VsZiwgV29ya2VyR2xvYmFsU2NvcGUgKi9cbmlmICh0eXBlb2YgV29ya2VyR2xvYmFsU2NvcGUgIT09ICd1bmRlZmluZWQnICYmXG4gICAgdHlwZW9mIHNlbGYgIT09ICd1bmRlZmluZWQnICYmXG4gICAgc2VsZiBpbnN0YW5jZW9mIFdvcmtlckdsb2JhbFNjb3BlKSB7XG4gICAgc2VsZi53b3JrZXIgPSBuZXcgV29ya2VyKHNlbGYpO1xufVxuIl0sIm5hbWVzIjpbInN0cmluZ2lmeSIsIm9iaiIsInVuZGVmaW5lZCIsIkpTT04iLCJBcnJheSIsImlzQXJyYXkiLCJzdHIiLCJ2YWwiLCJrZXkiLCJPYmplY3QiLCJrZXlzIiwic29ydCIsImdldEtleSIsImxheWVyIiwiayIsInJlZlByb3BlcnRpZXMiLCJncm91cEJ5TGF5b3V0IiwibGF5ZXJzIiwiY2FjaGVkS2V5cyIsImdyb3VwcyIsImkiLCJsZW5ndGgiLCJpZCIsImdyb3VwIiwicHVzaCIsInJlc3VsdCIsIlN0eWxlTGF5ZXJJbmRleCIsImNvbnN0cnVjdG9yIiwibGF5ZXJDb25maWdzIiwia2V5Q2FjaGUiLCJyZXBsYWNlIiwiX2xheWVyQ29uZmlncyIsIl9sYXllcnMiLCJ1cGRhdGUiLCJyZW1vdmVkSWRzIiwibGF5ZXJDb25maWciLCJjcmVhdGVTdHlsZUxheWVyIiwiY29tcGlsZUZpbHRlciIsImZhbWlsaWVzQnlTb3VyY2UiLCJ2YWx1ZXMiLCJtYXAiLCJ2aXNpYmlsaXR5Iiwic291cmNlSWQiLCJzb3VyY2UiLCJzb3VyY2VHcm91cCIsInNvdXJjZUxheWVySWQiLCJzb3VyY2VMYXllciIsInNvdXJjZUxheWVyRmFtaWxpZXMiLCJSYXN0ZXJERU1UaWxlV29ya2VyU291cmNlIiwibG9hZFRpbGUiLCJwYXJhbXMiLCJjYWxsYmFjayIsInVpZCIsImVuY29kaW5nIiwicmF3SW1hZ2VEYXRhIiwicGFkZGluZyIsImJ1aWxkUXVhZFRyZWUiLCJpbWFnZVBpeGVscyIsIndpbmRvdyIsIkltYWdlQml0bWFwIiwiZ2V0SW1hZ2VEYXRhIiwiZGVtIiwiREVNRGF0YSIsImltZ0JpdG1hcCIsIm9mZnNjcmVlbkNhbnZhcyIsIm9mZnNjcmVlbkNhbnZhc0NvbnRleHQiLCJPZmZzY3JlZW5DYW52YXMiLCJ3aWR0aCIsImhlaWdodCIsImdldENvbnRleHQiLCJkcmF3SW1hZ2UiLCJpbWdEYXRhIiwiY2xlYXJSZWN0IiwibW9kdWxlIiwicmV3aW5kIiwiZ2oiLCJvdXRlciIsInR5cGUiLCJmZWF0dXJlcyIsImdlb21ldHJpZXMiLCJnZW9tZXRyeSIsInJld2luZFJpbmdzIiwiY29vcmRpbmF0ZXMiLCJyaW5ncyIsInJld2luZFJpbmciLCJyaW5nIiwiZGlyIiwiYXJlYSIsImVyciIsImxlbiIsImoiLCJtIiwiTWF0aCIsImFicyIsInJldmVyc2UiLCJ0b0dlb0pTT04iLCJtdnQiLCJWZWN0b3JUaWxlRmVhdHVyZSIsInByb3RvdHlwZSIsIkZlYXR1cmVXcmFwcGVyIiwiZmVhdHVyZSIsIl9mZWF0dXJlIiwiZXh0ZW50IiwiRVhURU5UIiwicHJvcGVydGllcyIsInRhZ3MiLCJpc05hTiIsInBhcnNlSW50IiwibG9hZEdlb21ldHJ5IiwicG9pbnQiLCJQb2ludCIsIm5ld1JpbmciLCJ4IiwieSIsInoiLCJjYWxsIiwiR2VvSlNPTldyYXBwZXIiLCJuYW1lIiwiX2ZlYXR1cmVzIiwicmVxdWlyZSIsIm9wdGlvbnMiLCJyYXdHZW9tZXRyeSIsImJib3giLCJ4MSIsIkluZmluaXR5IiwieDIiLCJ5MSIsInkyIiwiY29vcmQiLCJtaW4iLCJtYXgiLCJQYmYiLCJmcm9tVmVjdG9yVGlsZUpzIiwidnRQYmYiLCJmcm9tR2VvanNvblZ0IiwidGlsZSIsIm91dCIsIndyaXRlVGlsZSIsImZpbmlzaCIsImwiLCJ2ZXJzaW9uIiwicGJmIiwid3JpdGVNZXNzYWdlIiwid3JpdGVMYXllciIsIndyaXRlVmFyaW50RmllbGQiLCJ3cml0ZVN0cmluZ0ZpZWxkIiwiY29udGV4dCIsImtleWNhY2hlIiwidmFsdWVjYWNoZSIsIndyaXRlRmVhdHVyZSIsIndyaXRlVmFsdWUiLCJ3cml0ZVByb3BlcnRpZXMiLCJ3cml0ZUdlb21ldHJ5IiwidmFsdWUiLCJrZXlJbmRleCIsIndyaXRlVmFyaW50IiwidmFsdWVLZXkiLCJ2YWx1ZUluZGV4IiwiY29tbWFuZCIsImNtZCIsInppZ3phZyIsIm51bSIsInIiLCJjb3VudCIsImxpbmVDb3VudCIsImR4IiwiZHkiLCJ3cml0ZUJvb2xlYW5GaWVsZCIsIndyaXRlRG91YmxlRmllbGQiLCJ3cml0ZVNWYXJpbnRGaWVsZCIsInNvcnRLRCIsImlkcyIsImNvb3JkcyIsIm5vZGVTaXplIiwibGVmdCIsInJpZ2h0IiwiZGVwdGgiLCJzZWxlY3QiLCJpbmMiLCJuIiwibG9nIiwicyIsImV4cCIsInNkIiwic3FydCIsIm5ld0xlZnQiLCJmbG9vciIsIm5ld1JpZ2h0IiwidCIsInN3YXBJdGVtIiwic3dhcCIsImFyciIsInRtcCIsInJhbmdlIiwibWluWCIsIm1pblkiLCJtYXhYIiwibWF4WSIsInN0YWNrIiwiYXhpcyIsInBvcCIsIm5leHRBeGlzIiwid2l0aGluIiwicXgiLCJxeSIsInIyIiwic3FEaXN0IiwiYXgiLCJheSIsImJ4IiwiYnkiLCJkZWZhdWx0R2V0WCIsInAiLCJkZWZhdWx0R2V0WSIsIktEQnVzaCIsInBvaW50cyIsImdldFgiLCJnZXRZIiwiQXJyYXlUeXBlIiwiRmxvYXQ2NEFycmF5IiwiSW5kZXhBcnJheVR5cGUiLCJVaW50MTZBcnJheSIsIlVpbnQzMkFycmF5IiwiZGVmYXVsdE9wdGlvbnMiLCJtaW5ab29tIiwibWF4Wm9vbSIsIm1pblBvaW50cyIsInJhZGl1cyIsImdlbmVyYXRlSWQiLCJyZWR1Y2UiLCJwcm9wcyIsImZyb3VuZCIsIkZsb2F0MzJBcnJheSIsIlN1cGVyY2x1c3RlciIsImV4dGVuZCIsImNyZWF0ZSIsInRyZWVzIiwibG9hZCIsImNvbnNvbGUiLCJ0aW1lIiwidGltZXJJZCIsImNsdXN0ZXJzIiwiY3JlYXRlUG9pbnRDbHVzdGVyIiwidGltZUVuZCIsIm5vdyIsIkRhdGUiLCJfY2x1c3RlciIsImdldENsdXN0ZXJzIiwiem9vbSIsIm1pbkxuZyIsIm1pbkxhdCIsIm1heExuZyIsIm1heExhdCIsImVhc3Rlcm5IZW0iLCJ3ZXN0ZXJuSGVtIiwiY29uY2F0IiwidHJlZSIsIl9saW1pdFpvb20iLCJsbmdYIiwibGF0WSIsImMiLCJudW1Qb2ludHMiLCJnZXRDbHVzdGVySlNPTiIsImluZGV4IiwiZ2V0Q2hpbGRyZW4iLCJjbHVzdGVySWQiLCJvcmlnaW5JZCIsIl9nZXRPcmlnaW5JZCIsIm9yaWdpblpvb20iLCJfZ2V0T3JpZ2luWm9vbSIsImVycm9yTXNnIiwiRXJyb3IiLCJvcmlnaW4iLCJwb3ciLCJjaGlsZHJlbiIsInBhcmVudElkIiwiZ2V0TGVhdmVzIiwibGltaXQiLCJvZmZzZXQiLCJsZWF2ZXMiLCJfYXBwZW5kTGVhdmVzIiwiZ2V0VGlsZSIsInoyIiwidG9wIiwiYm90dG9tIiwiX2FkZFRpbGVGZWF0dXJlcyIsImdldENsdXN0ZXJFeHBhbnNpb25ab29tIiwiZXhwYW5zaW9uWm9vbSIsImNsdXN0ZXJfaWQiLCJza2lwcGVkIiwiY2hpbGQiLCJjbHVzdGVyIiwicG9pbnRfY291bnQiLCJpc0NsdXN0ZXIiLCJweCIsInB5IiwiZ2V0Q2x1c3RlclByb3BlcnRpZXMiLCJmIiwicm91bmQiLCJuZWlnaGJvcklkcyIsIm51bVBvaW50c09yaWdpbiIsIm5laWdoYm9ySWQiLCJiIiwid3giLCJ3eSIsImNsdXN0ZXJQcm9wZXJ0aWVzIiwiX21hcCIsIm51bVBvaW50czIiLCJjcmVhdGVDbHVzdGVyIiwiY2xvbmUiLCJvcmlnaW5hbCIsInhMbmciLCJ5TGF0IiwiYWJicmV2IiwicG9pbnRfY291bnRfYWJicmV2aWF0ZWQiLCJsbmciLCJsYXQiLCJzaW4iLCJQSSIsImF0YW4iLCJkZXN0Iiwic3JjIiwic2ltcGxpZnkiLCJmaXJzdCIsImxhc3QiLCJzcVRvbGVyYW5jZSIsIm1heFNxRGlzdCIsIm1pZCIsIm1pblBvc1RvTWlkIiwiZCIsImdldFNxU2VnRGlzdCIsInBvc1RvTWlkIiwiY3JlYXRlRmVhdHVyZSIsImdlb20iLCJjYWxjQkJveCIsImNhbGNMaW5lQkJveCIsImNvbnZlcnQiLCJkYXRhIiwiY29udmVydEZlYXR1cmUiLCJnZW9qc29uIiwidG9sZXJhbmNlIiwicHJvbW90ZUlkIiwiY29udmVydFBvaW50IiwiY29udmVydExpbmUiLCJsaW5lTWV0cmljcyIsImNvbnZlcnRMaW5lcyIsInBvbHlnb24iLCJwcm9qZWN0WCIsInByb2plY3RZIiwiaXNQb2x5Z29uIiwieDAiLCJ5MCIsInNpemUiLCJzdGFydCIsImVuZCIsImNsaXAiLCJzY2FsZSIsImsxIiwiazIiLCJtaW5BbGwiLCJtYXhBbGwiLCJjbGlwcGVkIiwibmV3R2VvbWV0cnkiLCJjbGlwUG9pbnRzIiwiY2xpcExpbmUiLCJjbGlwTGluZXMiLCJuZXdHZW9tIiwiYSIsInRyYWNrTWV0cmljcyIsInNsaWNlIiwibmV3U2xpY2UiLCJpbnRlcnNlY3QiLCJpbnRlcnNlY3RYIiwiaW50ZXJzZWN0WSIsInNlZ0xlbiIsImF6IiwiZXhpdGVkIiwiYWRkUG9pbnQiLCJsaW5lIiwid3JhcCIsImJ1ZmZlciIsIm1lcmdlZCIsInNoaWZ0RmVhdHVyZUNvb3JkcyIsIm5ld0ZlYXR1cmVzIiwic2hpZnRDb29yZHMiLCJuZXdQb2x5Z29uIiwibmV3UG9pbnRzIiwidHJhbnNmb3JtVGlsZSIsInRyYW5zZm9ybWVkIiwidHgiLCJ0eSIsInRyYW5zZm9ybVBvaW50IiwiY3JlYXRlVGlsZSIsIm51bVNpbXBsaWZpZWQiLCJudW1GZWF0dXJlcyIsImFkZEZlYXR1cmUiLCJzaW1wbGlmaWVkIiwiYWRkTGluZSIsInRpbGVGZWF0dXJlIiwiaXNPdXRlciIsImNsb2Nrd2lzZSIsImdlb2pzb252dCIsIkdlb0pTT05WVCIsImRlYnVnIiwidGlsZXMiLCJ0aWxlQ29vcmRzIiwiaW5kZXhNYXhab29tIiwiaW5kZXhNYXhQb2ludHMiLCJzdGF0cyIsInRvdGFsIiwic3BsaXRUaWxlIiwiY3oiLCJjeCIsImN5IiwidG9JRCIsImszIiwiazQiLCJ0bCIsImJsIiwidHIiLCJiciIsInRyYW5zZm9ybSIsInowIiwicGFyZW50IiwibG9hZEdlb0pTT05UaWxlIiwiY2Fub25pY2FsIiwidGlsZUlEIiwiX2dlb0pTT05JbmRleCIsImdlb0pTT05UaWxlIiwiZ2VvanNvbldyYXBwZXIiLCJ2dHBiZiIsImJ5dGVPZmZzZXQiLCJieXRlTGVuZ3RoIiwiVWludDhBcnJheSIsInZlY3RvclRpbGUiLCJyYXdEYXRhIiwiR2VvSlNPTldvcmtlclNvdXJjZSIsIlZlY3RvclRpbGVXb3JrZXJTb3VyY2UiLCJhY3RvciIsImxheWVySW5kZXgiLCJhdmFpbGFibGVJbWFnZXMiLCJpc1Nwcml0ZUxvYWRlZCIsImxvYWRHZW9KU09OIiwibG9hZERhdGEiLCJyZXF1ZXN0UGFyYW0iLCJyZXF1ZXN0IiwicGVyZiIsImNvbGxlY3RSZXNvdXJjZVRpbWluZyIsImZpbHRlciIsImNvbXBpbGVkIiwiY3JlYXRlRXhwcmVzc2lvbiIsIm92ZXJyaWRhYmxlIiwidHJhbnNpdGlvbiIsIm1lc3NhZ2UiLCJqb2luIiwiZXZhbHVhdGUiLCJnZXRTdXBlcmNsdXN0ZXJPcHRpb25zIiwiZ2VvanNvblZ0T3B0aW9ucyIsImxvYWRlZCIsInJlc291cmNlVGltaW5nRGF0YSIsImdldFBlcmZvcm1hbmNlTWVhc3VyZW1lbnQiLCJyZXNvdXJjZVRpbWluZyIsInBhcnNlIiwicmVsb2FkVGlsZSIsImdldEpTT04iLCJlIiwiZ2V0Q2x1c3RlckNoaWxkcmVuIiwiZ2V0Q2x1c3RlckxlYXZlcyIsInN1cGVyY2x1c3Rlck9wdGlvbnMiLCJtYXBFeHByZXNzaW9ucyIsInJlZHVjZUV4cHJlc3Npb25zIiwiZ2xvYmFscyIsImFjY3VtdWxhdGVkIiwicHJvcGVydHlOYW1lcyIsIm9wZXJhdG9yIiwibWFwRXhwcmVzc2lvbiIsIm1hcEV4cHJlc3Npb25QYXJzZWQiLCJyZWR1Y2VFeHByZXNzaW9uUGFyc2VkIiwicG9pbnRQcm9wZXJ0aWVzIiwiV29ya2VyIiwic2VsZiIsIkFjdG9yIiwibGF5ZXJJbmRleGVzIiwicHJvamVjdGlvbnMiLCJkZWZhdWx0UHJvamVjdGlvbiIsImdldFByb2plY3Rpb24iLCJ3b3JrZXJTb3VyY2VUeXBlcyIsInZlY3RvciIsIndvcmtlclNvdXJjZXMiLCJkZW1Xb3JrZXJTb3VyY2VzIiwicmVnaXN0ZXJXb3JrZXJTb3VyY2UiLCJXb3JrZXJTb3VyY2UiLCJyZWdpc3RlclJUTFRleHRQbHVnaW4iLCJydGxUZXh0UGx1Z2luIiwiZ2xvYmFsUlRMVGV4dFBsdWdpbiIsImlzUGFyc2VkIiwiYXBwbHlBcmFiaWNTaGFwaW5nIiwicHJvY2Vzc0JpZGlyZWN0aW9uYWxUZXh0IiwicHJvY2Vzc1N0eWxlZEJpZGlyZWN0aW9uYWxUZXh0IiwiY2xlYXJDYWNoZXMiLCJtYXBJZCIsInVudXNlZCIsImNoZWNrSWZSZWFkeSIsIm1hcElEIiwic2V0UmVmZXJyZXIiLCJyZWZlcnJlciIsInNwcml0ZUxvYWRlZCIsImJvb2wiLCJ3b3JrZXJTb3VyY2UiLCJ3cyIsImZpcmUiLCJFdmVudCIsInNldEltYWdlcyIsImltYWdlcyIsImVuYWJsZVRlcnJhaW4iLCJlbmFibGUiLCJ0ZXJyYWluIiwic2V0UHJvamVjdGlvbiIsImNvbmZpZyIsInNldExheWVycyIsImdldExheWVySW5kZXgiLCJ1cGRhdGVMYXllcnMiLCJwcm9qZWN0aW9uIiwiZ2V0V29ya2VyU291cmNlIiwibG9hZERFTVRpbGUiLCJnZXRERU1Xb3JrZXJTb3VyY2UiLCJhYm9ydFRpbGUiLCJyZW1vdmVUaWxlIiwicmVtb3ZlU291cmNlIiwid29ya2VyIiwibG9hZFdvcmtlclNvdXJjZSIsImltcG9ydFNjcmlwdHMiLCJ1cmwiLCJ0b1N0cmluZyIsInN5bmNSVExQbHVnaW5TdGF0ZSIsInN0YXRlIiwic2V0U3RhdGUiLCJwbHVnaW5VUkwiLCJnZXRQbHVnaW5VUkwiLCJpc0xvYWRlZCIsImNvbXBsZXRlIiwiZXJyb3IiLCJnZXRBdmFpbGFibGVJbWFnZXMiLCJzZW5kIiwiXyIsIm11c3RRdWV1ZSIsIm1ldGFkYXRhIiwic2NoZWR1bGVyIiwiZW5mb3JjZUNhY2hlU2l6ZUxpbWl0IiwiZ2V0V29ya2VyUGVyZm9ybWFuY2VNZXRyaWNzIiwiV29ya2VyR2xvYmFsU2NvcGUiXSwibWFwcGluZ3MiOiI7O0FBSUE7QUFFQSxTQUFTQSxTQUFULENBQW1CQyxHQUFuQixFQUF3QjtBQUFBLElBQ3BCLElBQUksT0FBT0EsR0FBUCxLQUFlLFFBQWYsSUFBMkIsT0FBT0EsR0FBUCxLQUFlLFNBQTFDLElBQXVELE9BQU9BLEdBQVAsS0FBZSxRQUF0RSxJQUFrRkEsR0FBQSxLQUFRQyxTQUExRixJQUF1R0QsR0FBQSxLQUFRLElBQW5IO0FBQUEsUUFDSSxPQUFPRSxJQUFBLENBQUtILFNBQUwsQ0FBZUMsR0FBZixDQUFQLENBRmdCO0FBQUEsSUFJcEIsSUFBSUcsS0FBQSxDQUFNQyxPQUFOLENBQWNKLEdBQWQsQ0FBSixFQUF3QjtBQUFBLFFBQ3BCLElBQUlLLEdBQUEsR0FBTSxHQUFWLENBRG9CO0FBQUEsUUFFcEIsS0FBV0MsTUFBQUEsR0FBWCxJQUFrQk4sR0FBbEIsRUFBdUI7QUFBQSxZQUNuQkssR0FBQSxJQUFPLENBQUEsR0FBR04sU0FBQSxDQUFVTyxHQUFWLENBQUgsRUFBa0IsQ0FBbEIsQ0FBUCxDQURtQjtBQUFBLFNBRkg7QUFBQSxRQUtwQixPQUFPLENBQUdELEdBQUFBLEdBQUgsRUFBTyxDQUFQLENBQVAsQ0FMb0I7QUFBQSxLQUpKO0FBQUEsSUFZcEIsSUFBSUEsR0FBQSxHQUFNLEdBQVYsQ0Fab0I7QUFBQSxJQWFwQixLQUFBLE1BQVdFLEdBQVgsSUFBa0JDLE1BQUEsQ0FBT0MsSUFBUCxDQUFZVCxHQUFaLENBQUEsQ0FBaUJVLElBQWpCLEVBQWxCLEVBQTJDO0FBQUEsUUFDdkNMLEdBQUEsSUFBTyxDQUFBLEdBQUdFLEdBQUgsRUFBTyxDQUFQLEdBQVVSLFNBQUEsQ0FBV0MsR0FBRCxDQUFXTyxHQUFYLENBQVYsQ0FBVixFQUFxQyxDQUFyQyxDQUFQLENBRHVDO0FBQUEsS0FidkI7QUFBQSxJQWdCcEIsT0FBTyxDQUFHRixHQUFBQSxHQUFILEVBQU8sQ0FBUCxDQUFQLENBaEJvQjtBQUFBLENBTnhCO0FBeUJBLFNBQVNNLE1BQVQsQ0FBZ0JDLEtBQWhCLEVBQXVCO0FBQUEsSUFDbkIsSUFBSUwsR0FBQSxHQUFNLEVBQVYsQ0FEbUI7QUFBQSxJQUVuQixLQUFXTSxNQUFBQSxDQUFYLElBQWdCQyw0QkFBaEIsRUFBK0I7QUFBQSxRQUMzQlAsR0FBQSxJQUFPLENBQUMsQ0FBRCxHQUFJUixTQUFBLENBQVdhLEtBQUQsQ0FBYUMsQ0FBYixDQUFWLENBQUosR0FBUCxDQUQyQjtBQUFBLEtBRlo7QUFBQSxJQUtuQixPQUFPTixHQUFQLENBTG1CO0FBQUEsQ0F6QnZCO0FBZ0RBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFlLFNBQVNRLGFBQVQsQ0FBdUJDLE1BQXZCLEVBQTBEQyxVQUExRCxFQUFnSTtBQUFBLElBQzNJLE1BQU1DLE1BQUEsR0FBUyxFQUFmLENBRDJJO0FBQUEsSUFHM0ksS0FBSyxJQUFJQyxDQUFBLEdBQUksQ0FBUixFQUFXQSxDQUFBLEdBQUlILE1BQUEsQ0FBT0ksTUFBM0IsRUFBbUNELENBQUEsRUFBbkMsRUFBd0M7QUFBQSxRQUVwQyxNQUFNTixDQUFBLEdBQUtJLFVBQUEsSUFBY0EsVUFBQSxDQUFXRCxNQUFBLENBQU9HLENBQVAsRUFBVUUsRUFBckIsQ0FBZixJQUE0Q1YsTUFBQSxDQUFPSyxNQUFBLENBQU9HLENBQVAsQ0FBUCxDQUF0RCxDQUZvQztBQUFBO0FBSXBDLFFBQUEsSUFBSUYsVUFBSjtBQUFBLFlBQ0lBLFVBQUEsQ0FBV0QsTUFBQSxDQUFPRyxDQUFQLENBQVVFLENBQUFBLEVBQXJCLENBQTJCUixHQUFBQSxDQUEzQixDQUxnQztBQUFBLFFBT3BDLElBQUlTLEtBQUEsR0FBUUosTUFBQSxDQUFPTCxDQUFQLENBQVosQ0FQb0M7QUFBQSxRQVFwQyxJQUFJLENBQUNTLEtBQUwsRUFBWTtBQUFBLFlBQ1JBLEtBQUEsR0FBUUosTUFBQSxDQUFPTCxDQUFQLENBQUEsR0FBWSxFQUFwQixDQURRO0FBQUEsU0FSd0I7QUFBQSxRQVdwQ1MsS0FBQSxDQUFNQyxJQUFOLENBQVdQLE1BQUEsQ0FBT0csQ0FBUCxDQUFYLENBWG9DLENBQUE7QUFBQSxLQUhtRztBQUFBLElBaUIzSSxNQUFNSyxNQUFBLEdBQVMsRUFBZixDQWpCMkk7QUFBQSxJQW1CM0ksS0FBV1gsTUFBQUEsQ0FBWCxJQUFnQkssTUFBaEIsRUFBd0I7QUFBQSxRQUNwQk0sTUFBQSxDQUFPRCxJQUFQLENBQVlMLE1BQUEsQ0FBT0wsQ0FBUCxDQUFaLENBRG9CLENBQUE7QUFBQSxLQW5CbUg7QUFBQSxJQXVCM0ksT0FBT1csTUFBUCxDQXZCMkk7QUFBQTs7QUM5Qy9JO0FBV0EsTUFBTUMsZUFBTixDQUFzQjtBQUFBLElBT2xCQyxXQUFBLENBQVlDLFlBQVosRUFBc0Q7QUFBQSxRQUNsRCxJQUFLQyxDQUFBQSxRQUFMLEdBQWdCLEVBQWhCLENBRGtEO0FBQUEsUUFFbEQsSUFBSUQsWUFBSixFQUFrQjtBQUFBLFlBQ2QsSUFBS0UsQ0FBQUEsT0FBTCxDQUFhRixZQUFiLENBRGMsQ0FBQTtBQUFBLFNBRmdDO0FBQUEsS0FQcEM7QUFBQSxJQWNsQkUsT0FBQSxDQUFRRixZQUFSLEVBQWlEO0FBQUEsUUFDN0MsSUFBS0csQ0FBQUEsYUFBTCxHQUFxQixFQUFyQixDQUQ2QztBQUFBLFFBRTdDLElBQUtDLENBQUFBLE9BQUwsR0FBZSxFQUFmLENBRjZDO0FBQUEsUUFHN0MsS0FBS0MsTUFBTCxDQUFZTCxZQUFaLEVBQTBCLEVBQTFCLENBSDZDLENBQUE7QUFBQSxLQWQvQjtBQUFBLElBb0JsQkssTUFBQSxDQUFPTCxZQUFQLEVBQWdETSxVQUFoRCxFQUEyRTtBQUFBLFFBQ3ZFLEtBQVdDLE1BQUFBLFdBQVgsSUFBMEJQLFlBQTFCLEVBQXdDO0FBQUEsWUFDcEMsSUFBQSxDQUFLRyxhQUFMLENBQW1CSSxXQUFBLENBQVliLEVBQS9CLENBQUEsR0FBcUNhLFdBQXJDLENBRG9DO0FBQUEsWUFHcEMsTUFBTXRCLEtBQUEsR0FBUSxJQUFBLENBQUttQixPQUFMLENBQWFHLFdBQUEsQ0FBWWIsRUFBekIsQ0FBaUNjLEdBQUFBLCtCQUFBLENBQWlCRCxXQUFqQixDQUEvQyxDQUhvQztBQUFBLFlBSXBDdEIsS0FBQSxDQUFNd0IsYUFBTixFQUpvQyxDQUFBO0FBQUEsWUFLcEMsSUFBSSxJQUFLUixDQUFBQSxRQUFMLENBQWNNLFdBQUEsQ0FBWWIsRUFBMUIsQ0FBSjtBQUFBLGdCQUNJLE9BQU8sS0FBS08sUUFBTCxDQUFjTSxXQUFBLENBQVliLEVBQTFCLENBQVAsQ0FOZ0M7QUFBQSxTQUQrQjtBQUFBLFFBU3ZFLEtBQVdBLE1BQUFBLEVBQVgsSUFBaUJZLFVBQWpCLEVBQTZCO0FBQUEsWUFDekIsT0FBTyxJQUFLTCxDQUFBQSxRQUFMLENBQWNQLEVBQWQsQ0FBUCxDQUR5QjtBQUFBLFlBRXpCLE9BQU8sSUFBS1MsQ0FBQUEsYUFBTCxDQUFtQlQsRUFBbkIsQ0FBUCxDQUZ5QjtBQUFBLFlBR3pCLE9BQU8sSUFBS1UsQ0FBQUEsT0FBTCxDQUFhVixFQUFiLENBQVAsQ0FIeUI7QUFBQSxTQVQwQztBQUFBLFFBZXZFLElBQUtnQixDQUFBQSxnQkFBTCxHQUF3QixFQUF4QixDQWZ1RTtBQUFBLFFBaUJ2RSxNQUFNbkIsTUFBQSxHQUFTSCxhQUFBLENBQWN1QixxQkFBQSxDQUFPLElBQUtSLENBQUFBLGFBQVosQ0FBZCxFQUEwQyxJQUFLRixDQUFBQSxRQUEvQyxDQUFmLENBakJ1RTtBQUFBLFFBbUJ2RSxLQUFXRCxNQUFBQSxZQUFYLElBQTJCVCxNQUEzQixFQUFtQztBQUFBLFlBQy9CLE1BQU1GLE1BQUEsR0FBU1csWUFBQSxDQUFhWSxHQUFiLENBQWtCTCxXQUFELElBQWlCLElBQUEsQ0FBS0gsT0FBTCxDQUFhRyxXQUFBLENBQVliLEVBQXpCLENBQWxDLENBQWYsQ0FEK0I7QUFBQSxZQUcvQixNQUFNVCxLQUFBLEdBQVFJLE1BQUEsQ0FBTyxDQUFQLENBQWQsQ0FIK0I7QUFBQSxZQUkvQixJQUFJSixLQUFBLENBQU00QixVQUFOLEtBQXFCLE1BQXpCLEVBQWlDO0FBQUEsZ0JBQzdCLFNBRDZCO0FBQUEsYUFKRjtBQUFBLFlBUS9CLE1BQU1DLFFBQUEsR0FBVzdCLEtBQUEsQ0FBTThCLE1BQU4sSUFBZ0IsRUFBakMsQ0FSK0I7QUFBQSxZQVMvQixJQUFJQyxXQUFBLEdBQWMsSUFBQSxDQUFLTixnQkFBTCxDQUFzQkksUUFBdEIsQ0FBbEIsQ0FUK0I7QUFBQSxZQVUvQixJQUFJLENBQUNFLFdBQUwsRUFBa0I7QUFBQSxnQkFDZEEsV0FBQSxHQUFjLElBQUtOLENBQUFBLGdCQUFMLENBQXNCSSxRQUF0QixDQUFBLEdBQWtDLEVBQWhELENBRGM7QUFBQSxhQVZhO0FBQUEsWUFjL0IsTUFBTUcsYUFBQSxHQUFnQmhDLEtBQUEsQ0FBTWlDLFdBQU4sSUFBcUIsbUJBQTNDLENBZCtCO0FBQUEsWUFlL0IsSUFBSUMsbUJBQUEsR0FBc0JILFdBQUEsQ0FBWUMsYUFBWixDQUExQixDQWYrQjtBQUFBLFlBZ0IvQixJQUFJLENBQUNFLG1CQUFMLEVBQTBCO0FBQUEsZ0JBQ3RCQSxtQkFBQSxHQUFzQkgsV0FBQSxDQUFZQyxhQUFaLENBQUEsR0FBNkIsRUFBbkQsQ0FEc0I7QUFBQSxhQWhCSztBQUFBLFlBb0IvQkUsbUJBQUEsQ0FBb0J2QixJQUFwQixDQUF5QlAsTUFBekIsQ0FwQitCLENBQUE7QUFBQSxTQW5Cb0M7QUFBQSxLQXBCekQ7QUFBQTs7QUNYdEI7QUFNQSxNQUFNK0IseUJBQU4sQ0FBZ0M7QUFBQSxJQUs1QkMsUUFBQSxDQUFTQyxNQUFULEVBQTBDQyxRQUExQyxFQUEyRTtBQUFBLFFBQ3ZFLE1BQU0sQ0FBQ0MsR0FBRCxFQUFNQyxRQUFOLEVBQWdCQyxZQUFoQixFQUE4QkMsT0FBOUIsRUFBdUNDLGFBQXZDLENBQUEsR0FBd0ROLE1BQTlELENBRHVFO0FBQUE7QUFJdkU7QUFBQSxRQUFNTyxNQUFBQSxXQUFBLEdBQWNDLHFCQUFBLENBQU9DLFdBQVAsSUFBc0JMLFlBQUEsWUFBd0JJLHFCQUFBLENBQU9DLFdBQXJELEdBQW1FLElBQUEsQ0FBS0MsWUFBTCxDQUFrQk4sWUFBbEIsRUFBZ0NDLE9BQWhDLENBQW5FLEdBQWdIRCxZQUFwSSxDQUp1RTtBQUFBLFFBS3ZFLE1BQU1PLEdBQUEsR0FBTSxJQUFJQyxzQkFBSixDQUFZVixHQUFaLEVBQWlCSyxXQUFqQixFQUE4QkosUUFBOUIsRUFBd0NFLE9BQUEsR0FBVSxDQUFsRCxFQUFxREMsYUFBckQsQ0FBWixDQUx1RTtBQUFBLFFBTXZFTCxRQUFBLENBQVMsSUFBVCxFQUFlVSxHQUFmLENBTnVFLENBQUE7QUFBQSxLQUwvQztBQUFBLElBYzVCRCxZQUFBLENBQWFHLFNBQWIsRUFBcUNSLE9BQXJDLEVBQWlFO0FBQUE7QUFFN0QsUUFBQSxJQUFJLENBQUMsSUFBS1MsQ0FBQUEsZUFBTixJQUF5QixDQUFDLElBQUEsQ0FBS0Msc0JBQW5DLEVBQTJEO0FBQUE7QUFFdkQsWUFBQSxJQUFBLENBQUtELGVBQUwsR0FBdUIsSUFBSUUsZUFBSixDQUFvQkgsU0FBQSxDQUFVSSxLQUE5QixFQUFxQ0osU0FBQSxDQUFVSyxNQUEvQyxDQUF2QixDQUZ1RDtBQUFBLFlBR3ZELElBQUEsQ0FBS0gsc0JBQUwsR0FBOEIsSUFBS0QsQ0FBQUEsZUFBTCxDQUFxQkssVUFBckIsQ0FBZ0MsSUFBaEMsQ0FBOUIsQ0FIdUQ7QUFBQSxTQUZFO0FBQUEsUUFRN0QsSUFBQSxDQUFLTCxlQUFMLENBQXFCRyxLQUFyQixHQUE2QkosU0FBQSxDQUFVSSxLQUF2QyxDQVI2RDtBQUFBLFFBUzdELElBQUEsQ0FBS0gsZUFBTCxDQUFxQkksTUFBckIsR0FBOEJMLFNBQUEsQ0FBVUssTUFBeEMsQ0FUNkQ7QUFBQSxRQVc3RCxLQUFLSCxzQkFBTCxDQUE0QkssU0FBNUIsQ0FBc0NQLFNBQXRDLEVBQWlELENBQWpELEVBQW9ELENBQXBELEVBQXVEQSxTQUFBLENBQVVJLEtBQWpFLEVBQXdFSixTQUFBLENBQVVLLE1BQWxGLENBWDZELENBQUE7QUFBQTtBQWE3RCxRQUFNRyxNQUFBQSxPQUFBLEdBQVUsSUFBQSxDQUFLTixzQkFBTCxDQUE0QkwsWUFBNUIsQ0FBeUMsQ0FBQ0wsT0FBMUMsRUFBbUQsQ0FBQ0EsT0FBcEQsRUFBNkRRLFNBQUEsQ0FBVUksS0FBVixHQUFrQixDQUFJWixHQUFBQSxPQUFuRixFQUE0RlEsU0FBQSxDQUFVSyxNQUFWLEdBQW1CLENBQUEsR0FBSWIsT0FBbkgsQ0FBaEIsQ0FiNkQ7QUFBQSxRQWM3RCxLQUFLVSxzQkFBTCxDQUE0Qk8sU0FBNUIsQ0FBc0MsQ0FBdEMsRUFBeUMsQ0FBekMsRUFBNEMsSUFBS1IsQ0FBQUEsZUFBTCxDQUFxQkcsS0FBakUsRUFBd0UsS0FBS0gsZUFBTCxDQUFxQkksTUFBN0YsQ0FkNkQsQ0FBQTtBQUFBLFFBZTdELE9BQU9HLE9BQVAsQ0FmNkQ7QUFBQSxLQWRyQztBQUFBOztJQ1BoQ0UsYUFBQSxHQUFpQkMsUUFBakIsQ0FEQTtBQUdBLFNBQVNBLFFBQVQsQ0FBZ0JDLEVBQWhCLEVBQW9CQyxLQUFwQixFQUEyQjtBQUFBLElBQ3ZCLElBQUlDLElBQUEsR0FBT0YsRUFBQSxJQUFNQSxFQUFBLENBQUdFLElBQXBCLEVBQTBCekQsQ0FBMUIsQ0FEdUI7QUFBQSxJQUd2QixJQUFJeUQsSUFBQSxLQUFTLG1CQUFiLEVBQWtDO0FBQUEsUUFDOUIsS0FBS3pELENBQUEsR0FBSSxDQUFULEVBQVlBLENBQUEsR0FBSXVELEVBQUEsQ0FBR0csUUFBSCxDQUFZekQsTUFBNUIsRUFBb0NELENBQUEsRUFBcEM7QUFBQSxZQUF5Q3NELFFBQUEsQ0FBT0MsRUFBQSxDQUFHRyxRQUFILENBQVkxRCxDQUFaLENBQVAsRUFBdUJ3RCxLQUF2QixDQURYLENBQUE7QUFBQSxLQUFsQyxNQUdPLElBQUlDLElBQUEsS0FBUyxvQkFBYixFQUFtQztBQUFBLFFBQ3RDLEtBQUt6RCxDQUFBLEdBQUksQ0FBVCxFQUFZQSxDQUFBLEdBQUl1RCxFQUFBLENBQUdJLFVBQUgsQ0FBYzFELE1BQTlCLEVBQXNDRCxDQUFBLEVBQXRDO0FBQUEsWUFBMkNzRCxRQUFBLENBQU9DLEVBQUEsQ0FBR0ksVUFBSCxDQUFjM0QsQ0FBZCxDQUFQLEVBQXlCd0QsS0FBekIsQ0FETCxDQUFBO0FBQUEsS0FBbkMsTUFHQSxJQUFJQyxJQUFBLEtBQVMsU0FBYixFQUF3QjtBQUFBLFFBQzNCSCxRQUFBLENBQU9DLEVBQUEsQ0FBR0ssUUFBVixFQUFvQkosS0FBcEIsQ0FEMkIsQ0FBQTtBQUFBLEtBQXhCLE1BR0EsSUFBSUMsSUFBQSxLQUFTLFNBQWIsRUFBd0I7QUFBQSxRQUMzQkksV0FBQSxDQUFZTixFQUFBLENBQUdPLFdBQWYsRUFBNEJOLEtBQTVCLENBRDJCLENBQUE7QUFBQSxLQUF4QixNQUdBLElBQUlDLElBQUEsS0FBUyxjQUFiLEVBQTZCO0FBQUEsUUFDaEMsS0FBS3pELENBQUEsR0FBSSxDQUFULEVBQVlBLENBQUEsR0FBSXVELEVBQUEsQ0FBR08sV0FBSCxDQUFlN0QsTUFBL0IsRUFBdUNELENBQUEsRUFBdkM7QUFBQSxZQUE0QzZELFdBQUEsQ0FBWU4sRUFBQSxDQUFHTyxXQUFILENBQWU5RCxDQUFmLENBQVosRUFBK0J3RCxLQUEvQixDQURaLENBQUE7QUFBQSxLQWZiO0FBQUEsSUFtQnZCLE9BQU9ELEVBQVAsQ0FuQnVCO0FBQUEsQ0FIM0I7QUF5QkEsU0FBU00sV0FBVCxDQUFxQkUsS0FBckIsRUFBNEJQLEtBQTVCLEVBQW1DO0FBQUEsSUFDL0IsSUFBSU8sS0FBQSxDQUFNOUQsTUFBTixLQUFpQixDQUFyQjtBQUFBLFFBQXdCLE9BRE87QUFBQSxJQUcvQitELFVBQUEsQ0FBV0QsS0FBQSxDQUFNLENBQU4sQ0FBWCxFQUFxQlAsS0FBckIsQ0FIK0IsQ0FBQTtBQUFBLElBSS9CLEtBQUssSUFBSXhELENBQUEsR0FBSSxDQUFSLEVBQVdBLENBQUEsR0FBSStELEtBQUEsQ0FBTTlELE1BQTFCLEVBQWtDRCxDQUFBLEVBQWxDLEVBQXVDO0FBQUEsUUFDbkNnRSxVQUFBLENBQVdELEtBQUEsQ0FBTS9ELENBQU4sQ0FBWCxFQUFxQixDQUFDd0QsS0FBdEIsQ0FEbUMsQ0FBQTtBQUFBLEtBSlI7QUFBQSxDQXpCbkM7QUFrQ0EsU0FBU1EsVUFBVCxDQUFvQkMsSUFBcEIsRUFBMEJDLEdBQTFCLEVBQStCO0FBQUEsSUFDM0IsSUFBSUMsSUFBQSxHQUFPLENBQVgsRUFBY0MsR0FBQSxHQUFNLENBQXBCLENBRDJCO0FBQUEsSUFFM0IsS0FBSyxJQUFJcEUsQ0FBQSxHQUFJLENBQVIsRUFBV3FFLEdBQUEsR0FBTUosSUFBQSxDQUFLaEUsTUFBdEIsRUFBOEJxRSxDQUFBLEdBQUlELEdBQUEsR0FBTSxDQUF4QyxFQUEyQ3JFLENBQUEsR0FBSXFFLEdBQXBELEVBQXlEQyxDQUFBLEdBQUl0RSxDQUFBLEVBQTdELEVBQWtFO0FBQUEsUUFDOUQsSUFBSU4sQ0FBQSxHQUFLLENBQUF1RSxJQUFBLENBQUtqRSxDQUFMLENBQVEsQ0FBQSxDQUFSLENBQWFpRSxHQUFBQSxJQUFBLENBQUtLLENBQUwsRUFBUSxDQUFSLENBQWIsS0FBNEJMLElBQUEsQ0FBS0ssQ0FBTCxFQUFRLENBQVIsQ0FBQSxHQUFhTCxJQUFBLENBQUtqRSxDQUFMLENBQUEsQ0FBUSxDQUFSLENBQWIsQ0FBckMsQ0FEOEQ7QUFBQSxRQUU5RCxJQUFJdUUsQ0FBQSxHQUFJSixJQUFBLEdBQU96RSxDQUFmLENBRjhEO0FBQUEsUUFHOUQwRSxHQUFBLElBQU9JLElBQUEsQ0FBS0MsR0FBTCxDQUFTTixJQUFULENBQWtCSyxJQUFBQSxJQUFBLENBQUtDLEdBQUwsQ0FBUy9FLENBQVQsQ0FBbEIsR0FBZ0N5RSxJQUFBLEdBQU9JLENBQVAsR0FBVzdFLENBQTNDLEdBQStDQSxDQUFBLEdBQUk2RSxDQUFKLEdBQVFKLElBQTlELENBSDhEO0FBQUEsUUFJOURBLElBQUEsR0FBT0ksQ0FBUCxDQUo4RDtBQUFBLEtBRnZDO0FBQUEsSUFRM0IsSUFBSUosSUFBQSxHQUFPQyxHQUFQLElBQWMsQ0FBZCxLQUFvQixDQUFDLENBQUNGLEdBQTFCO0FBQUEsUUFBK0JELElBQUEsQ0FBS1MsT0FBTCxFQVJKLENBQUE7QUFBQTs7QUNoQy9CO0FBR0EsTUFBTUMsU0FBQSxHQUFZQyx5QkFBQSxDQUFJQyxpQkFBSixDQUFzQkMsU0FBdEIsQ0FBZ0NILFNBQWxELENBTEE7QUFzQkE7QUFBQTtBQUFBLE1BQU1JLGdCQUFOLENBQWtEO0FBQUEsSUFROUN4RSxXQUFBLENBQVl5RSxPQUFaLEVBQThCO0FBQUEsUUFDMUIsSUFBS0MsQ0FBQUEsUUFBTCxHQUFnQkQsT0FBaEIsQ0FEMEI7QUFBQSxRQUcxQixJQUFLRSxDQUFBQSxNQUFMLEdBQWNDLHFCQUFkLENBSDBCO0FBQUEsUUFJMUIsS0FBSzFCLElBQUwsR0FBWXVCLE9BQUEsQ0FBUXZCLElBQXBCLENBSjBCO0FBQUEsUUFLMUIsS0FBSzJCLFVBQUwsR0FBa0JKLE9BQUEsQ0FBUUssSUFBMUIsQ0FMMEI7QUFBQTtBQWExQjtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUEsUUFBSSxJQUFBLElBQUEsSUFBUUwsT0FBUixJQUFtQixDQUFDTSxLQUFBLENBQU1OLE9BQUEsQ0FBUTlFLEVBQWQsQ0FBeEIsRUFBMkM7QUFBQSxZQUN2QyxJQUFBLENBQUtBLEVBQUwsR0FBVXFGLFFBQUEsQ0FBU1AsT0FBQSxDQUFROUUsRUFBakIsRUFBcUIsRUFBckIsQ0FBVixDQUR1QztBQUFBLFNBYmpCO0FBQUEsS0FSZ0I7QUFBQSxJQTBCOUNzRixZQUFBLEdBQWU7QUFBQSxRQUNYLElBQUksSUFBS1AsQ0FBQUEsUUFBTCxDQUFjeEIsSUFBZCxLQUF1QixDQUEzQixFQUE4QjtBQUFBLFlBQzFCLE1BQU1HLFFBQUEsR0FBVyxFQUFqQixDQUQwQjtBQUFBLFlBRTFCLFdBQVc2QixLQUFYLElBQW9CLEtBQUtSLFFBQUwsQ0FBY3JCLFFBQWxDLEVBQTRDO0FBQUEsZ0JBQ3hDQSxRQUFBLENBQVN4RCxJQUFULENBQWMsQ0FBQyxJQUFJc0YsNEJBQUosQ0FBVUQsS0FBQSxDQUFNLENBQU4sQ0FBVixFQUFvQkEsS0FBQSxDQUFNLENBQU4sQ0FBcEIsQ0FBRCxDQUFkLENBRHdDLENBQUE7QUFBQSxhQUZsQjtBQUFBLFlBSzFCLE9BQU83QixRQUFQLENBTDBCO0FBQUEsU0FBOUIsTUFNTztBQUFBLFlBQ0gsTUFBTUEsUUFBQSxHQUFXLEVBQWpCLENBREc7QUFBQSxZQUVILFdBQVdLLElBQVgsSUFBbUIsS0FBS2dCLFFBQUwsQ0FBY3JCLFFBQWpDLEVBQTJDO0FBQUEsZ0JBQ3ZDLE1BQU0rQixPQUFBLEdBQVUsRUFBaEIsQ0FEdUM7QUFBQSxnQkFFdkMsS0FBV0YsTUFBQUEsS0FBWCxJQUFvQnhCLElBQXBCLEVBQTBCO0FBQUEsb0JBQ3RCMEIsT0FBQSxDQUFRdkYsSUFBUixDQUFhLElBQUlzRiw0QkFBSixDQUFVRCxLQUFBLENBQU0sQ0FBTixDQUFWLEVBQW9CQSxLQUFBLENBQU0sQ0FBTixDQUFwQixDQUFiLENBRHNCLENBQUE7QUFBQSxpQkFGYTtBQUFBLGdCQUt2QzdCLFFBQUEsQ0FBU3hELElBQVQsQ0FBY3VGLE9BQWQsQ0FMdUMsQ0FBQTtBQUFBLGFBRnhDO0FBQUEsWUFTSCxPQUFPL0IsUUFBUCxDQVRHO0FBQUEsU0FQSTtBQUFBLEtBMUIrQjtBQUFBLElBOEM5Q2UsU0FBQSxDQUFVaUIsQ0FBVixFQUFxQkMsQ0FBckIsRUFBZ0NDLENBQWhDLEVBQTJDO0FBQUEsUUFDdkMsT0FBT25CLFNBQUEsQ0FBVW9CLElBQVYsQ0FBZSxJQUFmLEVBQXFCSCxDQUFyQixFQUF3QkMsQ0FBeEIsRUFBMkJDLENBQTNCLENBQVAsQ0FEdUM7QUFBQSxLQTlDRztBQUFBLENBdEJsRDtBQXlFQSxNQUFNRSxnQkFBTixDQUE0RDtBQUFBLElBT3hEekYsV0FBQSxDQUFZbUQsUUFBWixFQUFzQztBQUFBLFFBQ2xDLEtBQUs3RCxNQUFMLEdBQWMsRUFBQyxtQkFBcUIsRUFBQSxJQUF0QixFQUFkLENBRGtDO0FBQUEsUUFFbEMsSUFBS29HLENBQUFBLElBQUwsR0FBWSxtQkFBWixDQUZrQztBQUFBLFFBR2xDLElBQUtmLENBQUFBLE1BQUwsR0FBY0MscUJBQWQsQ0FIa0M7QUFBQSxRQUlsQyxLQUFLbEYsTUFBTCxHQUFjeUQsUUFBQSxDQUFTekQsTUFBdkIsQ0FKa0M7QUFBQSxRQUtsQyxJQUFLaUcsQ0FBQUEsU0FBTCxHQUFpQnhDLFFBQWpCLENBTGtDO0FBQUEsS0FQa0I7QUFBQSxJQWV4RHNCLE9BQUEsQ0FBUWhGLENBQVIsRUFBc0M7QUFBQSxRQUNsQyxPQUFPLElBQUkrRSxnQkFBSixDQUFtQixJQUFBLENBQUttQixTQUFMLENBQWVsRyxDQUFmLENBQW5CLENBQVAsQ0FEa0M7QUFBQSxLQWZrQjtBQUFBOzs7O0FDdkU1RCxJQUFJMEYsS0FBQSxHQUFRUyw0QkFBWixDQUZBO0FBR0EsSUFBSXRCLGlCQUFBLEdBQW9Cc0IseUJBQUEsQ0FBK0J0QixpQkFBdkQsQ0FIQTtJQUtBeEIsZUFBQSxHQUFpQjJDLGdCQUFqQixDQUxBO0FBUUE7QUFBQSxTQUFTQSxnQkFBVCxDQUF5QnRDLFFBQXpCLEVBQW1DMEMsT0FBbkMsRUFBNEM7QUFBQSxJQUMxQyxLQUFLQSxPQUFMLEdBQWVBLE9BQUEsSUFBVyxFQUExQixDQUQwQztBQUFBLElBRTFDLElBQUsxQyxDQUFBQSxRQUFMLEdBQWdCQSxRQUFoQixDQUYwQztBQUFBLElBRzFDLEtBQUt6RCxNQUFMLEdBQWN5RCxRQUFBLENBQVN6RCxNQUF2QixDQUgwQztBQUFBLENBUjVDO0FBY0ErRixnQkFBQSxDQUFlbEIsU0FBZixDQUF5QkUsT0FBekIsR0FBbUMsVUFBVWhGLENBQVYsRUFBYTtBQUFBLElBQzlDLE9BQU8sSUFBSStFLGNBQUosQ0FBbUIsS0FBS3JCLFFBQUwsQ0FBYzFELENBQWQsQ0FBbkIsRUFBcUMsSUFBS29HLENBQUFBLE9BQUwsQ0FBYWxCLE1BQWxELENBQVAsQ0FEOEM7QUFBQSxDQUFoRCxDQWRBO0FBa0JBLFNBQVNILGNBQVQsQ0FBeUJDLE9BQXpCLEVBQWtDRSxNQUFsQyxFQUEwQztBQUFBLElBQ3hDLElBQUtoRixDQUFBQSxFQUFMLEdBQVUsT0FBTzhFLE9BQUEsQ0FBUTlFLEVBQWYsS0FBc0IsUUFBdEIsR0FBaUM4RSxPQUFBLENBQVE5RSxFQUF6QyxHQUE4Q3BCLFNBQXhELENBRHdDO0FBQUEsSUFFeEMsS0FBSzJFLElBQUwsR0FBWXVCLE9BQUEsQ0FBUXZCLElBQXBCLENBRndDO0FBQUEsSUFHeEMsS0FBSzRDLFdBQUwsR0FBbUJyQixPQUFBLENBQVF2QixJQUFSLEtBQWlCLENBQWpCLEdBQXFCLENBQUN1QixPQUFBLENBQVFwQixRQUFULENBQXJCLEdBQTBDb0IsT0FBQSxDQUFRcEIsUUFBckUsQ0FId0M7QUFBQSxJQUl4QyxLQUFLd0IsVUFBTCxHQUFrQkosT0FBQSxDQUFRSyxJQUExQixDQUp3QztBQUFBLElBS3hDLEtBQUtILE1BQUwsR0FBY0EsTUFBQSxJQUFVLElBQXhCLENBTHdDO0FBQUEsQ0FsQjFDO0FBMEJBSCxjQUFBLENBQWVELFNBQWYsQ0FBeUJVLFlBQXpCLEdBQXdDLFlBQVk7QUFBQSxJQUNsRCxJQUFJekIsS0FBQSxHQUFRLElBQUEsQ0FBS3NDLFdBQWpCLENBRGtEO0FBQUEsSUFFbEQsSUFBS3pDLENBQUFBLFFBQUwsR0FBZ0IsRUFBaEIsQ0FGa0Q7QUFBQSxJQUlsRCxLQUFLLElBQUk1RCxDQUFBLEdBQUksQ0FBUixFQUFXQSxDQUFBLEdBQUkrRCxLQUFBLENBQU05RCxNQUExQixFQUFrQ0QsQ0FBQSxFQUFsQyxFQUF1QztBQUFBLFFBQ3JDLElBQUlpRSxJQUFBLEdBQU9GLEtBQUEsQ0FBTS9ELENBQU4sQ0FBWCxDQURxQztBQUFBLFFBRXJDLElBQUkyRixPQUFBLEdBQVUsRUFBZCxDQUZxQztBQUFBLFFBR3JDLEtBQUssSUFBSXJCLENBQUEsR0FBSSxDQUFSLEVBQVdBLENBQUEsR0FBSUwsSUFBQSxDQUFLaEUsTUFBekIsRUFBaUNxRSxDQUFBLEVBQWpDLEVBQXNDO0FBQUEsWUFDcENxQixPQUFBLENBQVF2RixJQUFSLENBQWEsSUFBSXNGLEtBQUosQ0FBVXpCLElBQUEsQ0FBS0ssQ0FBTCxFQUFRLENBQVIsQ0FBVixFQUFzQkwsSUFBQSxDQUFLSyxDQUFMLENBQVEsQ0FBQSxDQUFSLENBQXRCLENBQWIsQ0FEb0MsQ0FBQTtBQUFBLFNBSEQ7QUFBQSxRQU1yQyxLQUFLVixRQUFMLENBQWN4RCxJQUFkLENBQW1CdUYsT0FBbkIsQ0FOcUMsQ0FBQTtBQUFBLEtBSlc7QUFBQSxJQVlsRCxPQUFPLElBQUsvQixDQUFBQSxRQUFaLENBWmtEO0FBQUEsQ0FBcEQsQ0ExQkE7QUF5Q0FtQixjQUFBLENBQWVELFNBQWYsQ0FBeUJ3QixJQUF6QixHQUFnQyxZQUFZO0FBQUEsSUFDMUMsSUFBSSxDQUFDLElBQUEsQ0FBSzFDLFFBQVY7QUFBQSxRQUFvQixJQUFBLENBQUs0QixZQUFMLEVBRHNCLENBQUE7QUFBQSxJQUcxQyxJQUFJekIsS0FBQSxHQUFRLElBQUEsQ0FBS0gsUUFBakIsQ0FIMEM7QUFBQSxJQUkxQyxJQUFJMkMsRUFBQSxHQUFLQyxRQUFULENBSjBDO0FBQUEsSUFLMUMsSUFBSUMsRUFBQSxHQUFLLENBQUNELFFBQVYsQ0FMMEM7QUFBQSxJQU0xQyxJQUFJRSxFQUFBLEdBQUtGLFFBQVQsQ0FOMEM7QUFBQSxJQU8xQyxJQUFJRyxFQUFBLEdBQUssQ0FBQ0gsUUFBVixDQVAwQztBQUFBLElBUzFDLEtBQUssSUFBSXhHLENBQUEsR0FBSSxDQUFSLEVBQVdBLENBQUEsR0FBSStELEtBQUEsQ0FBTTlELE1BQTFCLEVBQWtDRCxDQUFBLEVBQWxDLEVBQXVDO0FBQUEsUUFDckMsSUFBSWlFLElBQUEsR0FBT0YsS0FBQSxDQUFNL0QsQ0FBTixDQUFYLENBRHFDO0FBQUEsUUFHckMsS0FBSyxJQUFJc0UsQ0FBQSxHQUFJLENBQVIsRUFBV0EsQ0FBQSxHQUFJTCxJQUFBLENBQUtoRSxNQUF6QixFQUFpQ3FFLENBQUEsRUFBakMsRUFBc0M7QUFBQSxZQUNwQyxJQUFJc0MsS0FBQSxHQUFRM0MsSUFBQSxDQUFLSyxDQUFMLENBQVosQ0FEb0M7QUFBQSxZQUdwQ2lDLEVBQUEsR0FBSy9CLElBQUEsQ0FBS3FDLEdBQUwsQ0FBU04sRUFBVCxFQUFhSyxLQUFBLENBQU1oQixDQUFuQixDQUFMLENBSG9DO0FBQUEsWUFJcENhLEVBQUEsR0FBS2pDLElBQUEsQ0FBS3NDLEdBQUwsQ0FBU0wsRUFBVCxFQUFhRyxLQUFBLENBQU1oQixDQUFuQixDQUFMLENBSm9DO0FBQUEsWUFLcENjLEVBQUEsR0FBS2xDLElBQUEsQ0FBS3FDLEdBQUwsQ0FBU0gsRUFBVCxFQUFhRSxLQUFBLENBQU1mLENBQW5CLENBQUwsQ0FMb0M7QUFBQSxZQU1wQ2MsRUFBQSxHQUFLbkMsSUFBQSxDQUFLc0MsR0FBTCxDQUFTSCxFQUFULEVBQWFDLEtBQUEsQ0FBTWYsQ0FBbkIsQ0FBTCxDQU5vQztBQUFBLFNBSEQ7QUFBQSxLQVRHO0FBQUEsSUFzQjFDLE9BQU87QUFBQSxRQUFDVSxFQUFEO0FBQUEsUUFBS0csRUFBTDtBQUFBLFFBQVNELEVBQVQ7QUFBQSxRQUFhRSxFQUFiO0FBQUEsS0FBUCxDQXRCMEM7QUFBQSxDQUE1QyxDQXpDQTtBQWtFQTVCLGNBQUEsQ0FBZUQsU0FBZixDQUF5QkgsU0FBekIsR0FBcUNFLGlCQUFBLENBQWtCQyxTQUFsQixDQUE0QkgsU0FBakU7O0FDbEVBLElBQUlvQyxHQUFBLEdBQU1aLGtCQUFWLENBQUE7QUFDQSxJQUFJSCxjQUFBLEdBQWlCRyxlQUFyQixDQURBO0FBR0E5QyxLQUFBLENBQUEsT0FBQSxHQUFpQjJELGdCQUFqQixDQUhBO0FBSUFDLGFBQUEsQ0FBQSxnQkFBQSxHQUFrQ0QsaUJBSmxDO0FBS0FDLGFBQUEsQ0FBQSxhQUFBLEdBQStCQyxjQUwvQjtBQU1BRCxhQUFBLENBQUEsY0FBQSxHQUFnQ2pCLGVBTmhDO0FBY0E7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUEsU0FBU2dCLGdCQUFULENBQTJCRyxJQUEzQixFQUFpQztBQUFBLElBQy9CLElBQUlDLEdBQUEsR0FBTSxJQUFJTCxHQUFKLEVBQVYsQ0FEK0I7QUFBQSxJQUUvQk0sU0FBQSxDQUFVRixJQUFWLEVBQWdCQyxHQUFoQixDQUYrQixDQUFBO0FBQUEsSUFHL0IsT0FBT0EsR0FBQSxDQUFJRSxNQUFKLEVBQVAsQ0FIK0I7QUFBQSxDQWRqQztBQTZCQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQSxTQUFTSixhQUFULENBQXdCckgsTUFBeEIsRUFBZ0N1RyxPQUFoQyxFQUF5QztBQUFBLElBQ3ZDQSxPQUFBLEdBQVVBLE9BQUEsSUFBVyxFQUFyQixDQUR1QztBQUFBLElBRXZDLElBQUltQixDQUFBLEdBQUksRUFBUixDQUZ1QztBQUFBLElBR3ZDLEtBQVM3SCxJQUFBQSxDQUFULElBQWNHLE1BQWQsRUFBc0I7QUFBQSxRQUNwQjBILENBQUEsQ0FBRTdILENBQUYsQ0FBQSxHQUFPLElBQUlzRyxjQUFKLENBQW1CbkcsTUFBQSxDQUFPSCxDQUFQLENBQVVnRSxDQUFBQSxRQUE3QixFQUF1QzBDLE9BQXZDLENBQVAsQ0FEb0I7QUFBQSxRQUVwQm1CLENBQUEsQ0FBRTdILENBQUYsRUFBS3VHLElBQUwsR0FBWXZHLENBQVosQ0FGb0I7QUFBQSxRQUdwQjZILENBQUEsQ0FBRTdILENBQUYsQ0FBQSxDQUFLOEgsT0FBTCxHQUFlcEIsT0FBQSxDQUFRb0IsT0FBdkIsQ0FIb0I7QUFBQSxRQUlwQkQsQ0FBQSxDQUFFN0gsQ0FBRixDQUFBLENBQUt3RixNQUFMLEdBQWNrQixPQUFBLENBQVFsQixNQUF0QixDQUpvQjtBQUFBLEtBSGlCO0FBQUEsSUFTdkMsT0FBTzhCLGdCQUFBLENBQWlCLEVBQUVuSCxNQUFBLEVBQVEwSCxDQUFWLEVBQWpCLENBQVAsQ0FUdUM7QUFBQSxDQTdCekM7QUF5Q0EsU0FBU0YsU0FBVCxDQUFvQkYsSUFBcEIsRUFBMEJNLEdBQTFCLEVBQStCO0FBQUEsSUFDN0IsU0FBU3JJLEdBQVQsSUFBZ0IrSCxJQUFBLENBQUt0SCxNQUFyQixFQUE2QjtBQUFBLFFBQzNCNEgsR0FBQSxDQUFJQyxZQUFKLENBQWlCLENBQWpCLEVBQW9CQyxVQUFwQixFQUFnQ1IsSUFBQSxDQUFLdEgsTUFBTCxDQUFZVCxHQUFaLENBQWhDLENBRDJCLENBQUE7QUFBQSxLQURBO0FBQUEsQ0F6Qy9CO0FBK0NBLFNBQVN1SSxVQUFULENBQXFCbEksS0FBckIsRUFBNEJnSSxHQUE1QixFQUFpQztBQUFBLElBQy9CQSxHQUFBLENBQUlHLGdCQUFKLENBQXFCLEVBQXJCLEVBQXlCbkksS0FBQSxDQUFNK0gsT0FBTixJQUFpQixDQUExQyxDQUQrQixDQUFBO0FBQUEsSUFFL0JDLEdBQUEsQ0FBSUksZ0JBQUosQ0FBcUIsQ0FBckIsRUFBd0JwSSxLQUFBLENBQU13RyxJQUFOLElBQWMsRUFBdEMsQ0FGK0IsQ0FBQTtBQUFBLElBRy9Cd0IsR0FBQSxDQUFJRyxnQkFBSixDQUFxQixDQUFyQixFQUF3Qm5JLEtBQUEsQ0FBTXlGLE1BQU4sSUFBZ0IsSUFBeEMsQ0FIK0IsQ0FBQTtBQUFBLElBSy9CLElBQUlsRixDQUFKLENBTCtCO0FBQUEsSUFNL0IsSUFBSThILE9BQUEsR0FBVTtBQUFBLFFBQ1p4SSxJQUFBLEVBQU0sRUFETTtBQUFBLFFBRVo2QixNQUFBLEVBQVEsRUFGSTtBQUFBLFFBR1o0RyxRQUFBLEVBQVUsRUFIRTtBQUFBLFFBSVpDLFVBQUEsRUFBWSxFQUpBO0FBQUEsS0FBZCxDQU4rQjtBQUFBLElBYS9CLEtBQUtoSSxDQUFBLEdBQUksQ0FBVCxFQUFZQSxDQUFBLEdBQUlQLEtBQUEsQ0FBTVEsTUFBdEIsRUFBOEJELENBQUEsRUFBOUIsRUFBbUM7QUFBQSxRQUNqQzhILE9BQUEsQ0FBUTlDLE9BQVIsR0FBa0J2RixLQUFBLENBQU11RixPQUFOLENBQWNoRixDQUFkLENBQWxCLENBRGlDO0FBQUEsUUFFakN5SCxHQUFBLENBQUlDLFlBQUosQ0FBaUIsQ0FBakIsRUFBb0JPLFlBQXBCLEVBQWtDSCxPQUFsQyxDQUZpQyxDQUFBO0FBQUEsS0FiSjtBQUFBLElBa0IvQixJQUFJeEksSUFBQSxHQUFPd0ksT0FBQSxDQUFReEksSUFBbkIsQ0FsQitCO0FBQUEsSUFtQi9CLEtBQUtVLENBQUEsR0FBSSxDQUFULEVBQVlBLENBQUEsR0FBSVYsSUFBQSxDQUFLVyxNQUFyQixFQUE2QkQsQ0FBQSxFQUE3QixFQUFrQztBQUFBLFFBQ2hDeUgsR0FBQSxDQUFJSSxnQkFBSixDQUFxQixDQUFyQixFQUF3QnZJLElBQUEsQ0FBS1UsQ0FBTCxDQUF4QixDQURnQyxDQUFBO0FBQUEsS0FuQkg7QUFBQSxJQXVCL0IsSUFBSW1CLE1BQUEsR0FBUzJHLE9BQUEsQ0FBUTNHLE1BQXJCLENBdkIrQjtBQUFBLElBd0IvQixLQUFLbkIsQ0FBQSxHQUFJLENBQVQsRUFBWUEsQ0FBQSxHQUFJbUIsTUFBQSxDQUFPbEIsTUFBdkIsRUFBK0JELENBQUEsRUFBL0IsRUFBb0M7QUFBQSxRQUNsQ3lILEdBQUEsQ0FBSUMsWUFBSixDQUFpQixDQUFqQixFQUFvQlEsVUFBcEIsRUFBZ0MvRyxNQUFBLENBQU9uQixDQUFQLENBQWhDLENBRGtDLENBQUE7QUFBQSxLQXhCTDtBQUFBLENBL0NqQztBQTRFQSxTQUFTaUksWUFBVCxDQUF1QkgsT0FBdkIsRUFBZ0NMLEdBQWhDLEVBQXFDO0FBQUEsSUFDbkMsSUFBSXpDLE9BQUEsR0FBVThDLE9BQUEsQ0FBUTlDLE9BQXRCLENBRG1DO0FBQUEsSUFHbkMsSUFBSUEsT0FBQSxDQUFROUUsRUFBUixLQUFlcEIsU0FBbkIsRUFBOEI7QUFBQSxRQUM1QjJJLEdBQUEsQ0FBSUcsZ0JBQUosQ0FBcUIsQ0FBckIsRUFBd0I1QyxPQUFBLENBQVE5RSxFQUFoQyxDQUQ0QixDQUFBO0FBQUEsS0FISztBQUFBLElBT25DdUgsR0FBQSxDQUFJQyxZQUFKLENBQWlCLENBQWpCLEVBQW9CUyxlQUFwQixFQUFxQ0wsT0FBckMsQ0FQbUMsQ0FBQTtBQUFBLElBUW5DTCxHQUFBLENBQUlHLGdCQUFKLENBQXFCLENBQXJCLEVBQXdCNUMsT0FBQSxDQUFRdkIsSUFBaEMsQ0FSbUMsQ0FBQTtBQUFBLElBU25DZ0UsR0FBQSxDQUFJQyxZQUFKLENBQWlCLENBQWpCLEVBQW9CVSxhQUFwQixFQUFtQ3BELE9BQW5DLENBVG1DLENBQUE7QUFBQSxDQTVFckM7QUF3RkEsU0FBU21ELGVBQVQsQ0FBMEJMLE9BQTFCLEVBQW1DTCxHQUFuQyxFQUF3QztBQUFBLElBQ3RDLElBQUl6QyxPQUFBLEdBQVU4QyxPQUFBLENBQVE5QyxPQUF0QixDQURzQztBQUFBLElBRXRDLElBQUkxRixJQUFBLEdBQU93SSxPQUFBLENBQVF4SSxJQUFuQixDQUZzQztBQUFBLElBR3RDLElBQUk2QixNQUFBLEdBQVMyRyxPQUFBLENBQVEzRyxNQUFyQixDQUhzQztBQUFBLElBSXRDLElBQUk0RyxRQUFBLEdBQVdELE9BQUEsQ0FBUUMsUUFBdkIsQ0FKc0M7QUFBQSxJQUt0QyxJQUFJQyxVQUFBLEdBQWFGLE9BQUEsQ0FBUUUsVUFBekIsQ0FMc0M7QUFBQSxJQU90QyxTQUFTNUksR0FBVCxJQUFnQjRGLE9BQUEsQ0FBUUksVUFBeEIsRUFBb0M7QUFBQSxRQUNsQyxJQUFJaUQsS0FBQSxHQUFRckQsT0FBQSxDQUFRSSxVQUFSLENBQW1CaEcsR0FBbkIsQ0FBWixDQURrQztBQUFBLFFBR2xDLElBQUlrSixRQUFBLEdBQVdQLFFBQUEsQ0FBUzNJLEdBQVQsQ0FBZixDQUhrQztBQUFBLFFBSWxDLElBQUlpSixLQUFBLEtBQVUsSUFBZDtBQUFBLFlBQW9CLFNBSmM7QUFBQTtBQU1sQyxRQUFBLElBQUksT0FBT0MsUUFBUCxLQUFvQixXQUF4QixFQUFxQztBQUFBLFlBQ25DaEosSUFBQSxDQUFLYyxJQUFMLENBQVVoQixHQUFWLENBRG1DLENBQUE7QUFBQSxZQUVuQ2tKLFFBQUEsR0FBV2hKLElBQUEsQ0FBS1csTUFBTCxHQUFjLENBQXpCLENBRm1DO0FBQUEsWUFHbkM4SCxRQUFBLENBQVMzSSxHQUFULENBQUEsR0FBZ0JrSixRQUFoQixDQUhtQztBQUFBLFNBTkg7QUFBQSxRQVdsQ2IsR0FBQSxDQUFJYyxXQUFKLENBQWdCRCxRQUFoQixDQVhrQyxDQUFBO0FBQUEsUUFhbEMsSUFBSTdFLElBQUEsR0FBTyxPQUFPNEUsS0FBbEIsQ0Fia0M7QUFBQSxRQWNsQyxJQUFJNUUsSUFBQSxLQUFTLFFBQVQsSUFBcUJBLElBQUEsS0FBUyxTQUE5QixJQUEyQ0EsSUFBQSxLQUFTLFFBQXhELEVBQWtFO0FBQUEsWUFDaEU0RSxLQUFBLEdBQVF0SixJQUFBLENBQUtILFNBQUwsQ0FBZXlKLEtBQWYsQ0FBUixDQURnRTtBQUFBLFNBZGhDO0FBQUEsUUFpQmxDLElBQUlHLFFBQUEsR0FBVy9FLElBQUEsR0FBTyxHQUFQLEdBQWE0RSxLQUE1QixDQWpCa0M7QUFBQSxRQWtCbEMsSUFBSUksVUFBQSxHQUFhVCxVQUFBLENBQVdRLFFBQVgsQ0FBakIsQ0FsQmtDO0FBQUEsUUFtQmxDLElBQUksT0FBT0MsVUFBUCxLQUFzQixXQUExQixFQUF1QztBQUFBLFlBQ3JDdEgsTUFBQSxDQUFPZixJQUFQLENBQVlpSSxLQUFaLENBRHFDLENBQUE7QUFBQSxZQUVyQ0ksVUFBQSxHQUFhdEgsTUFBQSxDQUFPbEIsTUFBUCxHQUFnQixDQUE3QixDQUZxQztBQUFBLFlBR3JDK0gsVUFBQSxDQUFXUSxRQUFYLENBQUEsR0FBdUJDLFVBQXZCLENBSHFDO0FBQUEsU0FuQkw7QUFBQSxRQXdCbENoQixHQUFBLENBQUljLFdBQUosQ0FBZ0JFLFVBQWhCLENBeEJrQyxDQUFBO0FBQUEsS0FQRTtBQUFBLENBeEZ4QztBQTJIQSxTQUFTQyxPQUFULENBQWtCQyxHQUFsQixFQUF1QjFJLE1BQXZCLEVBQStCO0FBQUEsSUFDN0IsT0FBUSxDQUFBQSxNQUFBLElBQVUsQ0FBVixLQUFnQjBJLEdBQUEsR0FBTSxDQUFOLENBQXhCLENBRDZCO0FBQUEsQ0EzSC9CO0FBK0hBLFNBQVNDLE1BQVQsQ0FBaUJDLEdBQWpCLEVBQXNCO0FBQUEsSUFDcEIsT0FBUUEsR0FBQSxJQUFPLENBQVIsR0FBY0EsR0FBQSxJQUFPLEVBQTVCLENBRG9CO0FBQUEsQ0EvSHRCO0FBbUlBLFNBQVNULGFBQVQsQ0FBd0JwRCxPQUF4QixFQUFpQ3lDLEdBQWpDLEVBQXNDO0FBQUEsSUFDcEMsSUFBSTdELFFBQUEsR0FBV29CLE9BQUEsQ0FBUVEsWUFBUixFQUFmLENBRG9DO0FBQUEsSUFFcEMsSUFBSS9CLElBQUEsR0FBT3VCLE9BQUEsQ0FBUXZCLElBQW5CLENBRm9DO0FBQUEsSUFHcEMsSUFBSW1DLENBQUEsR0FBSSxDQUFSLENBSG9DO0FBQUEsSUFJcEMsSUFBSUMsQ0FBQSxHQUFJLENBQVIsQ0FKb0M7QUFBQSxJQUtwQyxJQUFJOUIsS0FBQSxHQUFRSCxRQUFBLENBQVMzRCxNQUFyQixDQUxvQztBQUFBLElBTXBDLEtBQUssSUFBSTZJLENBQUEsR0FBSSxDQUFSLEVBQVdBLENBQUEsR0FBSS9FLEtBQXBCLEVBQTJCK0UsQ0FBQSxFQUEzQixFQUFnQztBQUFBLFFBQzlCLElBQUk3RSxJQUFBLEdBQU9MLFFBQUEsQ0FBU2tGLENBQVQsQ0FBWCxDQUQ4QjtBQUFBLFFBRTlCLElBQUlDLEtBQUEsR0FBUSxDQUFaLENBRjhCO0FBQUEsUUFHOUIsSUFBSXRGLElBQUEsS0FBUyxDQUFiLEVBQWdCO0FBQUEsWUFDZHNGLEtBQUEsR0FBUTlFLElBQUEsQ0FBS2hFLE1BQWIsQ0FEYztBQUFBLFNBSGM7QUFBQSxRQU05QndILEdBQUEsQ0FBSWMsV0FBSixDQUFnQkcsT0FBQSxDQUFRLENBQVIsRUFBV0ssS0FBWCxDQUFoQixDQU44QixDQUFBO0FBQUE7QUFROUI7QUFBQSxRQUFBLElBQUlDLFNBQUEsR0FBWXZGLElBQUEsS0FBUyxDQUFULEdBQWFRLElBQUEsQ0FBS2hFLE1BQUwsR0FBYyxDQUEzQixHQUErQmdFLElBQUEsQ0FBS2hFLE1BQXBELENBUjhCO0FBQUEsUUFTOUIsS0FBSyxJQUFJRCxDQUFBLEdBQUksQ0FBUixFQUFXQSxDQUFBLEdBQUlnSixTQUFwQixFQUErQmhKLENBQUEsRUFBL0IsRUFBb0M7QUFBQSxZQUNsQyxJQUFJQSxDQUFBLEtBQU0sQ0FBTixJQUFXeUQsSUFBQSxLQUFTLENBQXhCLEVBQTJCO0FBQUEsZ0JBQ3pCZ0UsR0FBQSxDQUFJYyxXQUFKLENBQWdCRyxPQUFBLENBQVEsQ0FBUixFQUFXTSxTQUFBLEdBQVksQ0FBdkIsQ0FBaEIsQ0FBQTtBQUFBLENBRHlCO0FBQUEsYUFETztBQUFBLFlBSWxDLElBQUlDLEVBQUEsR0FBS2hGLElBQUEsQ0FBS2pFLENBQUwsQ0FBUTRGLENBQUFBLENBQVIsR0FBWUEsQ0FBckIsQ0FKa0M7QUFBQSxZQUtsQyxJQUFJc0QsRUFBQSxHQUFLakYsSUFBQSxDQUFLakUsQ0FBTCxDQUFRNkYsQ0FBQUEsQ0FBUixHQUFZQSxDQUFyQixDQUxrQztBQUFBLFlBTWxDNEIsR0FBQSxDQUFJYyxXQUFKLENBQWdCSyxNQUFBLENBQU9LLEVBQVAsQ0FBaEIsQ0FOa0MsQ0FBQTtBQUFBLFlBT2xDeEIsR0FBQSxDQUFJYyxXQUFKLENBQWdCSyxNQUFBLENBQU9NLEVBQVAsQ0FBaEIsQ0FQa0MsQ0FBQTtBQUFBLFlBUWxDdEQsQ0FBQSxJQUFLcUQsRUFBTCxDQVJrQztBQUFBLFlBU2xDcEQsQ0FBQSxJQUFLcUQsRUFBTCxDQVRrQztBQUFBLFNBVE47QUFBQSxRQW9COUIsSUFBSXpGLElBQUEsS0FBUyxDQUFiLEVBQWdCO0FBQUEsWUFDZGdFLEdBQUEsQ0FBSWMsV0FBSixDQUFnQkcsT0FBQSxDQUFRLENBQVIsRUFBVyxDQUFYLENBQWhCLENBQUE7QUFBQSxDQURjO0FBQUEsU0FwQmM7QUFBQSxLQU5JO0FBQUEsQ0FuSXRDO0FBbUtBLFNBQVNSLFVBQVQsQ0FBcUJHLEtBQXJCLEVBQTRCWixHQUE1QixFQUFpQztBQUFBLElBQy9CLElBQUloRSxJQUFBLEdBQU8sT0FBTzRFLEtBQWxCLENBRCtCO0FBQUEsSUFFL0IsSUFBSTVFLElBQUEsS0FBUyxRQUFiLEVBQXVCO0FBQUEsUUFDckJnRSxHQUFBLENBQUlJLGdCQUFKLENBQXFCLENBQXJCLEVBQXdCUSxLQUF4QixDQURxQixDQUFBO0FBQUEsS0FBdkIsTUFFTyxJQUFJNUUsSUFBQSxLQUFTLFNBQWIsRUFBd0I7QUFBQSxRQUM3QmdFLEdBQUEsQ0FBSTBCLGlCQUFKLENBQXNCLENBQXRCLEVBQXlCZCxLQUF6QixDQUQ2QixDQUFBO0FBQUEsS0FBeEIsTUFFQSxJQUFJNUUsSUFBQSxLQUFTLFFBQWIsRUFBdUI7QUFBQSxRQUM1QixJQUFJNEUsS0FBQSxHQUFRLENBQVIsS0FBYyxDQUFsQixFQUFxQjtBQUFBLFlBQ25CWixHQUFBLENBQUkyQixnQkFBSixDQUFxQixDQUFyQixFQUF3QmYsS0FBeEIsQ0FEbUIsQ0FBQTtBQUFBLFNBQXJCLE1BRU8sSUFBSUEsS0FBQSxHQUFRLENBQVosRUFBZTtBQUFBLFlBQ3BCWixHQUFBLENBQUk0QixpQkFBSixDQUFzQixDQUF0QixFQUF5QmhCLEtBQXpCLENBRG9CLENBQUE7QUFBQSxTQUFmLE1BRUE7QUFBQSxZQUNMWixHQUFBLENBQUlHLGdCQUFKLENBQXFCLENBQXJCLEVBQXdCUyxLQUF4QixDQURLLENBQUE7QUFBQSxTQUxxQjtBQUFBLEtBTkM7QUFBQTs7QUNsS2xCLFNBQVNpQixNQUFULENBQWdCQyxHQUFoQixFQUFxQkMsTUFBckIsRUFBNkJDLFFBQTdCLEVBQXVDQyxJQUF2QyxFQUE2Q0MsS0FBN0MsRUFBb0RDLEtBQXBELEVBQTJEO0FBQUEsSUFDdEUsSUFBSUQsS0FBQSxHQUFRRCxJQUFSLElBQWdCRCxRQUFwQjtBQUFBLFFBQThCLE9BRHdDO0FBQUEsSUFHdEUsTUFBTWxGLENBQUEsR0FBS21GLElBQUEsR0FBT0MsS0FBUixJQUFrQixDQUE1QixDQUhzRTtBQUFBLElBS3RFRSxNQUFBLENBQU9OLEdBQVAsRUFBWUMsTUFBWixFQUFvQmpGLENBQXBCLEVBQXVCbUYsSUFBdkIsRUFBNkJDLEtBQTdCLEVBQW9DQyxLQUFBLEdBQVEsQ0FBNUMsQ0FMc0UsQ0FBQTtBQUFBLElBT3RFTixNQUFBLENBQU9DLEdBQVAsRUFBWUMsTUFBWixFQUFvQkMsUUFBcEIsRUFBOEJDLElBQTlCLEVBQW9DbkYsQ0FBQSxHQUFJLENBQXhDLEVBQTJDcUYsS0FBQSxHQUFRLENBQW5ELENBUHNFLENBQUE7QUFBQSxJQVF0RU4sTUFBQSxDQUFPQyxHQUFQLEVBQVlDLE1BQVosRUFBb0JDLFFBQXBCLEVBQThCbEYsQ0FBQSxHQUFJLENBQWxDLEVBQXFDb0YsS0FBckMsRUFBNENDLEtBQUEsR0FBUSxDQUFwRCxDQVJzRSxDQUFBO0FBQUEsQ0FEMUU7QUFZQSxTQUFTQyxNQUFULENBQWdCTixHQUFoQixFQUFxQkMsTUFBckIsRUFBNkI5SixDQUE3QixFQUFnQ2dLLElBQWhDLEVBQXNDQyxLQUF0QyxFQUE2Q0csR0FBN0MsRUFBa0Q7QUFBQSxJQUU5QyxPQUFPSCxLQUFBLEdBQVFELElBQWYsRUFBcUI7QUFBQSxRQUNqQixJQUFJQyxLQUFBLEdBQVFELElBQVIsR0FBZSxHQUFuQixFQUF3QjtBQUFBLFlBQ3BCLE1BQU1LLENBQUEsR0FBSUosS0FBQSxHQUFRRCxJQUFSLEdBQWUsQ0FBekIsQ0FEb0I7QUFBQSxZQUVwQixNQUFNbkYsQ0FBQSxHQUFJN0UsQ0FBQSxHQUFJZ0ssSUFBSixHQUFXLENBQXJCLENBRm9CO0FBQUEsWUFHcEIsTUFBTTVELENBQUEsR0FBSXRCLElBQUEsQ0FBS3dGLEdBQUwsQ0FBU0QsQ0FBVCxDQUFWLENBSG9CO0FBQUEsWUFJcEIsTUFBTUUsQ0FBQSxHQUFJLEdBQUEsR0FBTXpGLElBQUEsQ0FBSzBGLEdBQUwsQ0FBUyxDQUFJcEUsR0FBQUEsQ0FBSixHQUFRLENBQWpCLENBQWhCLENBSm9CO0FBQUEsWUFLcEIsTUFBTXFFLEVBQUEsR0FBSyxHQUFBLEdBQU0zRixJQUFBLENBQUs0RixJQUFMLENBQVV0RSxDQUFBLEdBQUltRSxDQUFKLElBQVNGLENBQUEsR0FBSUUsQ0FBSixDQUFULEdBQWtCRixDQUE1QixDQUFOLElBQXdDeEYsQ0FBQSxHQUFJd0YsQ0FBQSxHQUFJLENBQVIsR0FBWSxDQUFaLEdBQWdCLENBQUMsQ0FBakIsR0FBcUIsQ0FBckIsQ0FBbkQsQ0FMb0I7QUFBQSxZQU1wQixNQUFNTSxPQUFBLEdBQVU3RixJQUFBLENBQUtzQyxHQUFMLENBQVM0QyxJQUFULEVBQWVsRixJQUFBLENBQUs4RixLQUFMLENBQVc1SyxDQUFBLEdBQUk2RSxDQUFBLEdBQUkwRixDQUFKLEdBQVFGLENBQVosR0FBZ0JJLEVBQTNCLENBQWYsQ0FBaEIsQ0FOb0I7QUFBQSxZQU9wQixNQUFNSSxRQUFBLEdBQVcvRixJQUFBLENBQUtxQyxHQUFMLENBQVM4QyxLQUFULEVBQWdCbkYsSUFBQSxDQUFLOEYsS0FBTCxDQUFXNUssQ0FBQSxHQUFLLENBQUFxSyxDQUFBLEdBQUl4RixDQUFKLElBQVMwRixDQUFWLEdBQWNGLENBQWxCLEdBQXNCSSxFQUFqQyxDQUFoQixDQUFqQixDQVBvQjtBQUFBLFlBUXBCTixNQUFBLENBQU9OLEdBQVAsRUFBWUMsTUFBWixFQUFvQjlKLENBQXBCLEVBQXVCMkssT0FBdkIsRUFBZ0NFLFFBQWhDLEVBQTBDVCxHQUExQyxDQVJvQixDQUFBO0FBQUEsU0FEUDtBQUFBLFFBWWpCLE1BQU1VLENBQUEsR0FBSWhCLE1BQUEsQ0FBTyxJQUFJOUosQ0FBSixHQUFRb0ssR0FBZixDQUFWLENBWmlCO0FBQUEsUUFhakIsSUFBSTlKLENBQUEsR0FBSTBKLElBQVIsQ0FiaUI7QUFBQSxRQWNqQixJQUFJcEYsQ0FBQSxHQUFJcUYsS0FBUixDQWRpQjtBQUFBLFFBZ0JqQmMsUUFBQSxDQUFTbEIsR0FBVCxFQUFjQyxNQUFkLEVBQXNCRSxJQUF0QixFQUE0QmhLLENBQTVCLENBaEJpQixDQUFBO0FBQUEsUUFpQmpCLElBQUk4SixNQUFBLENBQU8sQ0FBQSxHQUFJRyxLQUFKLEdBQVlHLEdBQW5CLElBQTBCVSxDQUE5QjtBQUFBLFlBQWlDQyxRQUFBLENBQVNsQixHQUFULEVBQWNDLE1BQWQsRUFBc0JFLElBQXRCLEVBQTRCQyxLQUE1QixDQWpCaEIsQ0FBQTtBQUFBLFFBbUJqQixPQUFPM0osQ0FBQSxHQUFJc0UsQ0FBWCxFQUFjO0FBQUEsWUFDVm1HLFFBQUEsQ0FBU2xCLEdBQVQsRUFBY0MsTUFBZCxFQUFzQnhKLENBQXRCLEVBQXlCc0UsQ0FBekIsQ0FEVSxDQUFBO0FBQUEsWUFFVnRFLENBQUEsRUFGVSxDQUFBO0FBQUEsWUFHVnNFLENBQUEsRUFIVSxDQUFBO0FBQUEsWUFJVixPQUFPa0YsTUFBQSxDQUFPLENBQUEsR0FBSXhKLENBQUosR0FBUThKLEdBQWYsSUFBc0JVLENBQTdCO0FBQUEsZ0JBQWdDeEssQ0FBQSxFQUp0QixDQUFBO0FBQUEsWUFLVixPQUFPd0osTUFBQSxDQUFPLENBQUEsR0FBSWxGLENBQUosR0FBUXdGLEdBQWYsSUFBc0JVLENBQTdCO0FBQUEsZ0JBQWdDbEcsQ0FBQSxFQUx0QixDQUFBO0FBQUEsU0FuQkc7QUFBQSxRQTJCakIsSUFBSWtGLE1BQUEsQ0FBTyxDQUFBLEdBQUlFLElBQUosR0FBV0ksR0FBbEIsTUFBMkJVLENBQS9CO0FBQUEsWUFBa0NDLFFBQUEsQ0FBU2xCLEdBQVQsRUFBY0MsTUFBZCxFQUFzQkUsSUFBdEIsRUFBNEJwRixDQUE1QixDQUFsQyxDQUFBO0FBQUEsYUFDSztBQUFBLFlBQ0RBLENBQUEsRUFEQyxDQUFBO0FBQUEsWUFFRG1HLFFBQUEsQ0FBU2xCLEdBQVQsRUFBY0MsTUFBZCxFQUFzQmxGLENBQXRCLEVBQXlCcUYsS0FBekIsQ0FGQyxDQUFBO0FBQUEsU0E1Qlk7QUFBQSxRQWlDakIsSUFBSXJGLENBQUEsSUFBSzVFLENBQVQ7QUFBQSxZQUFZZ0ssSUFBQSxHQUFPcEYsQ0FBQSxHQUFJLENBQVgsQ0FqQ0s7QUFBQSxRQWtDakIsSUFBSTVFLENBQUEsSUFBSzRFLENBQVQ7QUFBQSxZQUFZcUYsS0FBQSxHQUFRckYsQ0FBQSxHQUFJLENBQVosQ0FsQ0s7QUFBQSxLQUZ5QjtBQUFBLENBWmxEO0FBb0RBLFNBQVNtRyxRQUFULENBQWtCbEIsR0FBbEIsRUFBdUJDLE1BQXZCLEVBQStCeEosQ0FBL0IsRUFBa0NzRSxDQUFsQyxFQUFxQztBQUFBLElBQ2pDb0csSUFBQSxDQUFLbkIsR0FBTCxFQUFVdkosQ0FBVixFQUFhc0UsQ0FBYixDQURpQyxDQUFBO0FBQUEsSUFFakNvRyxJQUFBLENBQUtsQixNQUFMLEVBQWEsSUFBSXhKLENBQWpCLEVBQW9CLENBQUlzRSxHQUFBQSxDQUF4QixDQUZpQyxDQUFBO0FBQUEsSUFHakNvRyxJQUFBLENBQUtsQixNQUFMLEVBQWEsQ0FBSXhKLEdBQUFBLENBQUosR0FBUSxDQUFyQixFQUF3QixDQUFBLEdBQUlzRSxDQUFKLEdBQVEsQ0FBaEMsQ0FIaUMsQ0FBQTtBQUFBLENBcERyQztBQTBEQSxTQUFTb0csSUFBVCxDQUFjQyxHQUFkLEVBQW1CM0ssQ0FBbkIsRUFBc0JzRSxDQUF0QixFQUF5QjtBQUFBLElBQ3JCLE1BQU1zRyxHQUFBLEdBQU1ELEdBQUEsQ0FBSTNLLENBQUosQ0FBWixDQURxQjtBQUFBLElBRXJCMkssR0FBQSxDQUFJM0ssQ0FBSixDQUFBLEdBQVMySyxHQUFBLENBQUlyRyxDQUFKLENBQVQsQ0FGcUI7QUFBQSxJQUdyQnFHLEdBQUEsQ0FBSXJHLENBQUosQ0FBQSxHQUFTc0csR0FBVCxDQUhxQjtBQUFBOztBQ3pEVixTQUFTQyxLQUFULENBQWV0QixHQUFmLEVBQW9CQyxNQUFwQixFQUE0QnNCLElBQTVCLEVBQWtDQyxJQUFsQyxFQUF3Q0MsSUFBeEMsRUFBOENDLElBQTlDLEVBQW9EeEIsUUFBcEQsRUFBOEQ7QUFBQSxJQUN6RSxNQUFNeUIsS0FBQSxHQUFRO0FBQUEsUUFBQyxDQUFEO0FBQUEsUUFBSTNCLEdBQUEsQ0FBSXRKLE1BQUosR0FBYSxDQUFqQjtBQUFBLFFBQW9CLENBQXBCO0FBQUEsS0FBZCxDQUR5RTtBQUFBLElBRXpFLE1BQU1JLE1BQUEsR0FBUyxFQUFmLENBRnlFO0FBQUEsSUFHekUsSUFBSXVGLENBQUosRUFBT0MsQ0FBUCxDQUh5RTtBQUFBLElBS3pFLE9BQU9xRixLQUFBLENBQU1qTCxNQUFiLEVBQXFCO0FBQUEsUUFDakIsTUFBTWtMLElBQUEsR0FBT0QsS0FBQSxDQUFNRSxHQUFOLEVBQWIsQ0FEaUI7QUFBQSxRQUVqQixNQUFNekIsS0FBQSxHQUFRdUIsS0FBQSxDQUFNRSxHQUFOLEVBQWQsQ0FGaUI7QUFBQSxRQUdqQixNQUFNMUIsSUFBQSxHQUFPd0IsS0FBQSxDQUFNRSxHQUFOLEVBQWIsQ0FIaUI7QUFBQSxRQUtqQixJQUFJekIsS0FBQSxHQUFRRCxJQUFSLElBQWdCRCxRQUFwQixFQUE4QjtBQUFBLFlBQzFCLEtBQUssSUFBSXpKLENBQUEsR0FBSTBKLElBQVIsRUFBYzFKLENBQUEsSUFBSzJKLEtBQXhCLEVBQStCM0osQ0FBQSxFQUEvQixFQUFvQztBQUFBLGdCQUNoQzRGLENBQUEsR0FBSTRELE1BQUEsQ0FBTyxDQUFJeEosR0FBQUEsQ0FBWCxDQUFKLENBRGdDO0FBQUEsZ0JBRWhDNkYsQ0FBQSxHQUFJMkQsTUFBQSxDQUFPLElBQUl4SixDQUFKLEdBQVEsQ0FBZixDQUFKLENBRmdDO0FBQUEsZ0JBR2hDLElBQUk0RixDQUFBLElBQUtrRixJQUFMLElBQWFsRixDQUFBLElBQUtvRixJQUFsQixJQUEwQm5GLENBQUEsSUFBS2tGLElBQS9CLElBQXVDbEYsQ0FBQSxJQUFLb0YsSUFBaEQ7QUFBQSxvQkFBc0Q1SyxNQUFBLENBQU9ELElBQVAsQ0FBWW1KLEdBQUEsQ0FBSXZKLENBQUosQ0FBWixDQUh0QixDQUFBO0FBQUEsYUFEVjtBQUFBLFlBTTFCLFNBTjBCO0FBQUEsU0FMYjtBQUFBLFFBY2pCLE1BQU11RSxDQUFBLEdBQUlDLElBQUEsQ0FBSzhGLEtBQUwsQ0FBWSxDQUFBWixJQUFBLEdBQU9DLEtBQVAsSUFBZ0IsQ0FBNUIsQ0FBVixDQWRpQjtBQUFBLFFBZ0JqQi9ELENBQUEsR0FBSTRELE1BQUEsQ0FBTyxDQUFJakYsR0FBQUEsQ0FBWCxDQUFKLENBaEJpQjtBQUFBLFFBaUJqQnNCLENBQUEsR0FBSTJELE1BQUEsQ0FBTyxJQUFJakYsQ0FBSixHQUFRLENBQWYsQ0FBSixDQWpCaUI7QUFBQSxRQW1CakIsSUFBSXFCLENBQUEsSUFBS2tGLElBQUwsSUFBYWxGLENBQUEsSUFBS29GLElBQWxCLElBQTBCbkYsQ0FBQSxJQUFLa0YsSUFBL0IsSUFBdUNsRixDQUFBLElBQUtvRixJQUFoRDtBQUFBLFlBQXNENUssTUFBQSxDQUFPRCxJQUFQLENBQVltSixHQUFBLENBQUloRixDQUFKLENBQVosQ0FuQnJDLENBQUE7QUFBQSxRQXFCakIsTUFBTThHLFFBQUEsR0FBWSxDQUFBRixJQUFBLEdBQU8sQ0FBUCxJQUFZLENBQTlCLENBckJpQjtBQUFBLFFBdUJqQixJQUFJQSxJQUFBLEtBQVMsQ0FBVCxHQUFhTCxJQUFBLElBQVFsRixDQUFyQixHQUF5Qm1GLElBQUEsSUFBUWxGLENBQXJDLEVBQXdDO0FBQUEsWUFDcENxRixLQUFBLENBQU05SyxJQUFOLENBQVdzSixJQUFYLENBRG9DLENBQUE7QUFBQSxZQUVwQ3dCLEtBQUEsQ0FBTTlLLElBQU4sQ0FBV21FLENBQUEsR0FBSSxDQUFmLENBRm9DLENBQUE7QUFBQSxZQUdwQzJHLEtBQUEsQ0FBTTlLLElBQU4sQ0FBV2lMLFFBQVgsQ0FIb0MsQ0FBQTtBQUFBLFNBdkJ2QjtBQUFBLFFBNEJqQixJQUFJRixJQUFBLEtBQVMsQ0FBVCxHQUFhSCxJQUFBLElBQVFwRixDQUFyQixHQUF5QnFGLElBQUEsSUFBUXBGLENBQXJDLEVBQXdDO0FBQUEsWUFDcENxRixLQUFBLENBQU05SyxJQUFOLENBQVdtRSxDQUFBLEdBQUksQ0FBZixDQURvQyxDQUFBO0FBQUEsWUFFcEMyRyxLQUFBLENBQU05SyxJQUFOLENBQVd1SixLQUFYLENBRm9DLENBQUE7QUFBQSxZQUdwQ3VCLEtBQUEsQ0FBTTlLLElBQU4sQ0FBV2lMLFFBQVgsQ0FIb0MsQ0FBQTtBQUFBLFNBNUJ2QjtBQUFBLEtBTG9EO0FBQUEsSUF3Q3pFLE9BQU9oTCxNQUFQLENBeEN5RTtBQUFBOztBQ0E5RCxTQUFTaUwsTUFBVCxDQUFnQi9CLEdBQWhCLEVBQXFCQyxNQUFyQixFQUE2QitCLEVBQTdCLEVBQWlDQyxFQUFqQyxFQUFxQzFDLENBQXJDLEVBQXdDVyxRQUF4QyxFQUFrRDtBQUFBLElBQzdELE1BQU15QixLQUFBLEdBQVE7QUFBQSxRQUFDLENBQUQ7QUFBQSxRQUFJM0IsR0FBQSxDQUFJdEosTUFBSixHQUFhLENBQWpCO0FBQUEsUUFBb0IsQ0FBcEI7QUFBQSxLQUFkLENBRDZEO0FBQUEsSUFFN0QsTUFBTUksTUFBQSxHQUFTLEVBQWYsQ0FGNkQ7QUFBQSxJQUc3RCxNQUFNb0wsRUFBQSxHQUFLM0MsQ0FBQSxHQUFJQSxDQUFmLENBSDZEO0FBQUEsSUFLN0QsT0FBT29DLEtBQUEsQ0FBTWpMLE1BQWIsRUFBcUI7QUFBQSxRQUNqQixNQUFNa0wsSUFBQSxHQUFPRCxLQUFBLENBQU1FLEdBQU4sRUFBYixDQURpQjtBQUFBLFFBRWpCLE1BQU16QixLQUFBLEdBQVF1QixLQUFBLENBQU1FLEdBQU4sRUFBZCxDQUZpQjtBQUFBLFFBR2pCLE1BQU0xQixJQUFBLEdBQU93QixLQUFBLENBQU1FLEdBQU4sRUFBYixDQUhpQjtBQUFBLFFBS2pCLElBQUl6QixLQUFBLEdBQVFELElBQVIsSUFBZ0JELFFBQXBCLEVBQThCO0FBQUEsWUFDMUIsS0FBSyxJQUFJekosQ0FBQSxHQUFJMEosSUFBUixFQUFjMUosQ0FBQSxJQUFLMkosS0FBeEIsRUFBK0IzSixDQUFBLEVBQS9CLEVBQW9DO0FBQUEsZ0JBQ2hDLElBQUkwTCxNQUFBLENBQU9sQyxNQUFBLENBQU8sQ0FBSXhKLEdBQUFBLENBQVgsQ0FBUCxFQUFzQndKLE1BQUEsQ0FBTyxDQUFJeEosR0FBQUEsQ0FBSixHQUFRLENBQWYsQ0FBdEIsRUFBeUN1TCxFQUF6QyxFQUE2Q0MsRUFBN0MsQ0FBQSxJQUFvREMsRUFBeEQ7QUFBQSxvQkFBNERwTCxNQUFBLENBQU9ELElBQVAsQ0FBWW1KLEdBQUEsQ0FBSXZKLENBQUosQ0FBWixDQUQ1QixDQUFBO0FBQUEsYUFEVjtBQUFBLFlBSTFCLFNBSjBCO0FBQUEsU0FMYjtBQUFBLFFBWWpCLE1BQU11RSxDQUFBLEdBQUlDLElBQUEsQ0FBSzhGLEtBQUwsQ0FBWSxDQUFBWixJQUFBLEdBQU9DLEtBQVAsSUFBZ0IsQ0FBNUIsQ0FBVixDQVppQjtBQUFBLFFBY2pCLE1BQU0vRCxDQUFBLEdBQUk0RCxNQUFBLENBQU8sQ0FBSWpGLEdBQUFBLENBQVgsQ0FBVixDQWRpQjtBQUFBLFFBZWpCLE1BQU1zQixDQUFBLEdBQUkyRCxNQUFBLENBQU8sSUFBSWpGLENBQUosR0FBUSxDQUFmLENBQVYsQ0FmaUI7QUFBQSxRQWlCakIsSUFBSW1ILE1BQUEsQ0FBTzlGLENBQVAsRUFBVUMsQ0FBVixFQUFhMEYsRUFBYixFQUFpQkMsRUFBakIsQ0FBQSxJQUF3QkMsRUFBNUI7QUFBQSxZQUFnQ3BMLE1BQUEsQ0FBT0QsSUFBUCxDQUFZbUosR0FBQSxDQUFJaEYsQ0FBSixDQUFaLENBakJmLENBQUE7QUFBQSxRQW1CakIsTUFBTThHLFFBQUEsR0FBWSxDQUFBRixJQUFBLEdBQU8sQ0FBUCxJQUFZLENBQTlCLENBbkJpQjtBQUFBLFFBcUJqQixJQUFJQSxJQUFBLEtBQVMsQ0FBVCxHQUFhSSxFQUFBLEdBQUt6QyxDQUFMLElBQVVsRCxDQUF2QixHQUEyQjRGLEVBQUEsR0FBSzFDLENBQUwsSUFBVWpELENBQXpDLEVBQTRDO0FBQUEsWUFDeENxRixLQUFBLENBQU05SyxJQUFOLENBQVdzSixJQUFYLENBRHdDLENBQUE7QUFBQSxZQUV4Q3dCLEtBQUEsQ0FBTTlLLElBQU4sQ0FBV21FLENBQUEsR0FBSSxDQUFmLENBRndDLENBQUE7QUFBQSxZQUd4QzJHLEtBQUEsQ0FBTTlLLElBQU4sQ0FBV2lMLFFBQVgsQ0FId0MsQ0FBQTtBQUFBLFNBckIzQjtBQUFBLFFBMEJqQixJQUFJRixJQUFBLEtBQVMsQ0FBVCxHQUFhSSxFQUFBLEdBQUt6QyxDQUFMLElBQVVsRCxDQUF2QixHQUEyQjRGLEVBQUEsR0FBSzFDLENBQUwsSUFBVWpELENBQXpDLEVBQTRDO0FBQUEsWUFDeENxRixLQUFBLENBQU05SyxJQUFOLENBQVdtRSxDQUFBLEdBQUksQ0FBZixDQUR3QyxDQUFBO0FBQUEsWUFFeEMyRyxLQUFBLENBQU05SyxJQUFOLENBQVd1SixLQUFYLENBRndDLENBQUE7QUFBQSxZQUd4Q3VCLEtBQUEsQ0FBTTlLLElBQU4sQ0FBV2lMLFFBQVgsQ0FId0MsQ0FBQTtBQUFBLFNBMUIzQjtBQUFBLEtBTHdDO0FBQUEsSUFzQzdELE9BQU9oTCxNQUFQLENBdEM2RDtBQUFBLENBRGpFO0FBMENBLFNBQVNxTCxNQUFULENBQWdCQyxFQUFoQixFQUFvQkMsRUFBcEIsRUFBd0JDLEVBQXhCLEVBQTRCQyxFQUE1QixFQUFnQztBQUFBLElBQzVCLE1BQU03QyxFQUFBLEdBQUswQyxFQUFBLEdBQUtFLEVBQWhCLENBRDRCO0FBQUEsSUFFNUIsTUFBTTNDLEVBQUEsR0FBSzBDLEVBQUEsR0FBS0UsRUFBaEIsQ0FGNEI7QUFBQSxJQUc1QixPQUFPN0MsRUFBQSxHQUFLQSxFQUFMLEdBQVVDLEVBQUEsR0FBS0EsRUFBdEIsQ0FINEI7QUFBQTs7QUNyQ2hDLE1BQU02QyxXQUFBLEdBQWNDLENBQUEsSUFBS0EsQ0FBQSxDQUFFLENBQUYsQ0FBekIsQ0FMQTtBQU1BLE1BQU1DLFdBQUEsR0FBY0QsQ0FBQSxJQUFLQSxDQUFBLENBQUUsQ0FBRixDQUF6QixDQU5BO0FBUWUsTUFBTUUsTUFBTixDQUFhO0FBQUEsSUFDeEIzTCxXQUFBLENBQVk0TCxNQUFaLEVBQW9CQyxJQUFBLEdBQU9MLFdBQTNCLEVBQXdDTSxJQUFBLEdBQU9KLFdBQS9DLEVBQTREeEMsUUFBQSxHQUFXLEVBQXZFLEVBQTJFNkMsU0FBQSxHQUFZQyxZQUF2RixFQUFxRztBQUFBLFFBQ2pHLElBQUs5QyxDQUFBQSxRQUFMLEdBQWdCQSxRQUFoQixDQURpRztBQUFBLFFBRWpHLElBQUswQyxDQUFBQSxNQUFMLEdBQWNBLE1BQWQsQ0FGaUc7QUFBQSxRQUlqRyxNQUFNSyxjQUFBLEdBQWlCTCxNQUFBLENBQU9sTSxNQUFQLEdBQWdCLEtBQWhCLEdBQXdCd00sV0FBeEIsR0FBc0NDLFdBQTdELENBSmlHO0FBQUEsUUFNakcsTUFBTW5ELEdBQUEsR0FBTSxJQUFBLENBQUtBLEdBQUwsR0FBVyxJQUFJaUQsY0FBSixDQUFtQkwsTUFBQSxDQUFPbE0sTUFBMUIsQ0FBdkIsQ0FOaUc7QUFBQSxRQU9qRyxNQUFNdUosTUFBQSxHQUFTLElBQUEsQ0FBS0EsTUFBTCxHQUFjLElBQUk4QyxTQUFKLENBQWNILE1BQUEsQ0FBT2xNLE1BQVAsR0FBZ0IsQ0FBOUIsQ0FBN0IsQ0FQaUc7QUFBQSxRQVNqRyxLQUFLLElBQUlELENBQUEsR0FBSSxDQUFSLEVBQVdBLENBQUEsR0FBSW1NLE1BQUEsQ0FBT2xNLE1BQTNCLEVBQW1DRCxDQUFBLEVBQW5DLEVBQXdDO0FBQUEsWUFDcEN1SixHQUFBLENBQUl2SixDQUFKLENBQUEsR0FBU0EsQ0FBVCxDQURvQztBQUFBLFlBRXBDd0osTUFBQSxDQUFPLENBQUl4SixHQUFBQSxDQUFYLENBQWdCb00sR0FBQUEsSUFBQSxDQUFLRCxNQUFBLENBQU9uTSxDQUFQLENBQUwsQ0FBaEIsQ0FGb0M7QUFBQSxZQUdwQ3dKLE1BQUEsQ0FBTyxDQUFJeEosR0FBQUEsQ0FBSixHQUFRLENBQWYsQ0FBQSxHQUFvQnFNLElBQUEsQ0FBS0YsTUFBQSxDQUFPbk0sQ0FBUCxDQUFMLENBQXBCLENBSG9DO0FBQUEsU0FUeUQ7QUFBQSxRQWVqR1QsTUFBQSxDQUFLZ0ssR0FBTCxFQUFVQyxNQUFWLEVBQWtCQyxRQUFsQixFQUE0QixDQUE1QixFQUErQkYsR0FBQSxDQUFJdEosTUFBSixHQUFhLENBQTVDLEVBQStDLENBQS9DLENBZmlHLENBQUE7QUFBQSxLQUQ3RTtBQUFBLElBbUJ4QjRLLEtBQUEsQ0FBTUMsSUFBTixFQUFZQyxJQUFaLEVBQWtCQyxJQUFsQixFQUF3QkMsSUFBeEIsRUFBOEI7QUFBQSxRQUMxQixPQUFPSixLQUFBLENBQU0sSUFBQSxDQUFLdEIsR0FBWCxFQUFnQixJQUFBLENBQUtDLE1BQXJCLEVBQTZCc0IsSUFBN0IsRUFBbUNDLElBQW5DLEVBQXlDQyxJQUF6QyxFQUErQ0MsSUFBL0MsRUFBcUQsSUFBS3hCLENBQUFBLFFBQTFELENBQVAsQ0FEMEI7QUFBQSxLQW5CTjtBQUFBLElBdUJ4QjZCLE1BQUEsQ0FBTzFGLENBQVAsRUFBVUMsQ0FBVixFQUFhaUQsQ0FBYixFQUFnQjtBQUFBLFFBQ1osT0FBT3dDLE1BQUEsQ0FBTyxLQUFLL0IsR0FBWixFQUFpQixLQUFLQyxNQUF0QixFQUE4QjVELENBQTlCLEVBQWlDQyxDQUFqQyxFQUFvQ2lELENBQXBDLEVBQXVDLElBQUtXLENBQUFBLFFBQTVDLENBQVAsQ0FEWTtBQUFBLEtBdkJRO0FBQUE7O0FDTDVCLE1BQU1rRCxjQUFBLEdBQWlCO0FBQUEsSUFDbkJDLE9BQUEsRUFBUyxDQURVO0FBQUE7QUFFbkIsSUFBQUMsT0FBQSxFQUFTLEVBRlU7QUFBQTtBQUduQixJQUFBQyxTQUFBLEVBQVcsQ0FIUTtBQUFBO0FBSW5CLElBQUFDLE1BQUEsRUFBUSxFQUpXO0FBQUE7QUFLbkIsSUFBQTdILE1BQUEsRUFBUSxHQUxXO0FBQUE7QUFNbkIsSUFBQXVFLFFBQUEsRUFBVSxFQU5TO0FBQUE7QUFPbkIsSUFBQU8sR0FBQSxFQUFLLEtBUGM7QUFBQTtBQVVuQjtBQUFBLElBQUFnRCxVQUFBLEVBQVksS0FWTztBQUFBO0FBYW5CLElBQUFDLE1BQUEsRUFBUSxJQWJXO0FBQUE7QUFnQm5CO0FBQUEsSUFBQTdMLEdBQUEsRUFBSzhMLEtBQUEsSUFBU0EsS0FBQUE7QUFoQkssQ0FBdkIsQ0FIQTtBQXNCQSxNQUFNQyxNQUFBLEdBQVMzSSxJQUFBLENBQUsySSxNQUFMLElBQWdCLENBQUF2QyxHQUFBLElBQVNoRixDQUFELElBQU87QUFBQSxJQUFFZ0YsR0FBQSxDQUFJLENBQUosQ0FBUyxHQUFBLENBQUNoRixDQUFWLENBQUY7QUFBQSxJQUFlLE9BQU9nRixHQUFBLENBQUksQ0FBSixDQUFQLENBQWY7QUFBQSxDQUFmLEVBQWlELElBQUl3QyxZQUFKLENBQWlCLENBQWpCLENBQWxELENBQTlCLENBdEJBO0FBd0JlLE1BQU1DLFlBQU4sQ0FBbUI7QUFBQSxJQUM5QjlNLFdBQUEsQ0FBWTZGLE9BQVosRUFBcUI7QUFBQSxRQUNqQixJQUFLQSxDQUFBQSxPQUFMLEdBQWVrSCxRQUFBLENBQU9qTyxNQUFBLENBQU9rTyxNQUFQLENBQWNaLGNBQWQsQ0FBUCxFQUFzQ3ZHLE9BQXRDLENBQWYsQ0FEaUI7QUFBQSxRQUVqQixJQUFLb0gsQ0FBQUEsS0FBTCxHQUFhLElBQUl4TyxLQUFKLENBQVUsSUFBS29ILENBQUFBLE9BQUwsQ0FBYXlHLE9BQWIsR0FBdUIsQ0FBakMsQ0FBYixDQUZpQjtBQUFBLEtBRFM7QUFBQSxJQU05QlksSUFBQSxDQUFLdEIsTUFBTCxFQUFhO0FBQUEsUUFDVCxNQUFNLENBQUNuQyxHQUFELEVBQU00QyxPQUFOLEVBQWVDLE9BQWYsRUFBd0JwRCxRQUF4QixDQUFBLEdBQW9DLElBQUtyRCxDQUFBQSxPQUEvQyxDQURTO0FBQUEsUUFHVCxJQUFJNEQsR0FBSjtBQUFBLFlBQVMwRCxPQUFBLENBQVFDLElBQVIsQ0FBYSxZQUFiLENBSEEsQ0FBQTtBQUFBLFFBS1QsTUFBTUMsT0FBQSxHQUFVLENBQUMsUUFBRCxHQUFhekIsTUFBQSxDQUFPbE0sTUFBcEIsRUFBNkIsT0FBN0IsQ0FBaEIsQ0FMUztBQUFBLFFBTVQsSUFBSStKLEdBQUo7QUFBQSxZQUFTMEQsT0FBQSxDQUFRQyxJQUFSLENBQWFDLE9BQWIsQ0FOQSxDQUFBO0FBQUEsUUFRVCxJQUFLekIsQ0FBQUEsTUFBTCxHQUFjQSxNQUFkLENBUlM7QUFBQTtBQVdULFFBQUkwQixJQUFBQSxRQUFBLEdBQVcsRUFBZixDQVhTO0FBQUEsUUFZVCxLQUFLLElBQUk3TixDQUFBLEdBQUksQ0FBUixFQUFXQSxDQUFBLEdBQUltTSxNQUFBLENBQU9sTSxNQUEzQixFQUFtQ0QsQ0FBQSxFQUFuQyxFQUF3QztBQUFBLFlBQ3BDLElBQUksQ0FBQ21NLE1BQUEsQ0FBT25NLENBQVAsRUFBVTRELFFBQWY7QUFBQSxnQkFBeUIsU0FEVztBQUFBLFlBRXBDaUssUUFBQSxDQUFTek4sSUFBVCxDQUFjME4sa0JBQUEsQ0FBbUIzQixNQUFBLENBQU9uTSxDQUFQLENBQW5CLEVBQThCQSxDQUE5QixDQUFkLENBRm9DLENBQUE7QUFBQSxTQVovQjtBQUFBLFFBZ0JULElBQUEsQ0FBS3dOLEtBQUwsQ0FBV1gsT0FBQSxHQUFVLENBQXJCLENBQUEsR0FBMEIsSUFBSVgsTUFBSixDQUFXMkIsUUFBWCxFQUFxQnpCLElBQXJCLEVBQTJCQyxJQUEzQixFQUFpQzVDLFFBQWpDLEVBQTJDMkQsWUFBM0MsQ0FBMUIsQ0FoQlM7QUFBQSxRQWtCVCxJQUFJcEQsR0FBSjtBQUFBLFlBQVMwRCxPQUFBLENBQVFLLE9BQVIsQ0FBZ0JILE9BQWhCLENBbEJBLENBQUE7QUFBQTtBQXNCVDtBQUFBLFFBQUEsS0FBSyxJQUFJOUgsQ0FBQSxHQUFJK0csT0FBUixFQUFpQi9HLENBQUEsSUFBSzhHLE9BQTNCLEVBQW9DOUcsQ0FBQSxFQUFwQyxFQUF5QztBQUFBLFlBQ3JDLE1BQU1rSSxHQUFBLEdBQU0sQ0FBQ0MsSUFBQSxDQUFLRCxHQUFMLEVBQWIsQ0FEcUM7QUFBQTtBQUlyQyxZQUFBSCxRQUFBLEdBQVcsSUFBS0ssQ0FBQUEsUUFBTCxDQUFjTCxRQUFkLEVBQXdCL0gsQ0FBeEIsQ0FBWCxDQUpxQztBQUFBLFlBS3JDLEtBQUswSCxLQUFMLENBQVcxSCxDQUFYLENBQWdCLEdBQUEsSUFBSW9HLE1BQUosQ0FBVzJCLFFBQVgsRUFBcUJ6QixJQUFyQixFQUEyQkMsSUFBM0IsRUFBaUM1QyxRQUFqQyxFQUEyQzJELFlBQTNDLENBQWhCLENBTHFDO0FBQUEsWUFPckMsSUFBSXBELEdBQUo7QUFBQSxnQkFBUzBELE9BQUEsQ0FBUTFELEdBQVIsQ0FBWSwwQkFBWixFQUF3Q2xFLENBQXhDLEVBQTJDK0gsUUFBQSxDQUFTNU4sTUFBcEQsRUFBNEQsQ0FBQ2dPLElBQUEsQ0FBS0QsR0FBTCxFQUFELEdBQWNBLEdBQTFFLENBUDRCLENBQUE7QUFBQSxTQXRCaEM7QUFBQSxRQWdDVCxJQUFJaEUsR0FBSjtBQUFBLFlBQVMwRCxPQUFBLENBQVFLLE9BQVIsQ0FBZ0IsWUFBaEIsQ0FoQ0EsQ0FBQTtBQUFBLFFBa0NULE9BQU8sSUFBUCxDQWxDUztBQUFBLEtBTmlCO0FBQUEsSUEyQzlCSSxXQUFBLENBQVk3SCxJQUFaLEVBQWtCOEgsSUFBbEIsRUFBd0I7QUFBQSxRQUNwQixJQUFJQyxNQUFBLEdBQVUsQ0FBQyxDQUFBL0gsSUFBQSxDQUFLLENBQUwsQ0FBVSxHQUFBLEdBQVYsSUFBaUIsR0FBbEIsR0FBd0IsR0FBeEIsSUFBK0IsR0FBaEMsR0FBc0MsR0FBbkQsQ0FEb0I7QUFBQSxRQUVwQixNQUFNZ0ksTUFBQSxHQUFTOUosSUFBQSxDQUFLc0MsR0FBTCxDQUFTLENBQUMsRUFBVixFQUFjdEMsSUFBQSxDQUFLcUMsR0FBTCxDQUFTLEVBQVQsRUFBYVAsSUFBQSxDQUFLLENBQUwsQ0FBYixDQUFkLENBQWYsQ0FGb0I7QUFBQSxRQUdwQixJQUFJaUksTUFBQSxHQUFTakksSUFBQSxDQUFLLENBQUwsQ0FBWSxLQUFBLEdBQVosR0FBa0IsR0FBbEIsR0FBeUIsQ0FBQyxDQUFBQSxJQUFBLENBQUssQ0FBTCxDQUFVLEdBQUEsR0FBVixJQUFpQixHQUFsQixHQUF3QixHQUF4QixJQUErQixHQUFoQyxHQUFzQyxHQUEzRSxDQUhvQjtBQUFBLFFBSXBCLE1BQU1rSSxNQUFBLEdBQVNoSyxJQUFBLENBQUtzQyxHQUFMLENBQVMsQ0FBQyxFQUFWLEVBQWN0QyxJQUFBLENBQUtxQyxHQUFMLENBQVMsRUFBVCxFQUFhUCxJQUFBLENBQUssQ0FBTCxDQUFiLENBQWQsQ0FBZixDQUpvQjtBQUFBLFFBTXBCLElBQUlBLElBQUEsQ0FBSyxDQUFMLENBQUEsR0FBVUEsSUFBQSxDQUFLLENBQUwsQ0FBVixJQUFxQixHQUF6QixFQUE4QjtBQUFBLFlBQzFCK0gsTUFBQSxHQUFTLENBQUMsR0FBVixDQUQwQjtBQUFBLFlBRTFCRSxNQUFBLEdBQVMsR0FBVCxDQUYwQjtBQUFBLFNBQTlCLE1BR08sSUFBSUYsTUFBQSxHQUFTRSxNQUFiLEVBQXFCO0FBQUEsWUFDeEIsTUFBTUUsVUFBQSxHQUFhLElBQUEsQ0FBS04sV0FBTCxDQUFpQjtBQUFBLGdCQUFDRSxNQUFEO0FBQUEsZ0JBQVNDLE1BQVQ7QUFBQSxnQkFBaUIsR0FBakI7QUFBQSxnQkFBc0JFLE1BQXRCO0FBQUEsYUFBakIsRUFBZ0RKLElBQWhELENBQW5CLENBRHdCO0FBQUEsWUFFeEIsTUFBTU0sVUFBQSxHQUFhLElBQUEsQ0FBS1AsV0FBTCxDQUFpQjtBQUFBLGdCQUFDLENBQUMsR0FBRjtBQUFBLGdCQUFPRyxNQUFQO0FBQUEsZ0JBQWVDLE1BQWY7QUFBQSxnQkFBdUJDLE1BQXZCO0FBQUEsYUFBakIsRUFBaURKLElBQWpELENBQW5CLENBRndCO0FBQUEsWUFHeEIsT0FBT0ssVUFBQSxDQUFXRSxNQUFYLENBQWtCRCxVQUFsQixDQUFQLENBSHdCO0FBQUEsU0FUUjtBQUFBLFFBZXBCLE1BQU1FLElBQUEsR0FBTyxJQUFBLENBQUtwQixLQUFMLENBQVcsSUFBS3FCLENBQUFBLFVBQUwsQ0FBZ0JULElBQWhCLENBQVgsQ0FBYixDQWZvQjtBQUFBLFFBZ0JwQixNQUFNN0UsR0FBQSxHQUFNcUYsSUFBQSxDQUFLL0QsS0FBTCxDQUFXaUUsSUFBQSxDQUFLVCxNQUFMLENBQVgsRUFBeUJVLElBQUEsQ0FBS1AsTUFBTCxDQUF6QixFQUF1Q00sSUFBQSxDQUFLUCxNQUFMLENBQXZDLEVBQXFEUSxJQUFBLENBQUtULE1BQUwsQ0FBckQsQ0FBWixDQWhCb0I7QUFBQSxRQWlCcEIsTUFBTVQsUUFBQSxHQUFXLEVBQWpCLENBakJvQjtBQUFBLFFBa0JwQixLQUFXM04sTUFBQUEsRUFBWCxJQUFpQnFKLEdBQWpCLEVBQXNCO0FBQUEsWUFDbEIsTUFBTXlGLENBQUEsR0FBSUosSUFBQSxDQUFLekMsTUFBTCxDQUFZak0sRUFBWixDQUFWLENBRGtCO0FBQUEsWUFFbEIyTixRQUFBLENBQVN6TixJQUFULENBQWM0TyxDQUFBLENBQUVDLFNBQUYsR0FBY0MsY0FBQSxDQUFlRixDQUFmLENBQWQsR0FBa0MsS0FBSzdDLE1BQUwsQ0FBWTZDLENBQUEsQ0FBRUcsS0FBZCxDQUFoRCxDQUZrQixDQUFBO0FBQUEsU0FsQkY7QUFBQSxRQXNCcEIsT0FBT3RCLFFBQVAsQ0F0Qm9CO0FBQUEsS0EzQ007QUFBQSxJQW9FOUJ1QixXQUFBLENBQVlDLFNBQVosRUFBdUI7QUFBQSxRQUNuQixNQUFNQyxRQUFBLEdBQVcsSUFBQSxDQUFLQyxZQUFMLENBQWtCRixTQUFsQixDQUFqQixDQURtQjtBQUFBLFFBRW5CLE1BQU1HLFVBQUEsR0FBYSxJQUFBLENBQUtDLGNBQUwsQ0FBb0JKLFNBQXBCLENBQW5CLENBRm1CO0FBQUEsUUFHbkIsTUFBTUssUUFBQSxHQUFXLG1DQUFqQixDQUhtQjtBQUFBLFFBS25CLE1BQU1QLEtBQUEsR0FBUSxJQUFBLENBQUszQixLQUFMLENBQVdnQyxVQUFYLENBQWQsQ0FMbUI7QUFBQSxRQU1uQixJQUFJLENBQUNMLEtBQUw7QUFBQSxZQUFZLE1BQU0sSUFBSVEsS0FBSixDQUFVRCxRQUFWLENBQU4sQ0FOTztBQUFBLFFBUW5CLE1BQU1FLE1BQUEsR0FBU1QsS0FBQSxDQUFNaEQsTUFBTixDQUFhbUQsUUFBYixDQUFmLENBUm1CO0FBQUEsUUFTbkIsSUFBSSxDQUFDTSxNQUFMO0FBQUEsWUFBYSxNQUFNLElBQUlELEtBQUosQ0FBVUQsUUFBVixDQUFOLENBVE07QUFBQSxRQVduQixNQUFNNUcsQ0FBQSxHQUFJLElBQUsxQyxDQUFBQSxPQUFMLENBQWEyRyxNQUFiLElBQTRCM0csSUFBQUEsQ0FBQUEsT0FBTCxDQUFhbEIsTUFBYixHQUFzQlYsSUFBQSxDQUFLcUwsR0FBTCxDQUFTLENBQVQsRUFBWUwsVUFBQSxHQUFhLENBQXpCLENBQXRCLENBQWpDLENBWG1CO0FBQUEsUUFZbkIsTUFBTWpHLEdBQUEsR0FBTTRGLEtBQUEsQ0FBTTdELE1BQU4sQ0FBYXNFLE1BQUEsQ0FBT2hLLENBQXBCLEVBQXVCZ0ssTUFBQSxDQUFPL0osQ0FBOUIsRUFBaUNpRCxDQUFqQyxDQUFaLENBWm1CO0FBQUEsUUFhbkIsTUFBTWdILFFBQUEsR0FBVyxFQUFqQixDQWJtQjtBQUFBLFFBY25CLEtBQVc1UCxNQUFBQSxFQUFYLElBQWlCcUosR0FBakIsRUFBc0I7QUFBQSxZQUNsQixNQUFNeUYsQ0FBQSxHQUFJRyxLQUFBLENBQU1oRCxNQUFOLENBQWFqTSxFQUFiLENBQVYsQ0FEa0I7QUFBQSxZQUVsQixJQUFJOE8sQ0FBQSxDQUFFZSxRQUFGLEtBQWVWLFNBQW5CLEVBQThCO0FBQUEsZ0JBQzFCUyxRQUFBLENBQVMxUCxJQUFULENBQWM0TyxDQUFBLENBQUVDLFNBQUYsR0FBY0MsY0FBQSxDQUFlRixDQUFmLENBQWQsR0FBa0MsS0FBSzdDLE1BQUwsQ0FBWTZDLENBQUEsQ0FBRUcsS0FBZCxDQUFoRCxDQUQwQixDQUFBO0FBQUEsYUFGWjtBQUFBLFNBZEg7QUFBQSxRQXFCbkIsSUFBSVcsUUFBQSxDQUFTN1AsTUFBVCxLQUFvQixDQUF4QjtBQUFBLFlBQTJCLE1BQU0sSUFBSTBQLEtBQUosQ0FBVUQsUUFBVixDQUFOLENBckJSO0FBQUEsUUF1Qm5CLE9BQU9JLFFBQVAsQ0F2Qm1CO0FBQUEsS0FwRU87QUFBQSxJQThGOUJFLFNBQUEsQ0FBVVgsU0FBVixFQUFxQlksS0FBckIsRUFBNEJDLE1BQTVCLEVBQW9DO0FBQUEsUUFDaENELEtBQUEsR0FBUUEsS0FBQSxJQUFTLEVBQWpCLENBRGdDO0FBQUEsUUFFaENDLE1BQUEsR0FBU0EsTUFBQSxJQUFVLENBQW5CLENBRmdDO0FBQUEsUUFJaEMsTUFBTUMsTUFBQSxHQUFTLEVBQWYsQ0FKZ0M7QUFBQSxRQUtoQyxJQUFLQyxDQUFBQSxhQUFMLENBQW1CRCxNQUFuQixFQUEyQmQsU0FBM0IsRUFBc0NZLEtBQXRDLEVBQTZDQyxNQUE3QyxFQUFxRCxDQUFyRCxDQUxnQyxDQUFBO0FBQUEsUUFPaEMsT0FBT0MsTUFBUCxDQVBnQztBQUFBLEtBOUZOO0FBQUEsSUF3RzlCRSxPQUFBLENBQVF2SyxDQUFSLEVBQVdGLENBQVgsRUFBY0MsQ0FBZCxFQUFpQjtBQUFBLFFBQ2IsTUFBTStJLElBQUEsR0FBTyxJQUFBLENBQUtwQixLQUFMLENBQVcsSUFBS3FCLENBQUFBLFVBQUwsQ0FBZ0IvSSxDQUFoQixDQUFYLENBQWIsQ0FEYTtBQUFBLFFBRWIsTUFBTXdLLEVBQUEsR0FBSzlMLElBQUEsQ0FBS3FMLEdBQUwsQ0FBUyxDQUFULEVBQVkvSixDQUFaLENBQVgsQ0FGYTtBQUFBLFFBR2IsTUFBTSxDQUFDWixNQUFELEVBQVM2SCxNQUFULENBQW1CLEdBQUEsSUFBQSxDQUFLM0csT0FBOUIsQ0FIYTtBQUFBLFFBSWIsTUFBTTRGLENBQUEsR0FBSWUsTUFBQSxHQUFTN0gsTUFBbkIsQ0FKYTtBQUFBLFFBS2IsTUFBTXFMLEdBQUEsR0FBTyxDQUFBMUssQ0FBQSxHQUFJbUcsQ0FBSixJQUFTc0UsRUFBdEIsQ0FMYTtBQUFBLFFBTWIsTUFBTUUsTUFBQSxHQUFVLENBQUEzSyxDQUFBLEdBQUksQ0FBSixHQUFRbUcsQ0FBUixJQUFhc0UsRUFBN0IsQ0FOYTtBQUFBLFFBUWIsTUFBTW5KLElBQUEsR0FBTyxFQUNUekQsUUFBQSxFQUFVLEVBREQsRUFBYixDQVJhO0FBQUEsUUFZYixLQUFLK00sZ0JBQUwsQ0FDSTdCLElBQUEsQ0FBSy9ELEtBQUwsQ0FBWSxDQUFBakYsQ0FBQSxHQUFJb0csQ0FBSixJQUFTc0UsRUFBckIsRUFBeUJDLEdBQXpCLEVBQStCLENBQUEzSyxDQUFBLEdBQUksQ0FBSixHQUFRb0csQ0FBUixJQUFhc0UsRUFBNUMsRUFBZ0RFLE1BQWhELENBREosRUFFSTVCLElBQUEsQ0FBS3pDLE1BRlQsRUFFaUJ2RyxDQUZqQixFQUVvQkMsQ0FGcEIsRUFFdUJ5SyxFQUZ2QixFQUUyQm5KLElBRjNCLENBWmEsQ0FBQTtBQUFBLFFBZ0JiLElBQUl2QixDQUFBLEtBQU0sQ0FBVixFQUFhO0FBQUEsWUFDVCxJQUFLNkssQ0FBQUEsZ0JBQUwsQ0FDSTdCLElBQUEsQ0FBSy9ELEtBQUwsQ0FBVyxDQUFBLEdBQUltQixDQUFBLEdBQUlzRSxFQUFuQixFQUF1QkMsR0FBdkIsRUFBNEIsQ0FBNUIsRUFBK0JDLE1BQS9CLENBREosRUFFSTVCLElBQUEsQ0FBS3pDLE1BRlQsRUFFaUJtRSxFQUZqQixFQUVxQnpLLENBRnJCLEVBRXdCeUssRUFGeEIsRUFFNEJuSixJQUY1QixDQURTLENBQUE7QUFBQSxTQWhCQTtBQUFBLFFBcUJiLElBQUl2QixDQUFBLEtBQU0wSyxFQUFBLEdBQUssQ0FBZixFQUFrQjtBQUFBLFlBQ2QsSUFBS0csQ0FBQUEsZ0JBQUwsQ0FDSTdCLElBQUEsQ0FBSy9ELEtBQUwsQ0FBVyxDQUFYLEVBQWMwRixHQUFkLEVBQW1CdkUsQ0FBQSxHQUFJc0UsRUFBdkIsRUFBMkJFLE1BQTNCLENBREosRUFFSTVCLElBQUEsQ0FBS3pDLE1BRlQsRUFFaUIsQ0FBQyxDQUZsQixFQUVxQnRHLENBRnJCLEVBRXdCeUssRUFGeEIsRUFFNEJuSixJQUY1QixDQURjLENBQUE7QUFBQSxTQXJCTDtBQUFBLFFBMkJiLE9BQU9BLElBQUEsQ0FBS3pELFFBQUwsQ0FBY3pELE1BQWQsR0FBdUJrSCxJQUF2QixHQUE4QixJQUFyQyxDQTNCYTtBQUFBLEtBeEdhO0FBQUEsSUFzSTlCdUosdUJBQUEsQ0FBd0JyQixTQUF4QixFQUFtQztBQUFBLFFBQy9CLElBQUlzQixhQUFBLEdBQWdCLElBQUEsQ0FBS2xCLGNBQUwsQ0FBb0JKLFNBQXBCLENBQWlDLEdBQUEsQ0FBckQsQ0FEK0I7QUFBQSxRQUUvQixPQUFPc0IsYUFBQSxJQUFpQixLQUFLdkssT0FBTCxDQUFheUcsT0FBckMsRUFBOEM7QUFBQSxZQUMxQyxNQUFNaUQsUUFBQSxHQUFXLElBQUEsQ0FBS1YsV0FBTCxDQUFpQkMsU0FBakIsQ0FBakIsQ0FEMEM7QUFBQSxZQUUxQ3NCLGFBQUEsRUFGMEMsQ0FBQTtBQUFBLFlBRzFDLElBQUliLFFBQUEsQ0FBUzdQLE1BQVQsS0FBb0IsQ0FBeEI7QUFBQSxnQkFBMkIsTUFIZTtBQUFBLFlBSTFDb1AsU0FBQSxHQUFZUyxRQUFBLENBQVMsQ0FBVCxDQUFZMUssQ0FBQUEsVUFBWixDQUF1QndMLFVBQW5DLENBSjBDO0FBQUEsU0FGZjtBQUFBLFFBUS9CLE9BQU9ELGFBQVAsQ0FSK0I7QUFBQSxLQXRJTDtBQUFBLElBaUo5QlAsYUFBQSxDQUFjL1AsTUFBZCxFQUFzQmdQLFNBQXRCLEVBQWlDWSxLQUFqQyxFQUF3Q0MsTUFBeEMsRUFBZ0RXLE9BQWhELEVBQXlEO0FBQUEsUUFDckQsTUFBTWYsUUFBQSxHQUFXLElBQUEsQ0FBS1YsV0FBTCxDQUFpQkMsU0FBakIsQ0FBakIsQ0FEcUQ7QUFBQSxRQUdyRCxLQUFXeUIsTUFBQUEsS0FBWCxJQUFvQmhCLFFBQXBCLEVBQThCO0FBQUEsWUFDMUIsTUFBTTVDLEtBQUEsR0FBUTRELEtBQUEsQ0FBTTFMLFVBQXBCLENBRDBCO0FBQUEsWUFHMUIsSUFBSThILEtBQUEsSUFBU0EsS0FBQSxDQUFNNkQsT0FBbkIsRUFBNEI7QUFBQSxnQkFDeEIsSUFBSUYsT0FBQSxHQUFVM0QsS0FBQSxDQUFNOEQsV0FBaEIsSUFBK0JkLE1BQW5DLEVBQTJDO0FBQUE7QUFFdkMsb0JBQUFXLE9BQUEsSUFBVzNELEtBQUEsQ0FBTThELFdBQWpCLENBRnVDO0FBQUEsaUJBQTNDLE1BR087QUFBQTtBQUVILG9CQUFBSCxPQUFBLEdBQVUsSUFBQSxDQUFLVCxhQUFMLENBQW1CL1AsTUFBbkIsRUFBMkI2TSxLQUFBLENBQU0wRCxVQUFqQyxFQUE2Q1gsS0FBN0MsRUFBb0RDLE1BQXBELEVBQTREVyxPQUE1RCxDQUFWLENBQUE7QUFGRyxpQkFKaUI7QUFBQSxhQUE1QixNQVNPLElBQUlBLE9BQUEsR0FBVVgsTUFBZCxFQUFzQjtBQUFBO0FBRXpCLGdCQUFBVyxPQUFBLEVBRnlCLENBQUE7QUFBQSxhQUF0QixNQUdBO0FBQUE7QUFFSCxnQkFBQXhRLE1BQUEsQ0FBT0QsSUFBUCxDQUFZMFEsS0FBWixDQUZHLENBQUE7QUFBQSxhQWZtQjtBQUFBLFlBbUIxQixJQUFJelEsTUFBQSxDQUFPSixNQUFQLEtBQWtCZ1EsS0FBdEI7QUFBQSxnQkFBNkIsTUFuQkg7QUFBQSxTQUh1QjtBQUFBLFFBeUJyRCxPQUFPWSxPQUFQLENBekJxRDtBQUFBLEtBakozQjtBQUFBLElBNks5QkosZ0JBQUEsQ0FBaUJsSCxHQUFqQixFQUFzQjRDLE1BQXRCLEVBQThCdkcsQ0FBOUIsRUFBaUNDLENBQWpDLEVBQW9DeUssRUFBcEMsRUFBd0NuSixJQUF4QyxFQUE4QztBQUFBLFFBQzFDLEtBQVduSCxNQUFBQSxDQUFYLElBQWdCdUosR0FBaEIsRUFBcUI7QUFBQSxZQUNqQixNQUFNeUYsQ0FBQSxHQUFJN0MsTUFBQSxDQUFPbk0sQ0FBUCxDQUFWLENBRGlCO0FBQUEsWUFFakIsTUFBTWlSLFNBQUEsR0FBWWpDLENBQUEsQ0FBRUMsU0FBcEIsQ0FGaUI7QUFBQSxZQUlqQixJQUFJNUosSUFBSixFQUFVNkwsRUFBVixFQUFjQyxFQUFkLENBSmlCO0FBQUEsWUFLakIsSUFBSUYsU0FBSixFQUFlO0FBQUEsZ0JBQ1g1TCxJQUFBLEdBQU8rTCxvQkFBQSxDQUFxQnBDLENBQXJCLENBQVAsQ0FEVztBQUFBLGdCQUVYa0MsRUFBQSxHQUFLbEMsQ0FBQSxDQUFFcEosQ0FBUCxDQUZXO0FBQUEsZ0JBR1h1TCxFQUFBLEdBQUtuQyxDQUFBLENBQUVuSixDQUFQLENBSFc7QUFBQSxhQUFmLE1BSU87QUFBQSxnQkFDSCxNQUFNbUcsQ0FBQSxHQUFJLElBQUtHLENBQUFBLE1BQUwsQ0FBWTZDLENBQUEsQ0FBRUcsS0FBZCxDQUFWLENBREc7QUFBQSxnQkFFSDlKLElBQUEsR0FBTzJHLENBQUEsQ0FBRTVHLFVBQVQsQ0FGRztBQUFBLGdCQUdIOEwsRUFBQSxHQUFLcEMsSUFBQSxDQUFLOUMsQ0FBQSxDQUFFcEksUUFBRixDQUFXRSxXQUFYLENBQXVCLENBQXZCLENBQUwsQ0FBTCxDQUhHO0FBQUEsZ0JBSUhxTixFQUFBLEdBQUtwQyxJQUFBLENBQUsvQyxDQUFBLENBQUVwSSxRQUFGLENBQVdFLFdBQVgsQ0FBdUIsQ0FBdkIsQ0FBTCxDQUFMLENBSkc7QUFBQSxhQVRVO0FBQUEsWUFnQmpCLE1BQU11TixDQUFBLEdBQUk7QUFBQSxnQkFDTjVOLElBQUEsRUFBTSxDQURBO0FBQUEsZ0JBRU5HLFFBQUEsRUFBVSxDQUFDO0FBQUEsd0JBQ1BZLElBQUEsQ0FBSzhNLEtBQUwsQ0FBVyxJQUFBLENBQUtsTCxPQUFMLENBQWFsQixNQUFiLElBQXVCZ00sRUFBQSxHQUFLWixFQUFMLEdBQVUxSyxDQUFWLENBQWxDLENBRE87QUFBQSx3QkFFUHBCLElBQUEsQ0FBSzhNLEtBQUwsQ0FBVyxJQUFBLENBQUtsTCxPQUFMLENBQWFsQixNQUFiLElBQXVCaU0sRUFBQSxHQUFLYixFQUFMLEdBQVV6SyxDQUFWLENBQWxDLENBRk87QUFBQSxxQkFBRCxDQUZKO0FBQUEsZ0JBTU5SLElBTk07QUFBQSxhQUFWLENBaEJpQjtBQUFBO0FBMEJqQixZQUFBLElBQUluRixFQUFKLENBMUJpQjtBQUFBLFlBMkJqQixJQUFJK1EsU0FBSixFQUFlO0FBQUEsZ0JBQ1gvUSxFQUFBLEdBQUs4TyxDQUFBLENBQUU5TyxFQUFQLENBRFc7QUFBQSxhQUFmLE1BRU8sSUFBSSxJQUFBLENBQUtrRyxPQUFMLENBQWE0RyxVQUFqQixFQUE2QjtBQUFBO0FBRWhDLGdCQUFBOU0sRUFBQSxHQUFLOE8sQ0FBQSxDQUFFRyxLQUFQLENBRmdDO0FBQUEsYUFBN0IsTUFHQSxJQUFJLElBQUEsQ0FBS2hELE1BQUwsQ0FBWTZDLENBQUEsQ0FBRUcsS0FBZCxDQUFxQmpQLENBQUFBLEVBQXpCLEVBQTZCO0FBQUE7QUFFaEMsZ0JBQUFBLEVBQUEsR0FBSyxJQUFBLENBQUtpTSxNQUFMLENBQVk2QyxDQUFBLENBQUVHLEtBQWQsQ0FBcUJqUCxDQUFBQSxFQUExQixDQUZnQztBQUFBLGFBaENuQjtBQUFBLFlBcUNqQixJQUFJQSxFQUFBLEtBQU9wQixTQUFYO0FBQUEsZ0JBQXNCdVMsQ0FBQSxDQUFFblIsRUFBRixHQUFPQSxFQUFQLENBckNMO0FBQUEsWUF1Q2pCaUgsSUFBQSxDQUFLekQsUUFBTCxDQUFjdEQsSUFBZCxDQUFtQmlSLENBQW5CLENBdkNpQixDQUFBO0FBQUEsU0FEcUI7QUFBQSxLQTdLaEI7QUFBQSxJQXlOOUJ4QyxVQUFBLENBQVcvSSxDQUFYLEVBQWM7QUFBQSxRQUNWLE9BQU90QixJQUFBLENBQUtzQyxHQUFMLENBQVMsSUFBQSxDQUFLVixPQUFMLENBQWF3RyxPQUF0QixFQUErQnBJLElBQUEsQ0FBS3FDLEdBQUwsQ0FBU3JDLElBQUEsQ0FBSzhGLEtBQUwsQ0FBVyxDQUFDeEUsQ0FBWixDQUFULEVBQXlCLElBQUEsQ0FBS00sT0FBTCxDQUFheUcsT0FBYixHQUF1QixDQUFoRCxDQUEvQixDQUFQLENBRFU7QUFBQSxLQXpOZ0I7QUFBQSxJQTZOOUJxQixRQUFBLENBQVMvQixNQUFULEVBQWlCaUMsSUFBakIsRUFBdUI7QUFBQSxRQUNuQixNQUFNUCxRQUFBLEdBQVcsRUFBakIsQ0FEbUI7QUFBQSxRQUVuQixNQUFNLENBQUNkLE1BQUQsRUFBUzdILE1BQVQsRUFBaUIrSCxNQUFqQixFQUF5QkgsU0FBekIsQ0FBQSxHQUFzQyxJQUFLMUcsQ0FBQUEsT0FBakQsQ0FGbUI7QUFBQSxRQUduQixNQUFNMEMsQ0FBQSxHQUFJaUUsTUFBQSxJQUFVN0gsTUFBQSxHQUFTVixJQUFBLENBQUtxTCxHQUFMLENBQVMsQ0FBVCxFQUFZekIsSUFBWixDQUFULENBQXBCLENBSG1CO0FBQUE7QUFNbkIsUUFBQSxLQUFLLElBQUlwTyxDQUFBLEdBQUksQ0FBUixFQUFXQSxDQUFBLEdBQUltTSxNQUFBLENBQU9sTSxNQUEzQixFQUFtQ0QsQ0FBQSxFQUFuQyxFQUF3QztBQUFBLFlBQ3BDLE1BQU1nTSxDQUFBLEdBQUlHLE1BQUEsQ0FBT25NLENBQVAsQ0FBVixDQURvQztBQUFBO0FBR3BDLFlBQUEsSUFBSWdNLENBQUEsQ0FBRW9DLElBQUYsSUFBVUEsSUFBZDtBQUFBLGdCQUFvQixTQUhnQjtBQUFBLFlBSXBDcEMsQ0FBQSxDQUFFb0MsSUFBRixHQUFTQSxJQUFULENBSm9DO0FBQUE7QUFPcEMsWUFBTVEsTUFBQUEsSUFBQSxHQUFPLElBQUtwQixDQUFBQSxLQUFMLENBQVdZLElBQUEsR0FBTyxDQUFsQixDQUFiLENBUG9DO0FBQUEsWUFRcEMsTUFBTW1ELFdBQUEsR0FBYzNDLElBQUEsQ0FBS3RELE1BQUwsQ0FBWVUsQ0FBQSxDQUFFcEcsQ0FBZCxFQUFpQm9HLENBQUEsQ0FBRW5HLENBQW5CLEVBQXNCaUQsQ0FBdEIsQ0FBcEIsQ0FSb0M7QUFBQSxZQVVwQyxNQUFNMEksZUFBQSxHQUFrQnhGLENBQUEsQ0FBRWlELFNBQUYsSUFBZSxDQUF2QyxDQVZvQztBQUFBLFlBV3BDLElBQUlBLFNBQUEsR0FBWXVDLGVBQWhCLENBWG9DO0FBQUE7QUFjcEMsWUFBV0MsS0FBQUEsTUFBQUEsVUFBWCxJQUF5QkYsV0FBekIsRUFBc0M7QUFBQSxnQkFDbEMsTUFBTUcsQ0FBQSxHQUFJOUMsSUFBQSxDQUFLekMsTUFBTCxDQUFZc0YsVUFBWixDQUFWLENBRGtDO0FBQUE7QUFHbEMsZ0JBQUEsSUFBSUMsQ0FBQSxDQUFFdEQsSUFBRixHQUFTQSxJQUFiO0FBQUEsb0JBQW1CYSxTQUFBLElBQWF5QyxDQUFBLENBQUV6QyxTQUFGLElBQWUsQ0FBNUIsQ0FIZTtBQUFBLGFBZEY7QUFBQTtBQXFCcEMsWUFBQSxJQUFJQSxTQUFBLEdBQVl1QyxlQUFaLElBQStCdkMsU0FBQSxJQUFhbkMsU0FBaEQsRUFBMkQ7QUFBQSxnQkFDdkQsSUFBSTZFLEVBQUEsR0FBSzNGLENBQUEsQ0FBRXBHLENBQUYsR0FBTTRMLGVBQWYsQ0FEdUQ7QUFBQSxnQkFFdkQsSUFBSUksRUFBQSxHQUFLNUYsQ0FBQSxDQUFFbkcsQ0FBRixHQUFNMkwsZUFBZixDQUZ1RDtBQUFBLGdCQUl2RCxJQUFJSyxpQkFBQSxHQUFvQjVFLE1BQUEsSUFBVXVFLGVBQUEsR0FBa0IsQ0FBNUIsR0FBZ0MsSUFBS00sQ0FBQUEsSUFBTCxDQUFVOUYsQ0FBVixFQUFhLElBQWIsQ0FBaEMsR0FBcUQsSUFBN0UsQ0FKdUQ7QUFBQTtBQU92RCxnQkFBQSxNQUFNOUwsRUFBQSxHQUFNLENBQUFGLENBQUEsSUFBSyxDQUFMLEtBQVdvTyxJQUFBLEdBQU8sQ0FBUCxDQUFaLEdBQXdCLEtBQUtqQyxNQUFMLENBQVlsTSxNQUEvQyxDQVB1RDtBQUFBLGdCQVN2RCxLQUFXd1IsTUFBQUEsVUFBWCxJQUF5QkYsV0FBekIsRUFBc0M7QUFBQSxvQkFDbEMsTUFBTUcsQ0FBQSxHQUFJOUMsSUFBQSxDQUFLekMsTUFBTCxDQUFZc0YsVUFBWixDQUFWLENBRGtDO0FBQUEsb0JBR2xDLElBQUlDLENBQUEsQ0FBRXRELElBQUYsSUFBVUEsSUFBZDtBQUFBLHdCQUFvQixTQUhjO0FBQUEsb0JBSWxDc0QsQ0FBQSxDQUFFdEQsSUFBRixHQUFTQSxJQUFULENBSmtDO0FBQUE7QUFNbEMsb0JBQUEsTUFBTTJELFVBQUEsR0FBYUwsQ0FBQSxDQUFFekMsU0FBRixJQUFlLENBQWxDLENBTmtDO0FBQUEsb0JBT2xDMEMsRUFBQSxJQUFNRCxDQUFBLENBQUU5TCxDQUFGLEdBQU1tTSxVQUFaLENBUGtDO0FBQUE7QUFRbEMsb0JBQUFILEVBQUEsSUFBTUYsQ0FBQSxDQUFFN0wsQ0FBRixHQUFNa00sVUFBWixDQVJrQztBQUFBLG9CQVVsQ0wsQ0FBQSxDQUFFM0IsUUFBRixHQUFhN1AsRUFBYixDQVZrQztBQUFBLG9CQVlsQyxJQUFJK00sTUFBSixFQUFZO0FBQUEsd0JBQ1IsSUFBSSxDQUFDNEUsaUJBQUw7QUFBQSw0QkFBd0JBLGlCQUFBLEdBQW9CLElBQUtDLENBQUFBLElBQUwsQ0FBVTlGLENBQVYsRUFBYSxJQUFiLENBQXBCLENBRGhCO0FBQUEsd0JBRVJpQixNQUFBLENBQU80RSxpQkFBUCxFQUEwQixLQUFLQyxJQUFMLENBQVVKLENBQVYsQ0FBMUIsQ0FGUSxDQUFBO0FBQUEscUJBWnNCO0FBQUEsaUJBVGlCO0FBQUEsZ0JBMkJ2RDFGLENBQUEsQ0FBRStELFFBQUYsR0FBYTdQLEVBQWIsQ0EzQnVEO0FBQUEsZ0JBNEJ2RDJOLFFBQUEsQ0FBU3pOLElBQVQsQ0FBYzRSLGFBQUEsQ0FBY0wsRUFBQSxHQUFLMUMsU0FBbkIsRUFBOEIyQyxFQUFBLEdBQUszQyxTQUFuQyxFQUE4Qy9PLEVBQTlDLEVBQWtEK08sU0FBbEQsRUFBNkQ0QyxpQkFBN0QsQ0FBZCxDQTVCdUQsQ0FBQTtBQUFBLGFBQTNELE1BOEJPO0FBQUE7QUFDSCxnQkFBQWhFLFFBQUEsQ0FBU3pOLElBQVQsQ0FBYzRMLENBQWQsQ0FERyxDQUFBO0FBQUEsZ0JBR0gsSUFBSWlELFNBQUEsR0FBWSxDQUFoQixFQUFtQjtBQUFBLG9CQUNmLEtBQVd3QyxNQUFBQSxVQUFYLElBQXlCRixXQUF6QixFQUFzQztBQUFBLHdCQUNsQyxNQUFNRyxDQUFBLEdBQUk5QyxJQUFBLENBQUt6QyxNQUFMLENBQVlzRixVQUFaLENBQVYsQ0FEa0M7QUFBQSx3QkFFbEMsSUFBSUMsQ0FBQSxDQUFFdEQsSUFBRixJQUFVQSxJQUFkO0FBQUEsNEJBQW9CLFNBRmM7QUFBQSx3QkFHbENzRCxDQUFBLENBQUV0RCxJQUFGLEdBQVNBLElBQVQsQ0FIa0M7QUFBQSx3QkFJbENQLFFBQUEsQ0FBU3pOLElBQVQsQ0FBY3NSLENBQWQsQ0FKa0MsQ0FBQTtBQUFBLHFCQUR2QjtBQUFBLGlCQUhoQjtBQUFBLGFBbkQ2QjtBQUFBLFNBTnJCO0FBQUEsUUF1RW5CLE9BQU83RCxRQUFQLENBdkVtQjtBQUFBLEtBN05PO0FBQUE7QUF3UzlCLElBQUEwQixZQUFBLENBQWFGLFNBQWIsRUFBd0I7QUFBQSxRQUNwQixPQUFRQSxTQUFBLEdBQVksSUFBQSxDQUFLbEQsTUFBTCxDQUFZbE0sTUFBekIsSUFBb0MsQ0FBM0MsQ0FEb0I7QUFBQSxLQXhTTTtBQUFBO0FBNlM5QixJQUFBd1AsY0FBQSxDQUFlSixTQUFmLEVBQTBCO0FBQUEsUUFDdEIsT0FBUSxDQUFBQSxTQUFBLEdBQVksSUFBQSxDQUFLbEQsTUFBTCxDQUFZbE0sTUFBeEIsSUFBa0MsRUFBMUMsQ0FEc0I7QUFBQSxLQTdTSTtBQUFBLElBaVQ5QjZSLElBQUEsQ0FBS3JNLEtBQUwsRUFBWXdNLEtBQVosRUFBbUI7QUFBQSxRQUNmLElBQUl4TSxLQUFBLENBQU13SixTQUFWLEVBQXFCO0FBQUEsWUFDakIsT0FBT2dELEtBQUEsR0FBUTNFLFFBQUEsQ0FBTyxFQUFQLEVBQVc3SCxLQUFBLENBQU1MLFVBQWpCLENBQVIsR0FBdUNLLEtBQUEsQ0FBTUwsVUFBcEQsQ0FEaUI7QUFBQSxTQUROO0FBQUEsUUFJZixNQUFNOE0sUUFBQSxHQUFXLElBQUsvRixDQUFBQSxNQUFMLENBQVkxRyxLQUFBLENBQU0wSixLQUFsQixDQUF5Qi9KLENBQUFBLFVBQTFDLENBSmU7QUFBQSxRQUtmLE1BQU0vRSxNQUFBLEdBQVMsSUFBSytGLENBQUFBLE9BQUwsQ0FBYWhGLEdBQWIsQ0FBaUI4USxRQUFqQixDQUFmLENBTGU7QUFBQSxRQU1mLE9BQU9ELEtBQUEsSUFBUzVSLE1BQUEsS0FBVzZSLFFBQXBCLEdBQStCNUUsUUFBQSxDQUFPLEVBQVAsRUFBV2pOLE1BQVgsQ0FBL0IsR0FBb0RBLE1BQTNELENBTmU7QUFBQSxLQWpUVztBQUFBLENBeEJsQztBQW1WQSxTQUFTMlIsYUFBVCxDQUF1QnBNLENBQXZCLEVBQTBCQyxDQUExQixFQUE2QjNGLEVBQTdCLEVBQWlDK08sU0FBakMsRUFBNEM3SixVQUE1QyxFQUF3RDtBQUFBLElBQ3BELE9BQU87QUFBQSxRQUNIUSxDQUFBLEVBQUd1SCxNQUFBLENBQU92SCxDQUFQLENBREE7QUFBQTtBQUVILFFBQUFDLENBQUEsRUFBR3NILE1BQUEsQ0FBT3RILENBQVAsQ0FGQTtBQUFBLFFBR0h1SSxJQUFBLEVBQU01SCxRQUhIO0FBQUE7QUFJSCxRQUFBdEcsRUFKRztBQUFBO0FBS0gsUUFBQTZQLFFBQUEsRUFBVSxDQUFDLENBTFI7QUFBQTtBQU1ILFFBQUFkLFNBTkc7QUFBQSxRQU9IN0osVUFQRztBQUFBLEtBQVAsQ0FEb0Q7QUFBQSxDQW5WeEQ7QUErVkEsU0FBUzBJLGtCQUFULENBQTRCOUIsQ0FBNUIsRUFBK0I5TCxFQUEvQixFQUFtQztBQUFBLElBQy9CLE1BQU0sQ0FBQzBGLENBQUQsRUFBSUMsQ0FBSixDQUFTbUcsR0FBQUEsQ0FBQSxDQUFFcEksUUFBRixDQUFXRSxXQUExQixDQUQrQjtBQUFBLElBRS9CLE9BQU87QUFBQSxRQUNIOEIsQ0FBQSxFQUFHdUgsTUFBQSxDQUFPMkIsSUFBQSxDQUFLbEosQ0FBTCxDQUFQLENBREE7QUFBQTtBQUVILFFBQUFDLENBQUEsRUFBR3NILE1BQUEsQ0FBTzRCLElBQUEsQ0FBS2xKLENBQUwsQ0FBUCxDQUZBO0FBQUEsUUFHSHVJLElBQUEsRUFBTTVILFFBSEg7QUFBQTtBQUlILFFBQUEySSxLQUFBLEVBQU9qUCxFQUpKO0FBQUE7QUFLSCxRQUFBNlAsUUFBQSxFQUFVLENBQUMsQ0FBQTtBQUxSLEtBQVAsQ0FGK0I7QUFBQSxDQS9WbkM7QUEwV0EsU0FBU2IsY0FBVCxDQUF3QjZCLE9BQXhCLEVBQWlDO0FBQUEsSUFDN0IsT0FBTztBQUFBLFFBQ0h0TixJQUFBLEVBQU0sU0FESDtBQUFBLFFBRUh2RCxFQUFBLEVBQUk2USxPQUFBLENBQVE3USxFQUZUO0FBQUEsUUFHSGtGLFVBQUEsRUFBWWdNLG9CQUFBLENBQXFCTCxPQUFyQixDQUhUO0FBQUEsUUFJSG5OLFFBQUEsRUFBVTtBQUFBLFlBQ05ILElBQUEsRUFBTSxPQURBO0FBQUEsWUFFTkssV0FBQSxFQUFhO0FBQUEsZ0JBQUNxTyxJQUFBLENBQUtwQixPQUFBLENBQVFuTCxDQUFiLENBQUQ7QUFBQSxnQkFBa0J3TSxJQUFBLENBQUtyQixPQUFBLENBQVFsTCxDQUFiLENBQWxCO0FBQUEsYUFGUDtBQUFBLFNBSlA7QUFBQSxLQUFQLENBRDZCO0FBQUEsQ0ExV2pDO0FBc1hBLFNBQVN1TCxvQkFBVCxDQUE4QkwsT0FBOUIsRUFBdUM7QUFBQSxJQUNuQyxNQUFNaEksS0FBQSxHQUFRZ0ksT0FBQSxDQUFROUIsU0FBdEIsQ0FEbUM7QUFBQSxJQUVuQyxNQUFNb0QsTUFBQSxHQUNGdEosS0FBQSxJQUFTLEtBQVQsR0FBaUIsQ0FBQSxHQUFHdkUsSUFBQSxDQUFLOE0sS0FBTCxDQUFXdkksS0FBQSxHQUFRLElBQW5CLENBQUgsRUFBOEIsQ0FBOUIsQ0FBakIsR0FDQUEsS0FBQSxJQUFTLElBQVQsR0FBZ0IsSUFBR3ZFLElBQUEsQ0FBSzhNLEtBQUwsQ0FBV3ZJLEtBQUEsR0FBUSxHQUFuQixJQUEwQixFQUE3QixFQUFrQyxDQUFsQyxDQUFoQixHQUF1REEsS0FGM0QsQ0FGbUM7QUFBQSxJQUtuQyxPQUFPdUUsUUFBQSxDQUFPQSxRQUFBLENBQU8sRUFBUCxFQUFXeUQsT0FBQSxDQUFRM0wsVUFBbkIsQ0FBUCxFQUF1QztBQUFBLFFBQzFDMkwsT0FBQSxFQUFTLElBRGlDO0FBQUEsUUFFMUNILFVBQUEsRUFBWUcsT0FBQSxDQUFRN1EsRUFGc0I7QUFBQSxRQUcxQzhRLFdBQUEsRUFBYWpJLEtBSDZCO0FBQUEsUUFJMUN1Six1QkFBQSxFQUF5QkQsTUFKaUI7QUFBQSxLQUF2QyxDQUFQLENBTG1DO0FBQUEsQ0F0WHZDO0FBb1lBO0FBQUEsU0FBU3ZELElBQVQsQ0FBY3lELEdBQWQsRUFBbUI7QUFBQSxJQUNmLE9BQU9BLEdBQUEsR0FBTSxHQUFOLEdBQVksR0FBbkIsQ0FEZTtBQUFBLENBcFluQjtBQXVZQSxTQUFTeEQsSUFBVCxDQUFjeUQsR0FBZCxFQUFtQjtBQUFBLElBQ2YsTUFBTUMsR0FBQSxHQUFNak8sSUFBQSxDQUFLaU8sR0FBTCxDQUFTRCxHQUFBLEdBQU1oTyxJQUFBLENBQUtrTyxFQUFYLEdBQWdCLEdBQXpCLENBQVosQ0FEZTtBQUFBLElBRWYsTUFBTTdNLENBQUEsR0FBSyxHQUFBLEdBQU0sT0FBT3JCLElBQUEsQ0FBS3dGLEdBQUwsQ0FBVSxDQUFBLENBQUEsR0FBSXlJLEdBQUosU0FBZ0JBLEdBQUosQ0FBdEIsQ0FBUCxHQUF5Q2pPLElBQUEsQ0FBS2tPLEVBQS9ELENBRmU7QUFBQSxJQUdmLE9BQU83TSxDQUFBLEdBQUksQ0FBSixHQUFRLENBQVIsR0FBWUEsQ0FBQSxHQUFJLENBQUosR0FBUSxDQUFSLEdBQVlBLENBQS9CLENBSGU7QUFBQSxDQXZZbkI7QUE4WUE7QUFBQSxTQUFTc00sSUFBVCxDQUFjdk0sQ0FBZCxFQUFpQjtBQUFBLElBQ2IsT0FBUSxDQUFBQSxDQUFBLEdBQUksR0FBSixJQUFXLEdBQW5CLENBRGE7QUFBQSxDQTlZakI7QUFpWkEsU0FBU3dNLElBQVQsQ0FBY3ZNLENBQWQsRUFBaUI7QUFBQSxJQUNiLE1BQU1jLEVBQUEsR0FBTSxDQUFBLEdBQUEsR0FBTWQsQ0FBQSxHQUFJLEdBQVYsSUFBaUJyQixJQUFBLENBQUtrTyxFQUF2QixHQUE0QixHQUF2QyxDQURhO0FBQUEsSUFFYixPQUFPLEdBQU1sTyxHQUFBQSxJQUFBLENBQUttTyxJQUFMLENBQVVuTyxJQUFBLENBQUswRixHQUFMLENBQVN2RCxFQUFULENBQVYsQ0FBTixHQUFnQ25DLElBQUEsQ0FBS2tPLEVBQXJDLEdBQTBDLEVBQWpELENBRmE7QUFBQSxDQWpaakI7QUFzWkEsU0FBU3BGLFFBQVQsQ0FBZ0JzRixJQUFoQixFQUFzQkMsR0FBdEIsRUFBMkI7QUFBQSxJQUN2QixLQUFBLE1BQVczUyxFQUFYLElBQWlCMlMsR0FBakI7QUFBQSxRQUFzQkQsSUFBQSxDQUFLMVMsRUFBTCxDQUFBLEdBQVcyUyxHQUFBLENBQUkzUyxFQUFKLENBQVgsQ0FEQztBQUFBLElBRXZCLE9BQU8wUyxJQUFQLENBRnVCO0FBQUEsQ0F0WjNCO0FBMlpBLFNBQVN4RyxJQUFULENBQWNKLENBQWQsRUFBaUI7QUFBQSxJQUNiLE9BQU9BLENBQUEsQ0FBRXBHLENBQVQsQ0FEYTtBQUFBLENBM1pqQjtBQThaQSxTQUFTeUcsSUFBVCxDQUFjTCxDQUFkLEVBQWlCO0FBQUEsSUFDYixPQUFPQSxDQUFBLENBQUVuRyxDQUFULENBRGE7QUFBQTs7QUMzWmpCO0FBQWUsU0FBU2lOLFFBQVQsQ0FBa0J0SixNQUFsQixFQUEwQnVKLEtBQTFCLEVBQWlDQyxJQUFqQyxFQUF1Q0MsV0FBdkMsRUFBb0Q7QUFBQSxJQUMvRCxJQUFJQyxTQUFBLEdBQVlELFdBQWhCLENBRCtEO0FBQUEsSUFFL0QsSUFBSUUsR0FBQSxHQUFPSCxJQUFBLEdBQU9ELEtBQVIsSUFBa0IsQ0FBNUIsQ0FGK0Q7QUFBQSxJQUcvRCxJQUFJSyxXQUFBLEdBQWNKLElBQUEsR0FBT0QsS0FBekIsQ0FIK0Q7QUFBQSxJQUkvRCxJQUFJNUQsS0FBSixDQUorRDtBQUFBLElBTS9ELElBQUl4RCxFQUFBLEdBQUtuQyxNQUFBLENBQU91SixLQUFQLENBQVQsQ0FOK0Q7QUFBQSxJQU8vRCxJQUFJbkgsRUFBQSxHQUFLcEMsTUFBQSxDQUFPdUosS0FBQSxHQUFRLENBQWYsQ0FBVCxDQVArRDtBQUFBLElBUS9ELElBQUlsSCxFQUFBLEdBQUtyQyxNQUFBLENBQU93SixJQUFQLENBQVQsQ0FSK0Q7QUFBQSxJQVMvRCxJQUFJbEgsRUFBQSxHQUFLdEMsTUFBQSxDQUFPd0osSUFBQSxHQUFPLENBQWQsQ0FBVCxDQVQrRDtBQUFBLElBVy9ELEtBQUssSUFBSWhULENBQUEsR0FBSStTLEtBQUEsR0FBUSxDQUFoQixFQUFtQi9TLENBQUEsR0FBSWdULElBQTVCLEVBQWtDaFQsQ0FBQSxJQUFLLENBQXZDLEVBQTBDO0FBQUEsUUFDdEMsSUFBSXFULENBQUEsR0FBSUMsWUFBQSxDQUFhOUosTUFBQSxDQUFPeEosQ0FBUCxDQUFiLEVBQXdCd0osTUFBQSxDQUFPeEosQ0FBQSxHQUFJLENBQVgsQ0FBeEIsRUFBdUMyTCxFQUF2QyxFQUEyQ0MsRUFBM0MsRUFBK0NDLEVBQS9DLEVBQW1EQyxFQUFuRCxDQUFSLENBRHNDO0FBQUEsUUFHdEMsSUFBSXVILENBQUEsR0FBSUgsU0FBUixFQUFtQjtBQUFBLFlBQ2YvRCxLQUFBLEdBQVFuUCxDQUFSLENBRGU7QUFBQSxZQUVma1QsU0FBQSxHQUFZRyxDQUFaLENBRmU7QUFBQSxTQUFuQixNQUlPLElBQUlBLENBQUEsS0FBTUgsU0FBVixFQUFxQjtBQUFBO0FBSXhCO0FBQUE7QUFBQSxZQUFJSyxJQUFBQSxRQUFBLEdBQVcvTyxJQUFBLENBQUtDLEdBQUwsQ0FBU3pFLENBQUEsR0FBSW1ULEdBQWIsQ0FBZixDQUp3QjtBQUFBLFlBS3hCLElBQUlJLFFBQUEsR0FBV0gsV0FBZixFQUE0QjtBQUFBLGdCQUN4QmpFLEtBQUEsR0FBUW5QLENBQVIsQ0FEd0I7QUFBQSxnQkFFeEJvVCxXQUFBLEdBQWNHLFFBQWQsQ0FGd0I7QUFBQSxhQUxKO0FBQUEsU0FQVTtBQUFBLEtBWHFCO0FBQUEsSUE4Qi9ELElBQUlMLFNBQUEsR0FBWUQsV0FBaEIsRUFBNkI7QUFBQSxRQUN6QixJQUFJOUQsS0FBQSxHQUFRNEQsS0FBUixHQUFnQixDQUFwQjtBQUFBLFlBQXVCRCxRQUFBLENBQVN0SixNQUFULEVBQWlCdUosS0FBakIsRUFBd0I1RCxLQUF4QixFQUErQjhELFdBQS9CLENBREUsQ0FBQTtBQUFBLFFBRXpCekosTUFBQSxDQUFPMkYsS0FBQSxHQUFRLENBQWYsQ0FBQSxHQUFvQitELFNBQXBCLENBRnlCO0FBQUEsUUFHekIsSUFBSUYsSUFBQSxHQUFPN0QsS0FBUCxHQUFlLENBQW5CO0FBQUEsWUFBc0IyRCxRQUFBLENBQVN0SixNQUFULEVBQWlCMkYsS0FBakIsRUFBd0I2RCxJQUF4QixFQUE4QkMsV0FBOUIsQ0FIRyxDQUFBO0FBQUEsS0E5QmtDO0FBQUEsQ0FIbkU7QUF5Q0E7QUFBQSxTQUFTSyxZQUFULENBQXNCcEMsRUFBdEIsRUFBMEJDLEVBQTFCLEVBQThCdkwsQ0FBOUIsRUFBaUNDLENBQWpDLEVBQW9DZ0csRUFBcEMsRUFBd0NDLEVBQXhDLEVBQTRDO0FBQUEsSUFFeEMsSUFBSTdDLEVBQUEsR0FBSzRDLEVBQUEsR0FBS2pHLENBQWQsQ0FGd0M7QUFBQSxJQUd4QyxJQUFJc0QsRUFBQSxHQUFLNEMsRUFBQSxHQUFLakcsQ0FBZCxDQUh3QztBQUFBLElBS3hDLElBQUlvRCxFQUFBLEtBQU8sQ0FBUCxJQUFZQyxFQUFBLEtBQU8sQ0FBdkIsRUFBMEI7QUFBQSxRQUV0QixJQUFJc0IsQ0FBQSxHQUFLLENBQUMsQ0FBQTBHLEVBQUEsR0FBS3RMLENBQUwsSUFBVXFELEVBQVgsR0FBaUIsQ0FBQWtJLEVBQUEsR0FBS3RMLENBQUwsSUFBVXFELEVBQTNCLEtBQWtDRCxFQUFBLEdBQUtBLEVBQUwsR0FBVUMsRUFBQSxHQUFLQSxFQUFmLENBQTNDLENBRnNCO0FBQUEsUUFJdEIsSUFBSXNCLENBQUEsR0FBSSxDQUFSLEVBQVc7QUFBQSxZQUNQNUUsQ0FBQSxHQUFJaUcsRUFBSixDQURPO0FBQUEsWUFFUGhHLENBQUEsR0FBSWlHLEVBQUosQ0FGTztBQUFBLFNBQVgsTUFJTyxJQUFJdEIsQ0FBQSxHQUFJLENBQVIsRUFBVztBQUFBLFlBQ2Q1RSxDQUFBLElBQUtxRCxFQUFBLEdBQUt1QixDQUFWLENBRGM7QUFBQSxZQUVkM0UsQ0FBQSxJQUFLcUQsRUFBQSxHQUFLc0IsQ0FBVixDQUZjO0FBQUEsU0FSSTtBQUFBLEtBTGM7QUFBQSxJQW1CeEN2QixFQUFBLEdBQUtpSSxFQUFBLEdBQUt0TCxDQUFWLENBbkJ3QztBQUFBLElBb0J4Q3NELEVBQUEsR0FBS2lJLEVBQUEsR0FBS3RMLENBQVYsQ0FwQndDO0FBQUEsSUFzQnhDLE9BQU9vRCxFQUFBLEdBQUtBLEVBQUwsR0FBVUMsRUFBQSxHQUFLQSxFQUF0QixDQXRCd0M7QUFBQTs7QUN4QzdCLFNBQVNzSyxhQUFULENBQXVCdFQsRUFBdkIsRUFBMkJ1RCxJQUEzQixFQUFpQ2dRLElBQWpDLEVBQXVDcE8sSUFBdkMsRUFBNkM7QUFBQSxJQUN4RCxJQUFJTCxPQUFBLEdBQVU7QUFBQSxRQUNWOUUsRUFBQSxFQUFJLE9BQU9BLEVBQVAsS0FBYyxXQUFkLEdBQTRCLElBQTVCLEdBQW1DQSxFQUQ3QjtBQUFBLFFBRVZ1RCxJQUFBLEVBQU1BLElBRkk7QUFBQSxRQUdWRyxRQUFBLEVBQVU2UCxJQUhBO0FBQUEsUUFJVnBPLElBQUEsRUFBTUEsSUFKSTtBQUFBLFFBS1Z5RixJQUFBLEVBQU10RSxRQUxJO0FBQUEsUUFNVnVFLElBQUEsRUFBTXZFLFFBTkk7QUFBQSxRQU9Wd0UsSUFBQSxFQUFNLENBQUN4RSxRQVBHO0FBQUEsUUFRVnlFLElBQUEsRUFBTSxDQUFDekUsUUFSRztBQUFBLEtBQWQsQ0FEd0Q7QUFBQSxJQVd4RGtOLFFBQUEsQ0FBUzFPLE9BQVQsQ0FYd0QsQ0FBQTtBQUFBLElBWXhELE9BQU9BLE9BQVAsQ0Fad0Q7QUFBQSxDQUQ1RDtBQWdCQSxTQUFTME8sUUFBVCxDQUFrQjFPLE9BQWxCLEVBQTJCO0FBQUEsSUFDdkIsSUFBSXlPLElBQUEsR0FBT3pPLE9BQUEsQ0FBUXBCLFFBQW5CLENBRHVCO0FBQUEsSUFFdkIsSUFBSUgsSUFBQSxHQUFPdUIsT0FBQSxDQUFRdkIsSUFBbkIsQ0FGdUI7QUFBQSxJQUl2QixJQUFJQSxJQUFBLEtBQVMsT0FBVCxJQUFvQkEsSUFBQSxLQUFTLFlBQTdCLElBQTZDQSxJQUFBLEtBQVMsWUFBMUQsRUFBd0U7QUFBQSxRQUNwRWtRLFlBQUEsQ0FBYTNPLE9BQWIsRUFBc0J5TyxJQUF0QixDQURvRSxDQUFBO0FBQUEsS0FBeEUsTUFHTyxJQUFJaFEsSUFBQSxLQUFTLFNBQVQsSUFBc0JBLElBQUEsS0FBUyxpQkFBbkMsRUFBc0Q7QUFBQSxRQUN6RCxLQUFLLElBQUl6RCxDQUFBLEdBQUksQ0FBUixFQUFXQSxDQUFBLEdBQUl5VCxJQUFBLENBQUt4VCxNQUF6QixFQUFpQ0QsQ0FBQSxFQUFqQyxFQUFzQztBQUFBLFlBQ2xDMlQsWUFBQSxDQUFhM08sT0FBYixFQUFzQnlPLElBQUEsQ0FBS3pULENBQUwsQ0FBdEIsQ0FEa0MsQ0FBQTtBQUFBLFNBRG1CO0FBQUEsS0FBdEQsTUFLQSxJQUFJeUQsSUFBQSxLQUFTLGNBQWIsRUFBNkI7QUFBQSxRQUNoQyxLQUFLekQsQ0FBQSxHQUFJLENBQVQsRUFBWUEsQ0FBQSxHQUFJeVQsSUFBQSxDQUFLeFQsTUFBckIsRUFBNkJELENBQUEsRUFBN0IsRUFBa0M7QUFBQSxZQUM5QixLQUFLLElBQUlzRSxDQUFBLEdBQUksQ0FBUixFQUFXQSxDQUFBLEdBQUltUCxJQUFBLENBQUt6VCxDQUFMLENBQUEsQ0FBUUMsTUFBNUIsRUFBb0NxRSxDQUFBLEVBQXBDLEVBQXlDO0FBQUEsZ0JBQ3JDcVAsWUFBQSxDQUFhM08sT0FBYixFQUFzQnlPLElBQUEsQ0FBS3pULENBQUwsQ0FBQSxDQUFRc0UsQ0FBUixDQUF0QixDQURxQyxDQUFBO0FBQUEsYUFEWDtBQUFBLFNBREY7QUFBQSxLQVpiO0FBQUEsQ0FoQjNCO0FBcUNBLFNBQVNxUCxZQUFULENBQXNCM08sT0FBdEIsRUFBK0J5TyxJQUEvQixFQUFxQztBQUFBLElBQ2pDLEtBQUssSUFBSXpULENBQUEsR0FBSSxDQUFSLEVBQVdBLENBQUEsR0FBSXlULElBQUEsQ0FBS3hULE1BQXpCLEVBQWlDRCxDQUFBLElBQUssQ0FBdEMsRUFBeUM7QUFBQSxRQUNyQ2dGLE9BQUEsQ0FBUThGLElBQVIsR0FBZXRHLElBQUEsQ0FBS3FDLEdBQUwsQ0FBUzdCLE9BQUEsQ0FBUThGLElBQWpCLEVBQXVCMkksSUFBQSxDQUFLelQsQ0FBTCxDQUF2QixDQUFmLENBRHFDO0FBQUEsUUFFckNnRixPQUFBLENBQVErRixJQUFSLEdBQWV2RyxJQUFBLENBQUtxQyxHQUFMLENBQVM3QixPQUFBLENBQVErRixJQUFqQixFQUF1QjBJLElBQUEsQ0FBS3pULENBQUEsR0FBSSxDQUFULENBQXZCLENBQWYsQ0FGcUM7QUFBQSxRQUdyQ2dGLE9BQUEsQ0FBUWdHLElBQVIsR0FBZXhHLElBQUEsQ0FBS3NDLEdBQUwsQ0FBUzlCLE9BQUEsQ0FBUWdHLElBQWpCLEVBQXVCeUksSUFBQSxDQUFLelQsQ0FBTCxDQUF2QixDQUFmLENBSHFDO0FBQUEsUUFJckNnRixPQUFBLENBQVFpRyxJQUFSLEdBQWV6RyxJQUFBLENBQUtzQyxHQUFMLENBQVM5QixPQUFBLENBQVFpRyxJQUFqQixFQUF1QndJLElBQUEsQ0FBS3pULENBQUEsR0FBSSxDQUFULENBQXZCLENBQWYsQ0FKcUM7QUFBQSxLQURSO0FBQUE7O0FDL0JyQztBQUFlLFNBQVM0VCxPQUFULENBQWlCQyxJQUFqQixFQUF1QnpOLE9BQXZCLEVBQWdDO0FBQUEsSUFDM0MsSUFBSTFDLFFBQUEsR0FBVyxFQUFmLENBRDJDO0FBQUEsSUFFM0MsSUFBSW1RLElBQUEsQ0FBS3BRLElBQUwsS0FBYyxtQkFBbEIsRUFBdUM7QUFBQSxRQUNuQyxLQUFLLElBQUl6RCxDQUFBLEdBQUksQ0FBUixFQUFXQSxDQUFBLEdBQUk2VCxJQUFBLENBQUtuUSxRQUFMLENBQWN6RCxNQUFsQyxFQUEwQ0QsQ0FBQSxFQUExQyxFQUErQztBQUFBLFlBQzNDOFQsY0FBQSxDQUFlcFEsUUFBZixFQUF5Qm1RLElBQUEsQ0FBS25RLFFBQUwsQ0FBYzFELENBQWQsQ0FBekIsRUFBMkNvRyxPQUEzQyxFQUFvRHBHLENBQXBELENBRDJDLENBQUE7QUFBQSxTQURaO0FBQUEsS0FBdkMsTUFLTyxJQUFJNlQsSUFBQSxDQUFLcFEsSUFBTCxLQUFjLFNBQWxCLEVBQTZCO0FBQUEsUUFDaENxUSxjQUFBLENBQWVwUSxRQUFmLEVBQXlCbVEsSUFBekIsRUFBK0J6TixPQUEvQixDQURnQyxDQUFBO0FBQUEsS0FBN0IsTUFHQTtBQUFBO0FBRUgsUUFBQTBOLGNBQUEsQ0FBZXBRLFFBQWYsRUFBeUIsRUFBQ0UsUUFBQSxFQUFVaVEsSUFBWCxFQUF6QixFQUEyQ3pOLE9BQTNDLENBRkcsQ0FBQTtBQUFBLEtBVm9DO0FBQUEsSUFlM0MsT0FBTzFDLFFBQVAsQ0FmMkM7QUFBQSxDQU4vQztBQXdCQSxTQUFTb1EsY0FBVCxDQUF3QnBRLFFBQXhCLEVBQWtDcVEsT0FBbEMsRUFBMkMzTixPQUEzQyxFQUFvRCtJLEtBQXBELEVBQTJEO0FBQUEsSUFDdkQsSUFBSSxDQUFDNEUsT0FBQSxDQUFRblEsUUFBYjtBQUFBLFFBQXVCLE9BRGdDO0FBQUEsSUFHdkQsSUFBSTRGLE1BQUEsR0FBU3VLLE9BQUEsQ0FBUW5RLFFBQVIsQ0FBaUJFLFdBQTlCLENBSHVEO0FBQUEsSUFJdkQsSUFBSUwsSUFBQSxHQUFPc1EsT0FBQSxDQUFRblEsUUFBUixDQUFpQkgsSUFBNUIsQ0FKdUQ7QUFBQSxJQUt2RCxJQUFJdVEsU0FBQSxHQUFZeFAsSUFBQSxDQUFLcUwsR0FBTCxDQUFTekosT0FBQSxDQUFRNE4sU0FBUixJQUFzQixNQUFLNU4sT0FBQSxDQUFReUcsT0FBYixJQUF3QnpHLE9BQUEsQ0FBUWxCLE1BQWpDLENBQTlCLEVBQXdFLENBQXhFLENBQWhCLENBTHVEO0FBQUEsSUFNdkQsSUFBSXRCLFFBQUEsR0FBVyxFQUFmLENBTnVEO0FBQUEsSUFPdkQsSUFBSTFELEVBQUEsR0FBSzZULE9BQUEsQ0FBUTdULEVBQWpCLENBUHVEO0FBQUEsSUFRdkQsSUFBSWtHLE9BQUEsQ0FBUTZOLFNBQVosRUFBdUI7QUFBQSxRQUNuQi9ULEVBQUEsR0FBSzZULE9BQUEsQ0FBUTNPLFVBQVIsQ0FBbUJnQixPQUFBLENBQVE2TixTQUEzQixDQUFMLENBRG1CO0FBQUEsS0FBdkIsTUFFTyxJQUFJN04sT0FBQSxDQUFRNEcsVUFBWixFQUF3QjtBQUFBLFFBQzNCOU0sRUFBQSxHQUFLaVAsS0FBQSxJQUFTLENBQWQsQ0FEMkI7QUFBQSxLQVZ3QjtBQUFBLElBYXZELElBQUkxTCxJQUFBLEtBQVMsT0FBYixFQUFzQjtBQUFBLFFBQ2xCeVEsWUFBQSxDQUFhMUssTUFBYixFQUFxQjVGLFFBQXJCLENBRGtCLENBQUE7QUFBQSxLQUF0QixNQUdPLElBQUlILElBQUEsS0FBUyxZQUFiLEVBQTJCO0FBQUEsUUFDOUIsS0FBSyxJQUFJekQsQ0FBQSxHQUFJLENBQVIsRUFBV0EsQ0FBQSxHQUFJd0osTUFBQSxDQUFPdkosTUFBM0IsRUFBbUNELENBQUEsRUFBbkMsRUFBd0M7QUFBQSxZQUNwQ2tVLFlBQUEsQ0FBYTFLLE1BQUEsQ0FBT3hKLENBQVAsQ0FBYixFQUF3QjRELFFBQXhCLENBRG9DLENBQUE7QUFBQSxTQURWO0FBQUEsS0FBM0IsTUFLQSxJQUFJSCxJQUFBLEtBQVMsWUFBYixFQUEyQjtBQUFBLFFBQzlCMFEsV0FBQSxDQUFZM0ssTUFBWixFQUFvQjVGLFFBQXBCLEVBQThCb1EsU0FBOUIsRUFBeUMsS0FBekMsQ0FEOEIsQ0FBQTtBQUFBLEtBQTNCLE1BR0EsSUFBSXZRLElBQUEsS0FBUyxpQkFBYixFQUFnQztBQUFBLFFBQ25DLElBQUkyQyxPQUFBLENBQVFnTyxXQUFaLEVBQXlCO0FBQUE7QUFFckIsWUFBQSxLQUFLcFUsQ0FBQSxHQUFJLENBQVQsRUFBWUEsQ0FBQSxHQUFJd0osTUFBQSxDQUFPdkosTUFBdkIsRUFBK0JELENBQUEsRUFBL0IsRUFBb0M7QUFBQSxnQkFDaEM0RCxRQUFBLEdBQVcsRUFBWCxDQURnQztBQUFBLGdCQUVoQ3VRLFdBQUEsQ0FBWTNLLE1BQUEsQ0FBT3hKLENBQVAsQ0FBWixFQUF1QjRELFFBQXZCLEVBQWlDb1EsU0FBakMsRUFBNEMsS0FBNUMsQ0FGZ0MsQ0FBQTtBQUFBLGdCQUdoQ3RRLFFBQUEsQ0FBU3RELElBQVQsQ0FBY29ULGFBQUEsQ0FBY3RULEVBQWQsRUFBa0IsWUFBbEIsRUFBZ0MwRCxRQUFoQyxFQUEwQ21RLE9BQUEsQ0FBUTNPLFVBQWxELENBQWQsQ0FIZ0MsQ0FBQTtBQUFBLGFBRmY7QUFBQSxZQU9yQixPQVBxQjtBQUFBLFNBQXpCLE1BUU87QUFBQSxZQUNIaVAsWUFBQSxDQUFhN0ssTUFBYixFQUFxQjVGLFFBQXJCLEVBQStCb1EsU0FBL0IsRUFBMEMsS0FBMUMsQ0FERyxDQUFBO0FBQUEsU0FUNEI7QUFBQSxLQUFoQyxNQWFBLElBQUl2USxJQUFBLEtBQVMsU0FBYixFQUF3QjtBQUFBLFFBQzNCNFEsWUFBQSxDQUFhN0ssTUFBYixFQUFxQjVGLFFBQXJCLEVBQStCb1EsU0FBL0IsRUFBMEMsSUFBMUMsQ0FEMkIsQ0FBQTtBQUFBLEtBQXhCLE1BR0EsSUFBSXZRLElBQUEsS0FBUyxjQUFiLEVBQTZCO0FBQUEsUUFDaEMsS0FBS3pELENBQUEsR0FBSSxDQUFULEVBQVlBLENBQUEsR0FBSXdKLE1BQUEsQ0FBT3ZKLE1BQXZCLEVBQStCRCxDQUFBLEVBQS9CLEVBQW9DO0FBQUEsWUFDaEMsSUFBSXNVLE9BQUEsR0FBVSxFQUFkLENBRGdDO0FBQUEsWUFFaENELFlBQUEsQ0FBYTdLLE1BQUEsQ0FBT3hKLENBQVAsQ0FBYixFQUF3QnNVLE9BQXhCLEVBQWlDTixTQUFqQyxFQUE0QyxJQUE1QyxDQUZnQyxDQUFBO0FBQUEsWUFHaENwUSxRQUFBLENBQVN4RCxJQUFULENBQWNrVSxPQUFkLENBSGdDLENBQUE7QUFBQSxTQURKO0FBQUEsS0FBN0IsTUFNQSxJQUFJN1EsSUFBQSxLQUFTLG9CQUFiLEVBQW1DO0FBQUEsUUFDdEMsS0FBS3pELENBQUEsR0FBSSxDQUFULEVBQVlBLENBQUEsR0FBSStULE9BQUEsQ0FBUW5RLFFBQVIsQ0FBaUJELFVBQWpCLENBQTRCMUQsTUFBNUMsRUFBb0RELENBQUEsRUFBcEQsRUFBeUQ7QUFBQSxZQUNyRDhULGNBQUEsQ0FBZXBRLFFBQWYsRUFBeUI7QUFBQSxnQkFDckJ4RCxFQUFBLEVBQUlBLEVBRGlCO0FBQUEsZ0JBRXJCMEQsUUFBQSxFQUFVbVEsT0FBQSxDQUFRblEsUUFBUixDQUFpQkQsVUFBakIsQ0FBNEIzRCxDQUE1QixDQUZXO0FBQUEsZ0JBR3JCb0YsVUFBQSxFQUFZMk8sT0FBQSxDQUFRM08sVUFIQztBQUFBLGFBQXpCLEVBSUdnQixPQUpILEVBSVkrSSxLQUpaLENBRHFELENBQUE7QUFBQSxTQURuQjtBQUFBLFFBUXRDLE9BUnNDO0FBQUEsS0FBbkMsTUFTQTtBQUFBLFFBQ0gsTUFBTSxJQUFJUSxLQUFKLENBQVUsMkNBQVYsQ0FBTixDQURHO0FBQUEsS0F2RGdEO0FBQUEsSUEyRHZEak0sUUFBQSxDQUFTdEQsSUFBVCxDQUFjb1QsYUFBQSxDQUFjdFQsRUFBZCxFQUFrQnVELElBQWxCLEVBQXdCRyxRQUF4QixFQUFrQ21RLE9BQUEsQ0FBUTNPLFVBQTFDLENBQWQsQ0EzRHVELENBQUE7QUFBQSxDQXhCM0Q7QUFzRkEsU0FBUzhPLFlBQVQsQ0FBc0IxSyxNQUF0QixFQUE4QnBDLEdBQTlCLEVBQW1DO0FBQUEsSUFDL0JBLEdBQUEsQ0FBSWhILElBQUosQ0FBU21VLFFBQUEsQ0FBUy9LLE1BQUEsQ0FBTyxDQUFQLENBQVQsQ0FBVCxDQUQrQixDQUFBO0FBQUEsSUFFL0JwQyxHQUFBLENBQUloSCxJQUFKLENBQVNvVSxRQUFBLENBQVNoTCxNQUFBLENBQU8sQ0FBUCxDQUFULENBQVQsQ0FGK0IsQ0FBQTtBQUFBLElBRy9CcEMsR0FBQSxDQUFJaEgsSUFBSixDQUFTLENBQVQsQ0FIK0IsQ0FBQTtBQUFBLENBdEZuQztBQTRGQSxTQUFTK1QsV0FBVCxDQUFxQmxRLElBQXJCLEVBQTJCbUQsR0FBM0IsRUFBZ0M0TSxTQUFoQyxFQUEyQ1MsU0FBM0MsRUFBc0Q7QUFBQSxJQUNsRCxJQUFJQyxFQUFKLEVBQVFDLEVBQVIsQ0FEa0Q7QUFBQSxJQUVsRCxJQUFJQyxJQUFBLEdBQU8sQ0FBWCxDQUZrRDtBQUFBLElBSWxELEtBQUssSUFBSXRRLENBQUEsR0FBSSxDQUFSLEVBQVdBLENBQUEsR0FBSUwsSUFBQSxDQUFLaEUsTUFBekIsRUFBaUNxRSxDQUFBLEVBQWpDLEVBQXNDO0FBQUEsUUFDbEMsSUFBSXNCLENBQUEsR0FBSTJPLFFBQUEsQ0FBU3RRLElBQUEsQ0FBS0ssQ0FBTCxDQUFRLENBQUEsQ0FBUixDQUFULENBQVIsQ0FEa0M7QUFBQSxRQUVsQyxJQUFJdUIsQ0FBQSxHQUFJMk8sUUFBQSxDQUFTdlEsSUFBQSxDQUFLSyxDQUFMLENBQVEsQ0FBQSxDQUFSLENBQVQsQ0FBUixDQUZrQztBQUFBLFFBSWxDOEMsR0FBQSxDQUFJaEgsSUFBSixDQUFTd0YsQ0FBVCxDQUprQyxDQUFBO0FBQUEsUUFLbEN3QixHQUFBLENBQUloSCxJQUFKLENBQVN5RixDQUFULENBTGtDLENBQUE7QUFBQSxRQU1sQ3VCLEdBQUEsQ0FBSWhILElBQUosQ0FBUyxDQUFULENBTmtDLENBQUE7QUFBQSxRQVFsQyxJQUFJa0UsQ0FBQSxHQUFJLENBQVIsRUFBVztBQUFBLFlBQ1AsSUFBSW1RLFNBQUosRUFBZTtBQUFBLGdCQUNYRyxJQUFBLElBQVMsQ0FBQUYsRUFBQSxHQUFLN08sQ0FBTCxHQUFTRCxDQUFBLEdBQUkrTyxFQUFiLElBQW1CLENBQTVCLENBQUE7QUFEVyxhQUFmLE1BRU87QUFBQSxnQkFDSEMsSUFBQSxJQUFRcFEsSUFBQSxDQUFLNEYsSUFBTCxDQUFVNUYsSUFBQSxDQUFLcUwsR0FBTCxDQUFTakssQ0FBQSxHQUFJOE8sRUFBYixFQUFpQixDQUFqQixDQUFBLEdBQXNCbFEsSUFBQSxDQUFLcUwsR0FBTCxDQUFTaEssQ0FBQSxHQUFJOE8sRUFBYixFQUFpQixDQUFqQixDQUFoQyxDQUFSLENBQUE7QUFERyxhQUhBO0FBQUEsU0FSdUI7QUFBQSxRQWVsQ0QsRUFBQSxHQUFLOU8sQ0FBTCxDQWZrQztBQUFBLFFBZ0JsQytPLEVBQUEsR0FBSzlPLENBQUwsQ0FoQmtDO0FBQUEsS0FKWTtBQUFBLElBdUJsRCxJQUFJbU4sSUFBQSxHQUFPNUwsR0FBQSxDQUFJbkgsTUFBSixHQUFhLENBQXhCLENBdkJrRDtBQUFBLElBd0JsRG1ILEdBQUEsQ0FBSSxDQUFKLENBQUEsR0FBUyxDQUFULENBeEJrRDtBQUFBLElBeUJsRDBMLFFBQUEsQ0FBUzFMLEdBQVQsRUFBYyxDQUFkLEVBQWlCNEwsSUFBakIsRUFBdUJnQixTQUF2QixDQXpCa0QsQ0FBQTtBQUFBLElBMEJsRDVNLEdBQUEsQ0FBSTRMLElBQUEsR0FBTyxDQUFYLENBQUEsR0FBZ0IsQ0FBaEIsQ0ExQmtEO0FBQUEsSUE0QmxENUwsR0FBQSxDQUFJd04sSUFBSixHQUFXcFEsSUFBQSxDQUFLQyxHQUFMLENBQVNtUSxJQUFULENBQVgsQ0E1QmtEO0FBQUEsSUE2QmxEeE4sR0FBQSxDQUFJeU4sS0FBSixHQUFZLENBQVosQ0E3QmtEO0FBQUEsSUE4QmxEek4sR0FBQSxDQUFJME4sR0FBSixHQUFVMU4sR0FBQSxDQUFJd04sSUFBZCxDQTlCa0Q7QUFBQSxDQTVGdEQ7QUE2SEEsU0FBU1AsWUFBVCxDQUFzQnRRLEtBQXRCLEVBQTZCcUQsR0FBN0IsRUFBa0M0TSxTQUFsQyxFQUE2Q1MsU0FBN0MsRUFBd0Q7QUFBQSxJQUNwRCxLQUFLLElBQUl6VSxDQUFBLEdBQUksQ0FBUixFQUFXQSxDQUFBLEdBQUkrRCxLQUFBLENBQU05RCxNQUExQixFQUFrQ0QsQ0FBQSxFQUFsQyxFQUF1QztBQUFBLFFBQ25DLElBQUl5VCxJQUFBLEdBQU8sRUFBWCxDQURtQztBQUFBLFFBRW5DVSxXQUFBLENBQVlwUSxLQUFBLENBQU0vRCxDQUFOLENBQVosRUFBc0J5VCxJQUF0QixFQUE0Qk8sU0FBNUIsRUFBdUNTLFNBQXZDLENBRm1DLENBQUE7QUFBQSxRQUduQ3JOLEdBQUEsQ0FBSWhILElBQUosQ0FBU3FULElBQVQsQ0FIbUMsQ0FBQTtBQUFBLEtBRGE7QUFBQSxDQTdIeEQ7QUFxSUEsU0FBU2MsUUFBVCxDQUFrQjNPLENBQWxCLEVBQXFCO0FBQUEsSUFDakIsT0FBT0EsQ0FBQSxHQUFJLEdBQUosR0FBVSxHQUFqQixDQURpQjtBQUFBLENBcklyQjtBQXlJQSxTQUFTNE8sUUFBVCxDQUFrQjNPLENBQWxCLEVBQXFCO0FBQUEsSUFDakIsSUFBSTRNLEdBQUEsR0FBTWpPLElBQUEsQ0FBS2lPLEdBQUwsQ0FBUzVNLENBQUEsR0FBSXJCLElBQUEsQ0FBS2tPLEVBQVQsR0FBYyxHQUF2QixDQUFWLENBRGlCO0FBQUEsSUFFakIsSUFBSS9MLEVBQUEsR0FBSyxHQUFBLEdBQU0sT0FBT25DLElBQUEsQ0FBS3dGLEdBQUwsQ0FBVSxDQUFBLENBQUEsR0FBSXlJLEdBQUosU0FBZ0JBLEdBQUosQ0FBdEIsQ0FBUCxHQUF5Q2pPLElBQUEsQ0FBS2tPLEVBQTdELENBRmlCO0FBQUEsSUFHakIsT0FBTy9MLEVBQUEsR0FBSyxDQUFMLEdBQVMsQ0FBVCxHQUFhQSxFQUFBLEdBQUssQ0FBTCxHQUFTLENBQVQsR0FBYUEsRUFBakMsQ0FIaUI7QUFBQTs7QUMvSHJCO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFlLFNBQVNvTyxJQUFULENBQWNyUixRQUFkLEVBQXdCc1IsS0FBeEIsRUFBK0JDLEVBQS9CLEVBQW1DQyxFQUFuQyxFQUF1Qy9KLElBQXZDLEVBQTZDZ0ssTUFBN0MsRUFBcURDLE1BQXJELEVBQTZEaFAsT0FBN0QsRUFBc0U7QUFBQSxJQUVqRjZPLEVBQUEsSUFBTUQsS0FBTixDQUZpRjtBQUFBLElBR2pGRSxFQUFBLElBQU1GLEtBQU4sQ0FIaUY7QUFBQSxJQUtqRixJQUFJRyxNQUFBLElBQVVGLEVBQVYsSUFBZ0JHLE1BQUEsR0FBU0YsRUFBN0I7QUFBQSxRQUFpQyxPQUFPeFIsUUFBUCxDQUFBO0FBQWpDLFNBQ0ssSUFBSTBSLE1BQUEsR0FBU0gsRUFBVCxJQUFlRSxNQUFBLElBQVVELEVBQTdCO0FBQUEsUUFBaUMsT0FBTyxJQUFQLENBTjJDO0FBQUE7QUFRakYsSUFBSUcsSUFBQUEsT0FBQSxHQUFVLEVBQWQsQ0FSaUY7QUFBQSxJQVVqRixLQUFLLElBQUlyVixDQUFBLEdBQUksQ0FBUixFQUFXQSxDQUFBLEdBQUkwRCxRQUFBLENBQVN6RCxNQUE3QixFQUFxQ0QsQ0FBQSxFQUFyQyxFQUEwQztBQUFBLFFBRXRDLElBQUlnRixPQUFBLEdBQVV0QixRQUFBLENBQVMxRCxDQUFULENBQWQsQ0FGc0M7QUFBQSxRQUd0QyxJQUFJNEQsUUFBQSxHQUFXb0IsT0FBQSxDQUFRcEIsUUFBdkIsQ0FIc0M7QUFBQSxRQUl0QyxJQUFJSCxJQUFBLEdBQU91QixPQUFBLENBQVF2QixJQUFuQixDQUpzQztBQUFBLFFBTXRDLElBQUlvRCxHQUFBLEdBQU1zRSxJQUFBLEtBQVMsQ0FBVCxHQUFhbkcsT0FBQSxDQUFROEYsSUFBckIsR0FBNEI5RixPQUFBLENBQVErRixJQUE5QyxDQU5zQztBQUFBLFFBT3RDLElBQUlqRSxHQUFBLEdBQU1xRSxJQUFBLEtBQVMsQ0FBVCxHQUFhbkcsT0FBQSxDQUFRZ0csSUFBckIsR0FBNEJoRyxPQUFBLENBQVFpRyxJQUE5QyxDQVBzQztBQUFBLFFBU3RDLElBQUlwRSxHQUFBLElBQU9vTyxFQUFQLElBQWFuTyxHQUFBLEdBQU1vTyxFQUF2QixFQUEyQjtBQUFBO0FBQ3ZCLFlBQUFHLE9BQUEsQ0FBUWpWLElBQVIsQ0FBYTRFLE9BQWIsQ0FEdUIsQ0FBQTtBQUFBLFlBRXZCLFNBRnVCO0FBQUEsU0FBM0IsTUFHTyxJQUFJOEIsR0FBQSxHQUFNbU8sRUFBTixJQUFZcE8sR0FBQSxJQUFPcU8sRUFBdkIsRUFBMkI7QUFBQTtBQUM5QixZQUQ4QixTQUFBO0FBQUEsU0FaSTtBQUFBLFFBZ0J0QyxJQUFJSSxXQUFBLEdBQWMsRUFBbEIsQ0FoQnNDO0FBQUEsUUFrQnRDLElBQUk3UixJQUFBLEtBQVMsT0FBVCxJQUFvQkEsSUFBQSxLQUFTLFlBQWpDLEVBQStDO0FBQUEsWUFDM0M4UixVQUFBLENBQVczUixRQUFYLEVBQXFCMFIsV0FBckIsRUFBa0NMLEVBQWxDLEVBQXNDQyxFQUF0QyxFQUEwQy9KLElBQTFDLENBRDJDLENBQUE7QUFBQSxTQUEvQyxNQUdPLElBQUkxSCxJQUFBLEtBQVMsWUFBYixFQUEyQjtBQUFBLFlBQzlCK1IsUUFBQSxDQUFTNVIsUUFBVCxFQUFtQjBSLFdBQW5CLEVBQWdDTCxFQUFoQyxFQUFvQ0MsRUFBcEMsRUFBd0MvSixJQUF4QyxFQUE4QyxLQUE5QyxFQUFxRC9FLE9BQUEsQ0FBUWdPLFdBQTdELENBRDhCLENBQUE7QUFBQSxTQUEzQixNQUdBLElBQUkzUSxJQUFBLEtBQVMsaUJBQWIsRUFBZ0M7QUFBQSxZQUNuQ2dTLFNBQUEsQ0FBVTdSLFFBQVYsRUFBb0IwUixXQUFwQixFQUFpQ0wsRUFBakMsRUFBcUNDLEVBQXJDLEVBQXlDL0osSUFBekMsRUFBK0MsS0FBL0MsQ0FEbUMsQ0FBQTtBQUFBLFNBQWhDLE1BR0EsSUFBSTFILElBQUEsS0FBUyxTQUFiLEVBQXdCO0FBQUEsWUFDM0JnUyxTQUFBLENBQVU3UixRQUFWLEVBQW9CMFIsV0FBcEIsRUFBaUNMLEVBQWpDLEVBQXFDQyxFQUFyQyxFQUF5Qy9KLElBQXpDLEVBQStDLElBQS9DLENBRDJCLENBQUE7QUFBQSxTQUF4QixNQUdBLElBQUkxSCxJQUFBLEtBQVMsY0FBYixFQUE2QjtBQUFBLFlBQ2hDLEtBQUssSUFBSWEsQ0FBQSxHQUFJLENBQVIsRUFBV0EsQ0FBQSxHQUFJVixRQUFBLENBQVMzRCxNQUE3QixFQUFxQ3FFLENBQUEsRUFBckMsRUFBMEM7QUFBQSxnQkFDdEMsSUFBSWdRLE9BQUEsR0FBVSxFQUFkLENBRHNDO0FBQUEsZ0JBRXRDbUIsU0FBQSxDQUFVN1IsUUFBQSxDQUFTVSxDQUFULENBQVYsRUFBdUJnUSxPQUF2QixFQUFnQ1csRUFBaEMsRUFBb0NDLEVBQXBDLEVBQXdDL0osSUFBeEMsRUFBOEMsSUFBOUMsQ0FGc0MsQ0FBQTtBQUFBLGdCQUd0QyxJQUFJbUosT0FBQSxDQUFRclUsTUFBWixFQUFvQjtBQUFBLG9CQUNoQnFWLFdBQUEsQ0FBWWxWLElBQVosQ0FBaUJrVSxPQUFqQixDQURnQixDQUFBO0FBQUEsaUJBSGtCO0FBQUEsYUFEVjtBQUFBLFNBOUJFO0FBQUEsUUF3Q3RDLElBQUlnQixXQUFBLENBQVlyVixNQUFoQixFQUF3QjtBQUFBLFlBQ3BCLElBQUltRyxPQUFBLENBQVFnTyxXQUFSLElBQXVCM1EsSUFBQSxLQUFTLFlBQXBDLEVBQWtEO0FBQUEsZ0JBQzlDLEtBQUthLENBQUEsR0FBSSxDQUFULEVBQVlBLENBQUEsR0FBSWdSLFdBQUEsQ0FBWXJWLE1BQTVCLEVBQW9DcUUsQ0FBQSxFQUFwQyxFQUF5QztBQUFBLG9CQUNyQytRLE9BQUEsQ0FBUWpWLElBQVIsQ0FBYW9ULGFBQUEsQ0FBY3hPLE9BQUEsQ0FBUTlFLEVBQXRCLEVBQTBCdUQsSUFBMUIsRUFBZ0M2UixXQUFBLENBQVloUixDQUFaLENBQWhDLEVBQWdEVSxPQUFBLENBQVFLLElBQXhELENBQWIsQ0FEcUMsQ0FBQTtBQUFBLGlCQURLO0FBQUEsZ0JBSTlDLFNBSjhDO0FBQUEsYUFEOUI7QUFBQSxZQVFwQixJQUFJNUIsSUFBQSxLQUFTLFlBQVQsSUFBeUJBLElBQUEsS0FBUyxpQkFBdEMsRUFBeUQ7QUFBQSxnQkFDckQsSUFBSTZSLFdBQUEsQ0FBWXJWLE1BQVosS0FBdUIsQ0FBM0IsRUFBOEI7QUFBQSxvQkFDMUJ3RCxJQUFBLEdBQU8sWUFBUCxDQUQwQjtBQUFBLG9CQUUxQjZSLFdBQUEsR0FBY0EsV0FBQSxDQUFZLENBQVosQ0FBZCxDQUYwQjtBQUFBLGlCQUE5QixNQUdPO0FBQUEsb0JBQ0g3UixJQUFBLEdBQU8saUJBQVAsQ0FERztBQUFBLGlCQUo4QztBQUFBLGFBUnJDO0FBQUEsWUFnQnBCLElBQUlBLElBQUEsS0FBUyxPQUFULElBQW9CQSxJQUFBLEtBQVMsWUFBakMsRUFBK0M7QUFBQSxnQkFDM0NBLElBQUEsR0FBTzZSLFdBQUEsQ0FBWXJWLE1BQVosS0FBdUIsQ0FBdkIsR0FBMkIsT0FBM0IsR0FBcUMsWUFBNUMsQ0FEMkM7QUFBQSxhQWhCM0I7QUFBQSxZQW9CcEJvVixPQUFBLENBQVFqVixJQUFSLENBQWFvVCxhQUFBLENBQWN4TyxPQUFBLENBQVE5RSxFQUF0QixFQUEwQnVELElBQTFCLEVBQWdDNlIsV0FBaEMsRUFBNkN0USxPQUFBLENBQVFLLElBQXJELENBQWIsQ0FwQm9CLENBQUE7QUFBQSxTQXhDYztBQUFBLEtBVnVDO0FBQUEsSUEwRWpGLE9BQU9nUSxPQUFBLENBQVFwVixNQUFSLEdBQWlCb1YsT0FBakIsR0FBMkIsSUFBbEMsQ0ExRWlGO0FBQUEsQ0FWckY7QUF1RkEsU0FBU0UsVUFBVCxDQUFvQjlCLElBQXBCLEVBQTBCaUMsT0FBMUIsRUFBbUNULEVBQW5DLEVBQXVDQyxFQUF2QyxFQUEyQy9KLElBQTNDLEVBQWlEO0FBQUEsSUFDN0MsS0FBSyxJQUFJbkwsQ0FBQSxHQUFJLENBQVIsRUFBV0EsQ0FBQSxHQUFJeVQsSUFBQSxDQUFLeFQsTUFBekIsRUFBaUNELENBQUEsSUFBSyxDQUF0QyxFQUF5QztBQUFBLFFBQ3JDLElBQUkyVixDQUFBLEdBQUlsQyxJQUFBLENBQUt6VCxDQUFBLEdBQUltTCxJQUFULENBQVIsQ0FEcUM7QUFBQSxRQUdyQyxJQUFJd0ssQ0FBQSxJQUFLVixFQUFMLElBQVdVLENBQUEsSUFBS1QsRUFBcEIsRUFBd0I7QUFBQSxZQUNwQlEsT0FBQSxDQUFRdFYsSUFBUixDQUFhcVQsSUFBQSxDQUFLelQsQ0FBTCxDQUFiLENBRG9CLENBQUE7QUFBQSxZQUVwQjBWLE9BQUEsQ0FBUXRWLElBQVIsQ0FBYXFULElBQUEsQ0FBS3pULENBQUEsR0FBSSxDQUFULENBQWIsQ0FGb0IsQ0FBQTtBQUFBLFlBR3BCMFYsT0FBQSxDQUFRdFYsSUFBUixDQUFhcVQsSUFBQSxDQUFLelQsQ0FBQSxHQUFJLENBQVQsQ0FBYixDQUhvQixDQUFBO0FBQUEsU0FIYTtBQUFBLEtBREk7QUFBQSxDQXZGakQ7QUFtR0EsU0FBU3dWLFFBQVQsQ0FBa0IvQixJQUFsQixFQUF3QmlDLE9BQXhCLEVBQWlDVCxFQUFqQyxFQUFxQ0MsRUFBckMsRUFBeUMvSixJQUF6QyxFQUErQ3NKLFNBQS9DLEVBQTBEbUIsWUFBMUQsRUFBd0U7QUFBQSxJQUVwRSxJQUFJQyxLQUFBLEdBQVFDLFFBQUEsQ0FBU3JDLElBQVQsQ0FBWixDQUZvRTtBQUFBLElBR3BFLElBQUlzQyxTQUFBLEdBQVk1SyxJQUFBLEtBQVMsQ0FBVCxHQUFhNkssVUFBYixHQUEwQkMsVUFBMUMsQ0FIb0U7QUFBQSxJQUlwRSxJQUFJNVIsR0FBQSxHQUFNb1AsSUFBQSxDQUFLb0IsS0FBZixDQUpvRTtBQUFBLElBS3BFLElBQUlxQixNQUFKLEVBQVkxTCxDQUFaLENBTG9FO0FBQUEsSUFPcEUsS0FBSyxJQUFJeEssQ0FBQSxHQUFJLENBQVIsRUFBV0EsQ0FBQSxHQUFJeVQsSUFBQSxDQUFLeFQsTUFBTCxHQUFjLENBQWxDLEVBQXFDRCxDQUFBLElBQUssQ0FBMUMsRUFBNkM7QUFBQSxRQUN6QyxJQUFJMkwsRUFBQSxHQUFLOEgsSUFBQSxDQUFLelQsQ0FBTCxDQUFULENBRHlDO0FBQUEsUUFFekMsSUFBSTRMLEVBQUEsR0FBSzZILElBQUEsQ0FBS3pULENBQUEsR0FBSSxDQUFULENBQVQsQ0FGeUM7QUFBQSxRQUd6QyxJQUFJbVcsRUFBQSxHQUFLMUMsSUFBQSxDQUFLelQsQ0FBQSxHQUFJLENBQVQsQ0FBVCxDQUh5QztBQUFBLFFBSXpDLElBQUk2TCxFQUFBLEdBQUs0SCxJQUFBLENBQUt6VCxDQUFBLEdBQUksQ0FBVCxDQUFULENBSnlDO0FBQUEsUUFLekMsSUFBSThMLEVBQUEsR0FBSzJILElBQUEsQ0FBS3pULENBQUEsR0FBSSxDQUFULENBQVQsQ0FMeUM7QUFBQSxRQU16QyxJQUFJMlYsQ0FBQSxHQUFJeEssSUFBQSxLQUFTLENBQVQsR0FBYVEsRUFBYixHQUFrQkMsRUFBMUIsQ0FOeUM7QUFBQSxRQU96QyxJQUFJOEYsQ0FBQSxHQUFJdkcsSUFBQSxLQUFTLENBQVQsR0FBYVUsRUFBYixHQUFrQkMsRUFBMUIsQ0FQeUM7QUFBQSxRQVF6QyxJQUFJc0ssTUFBQSxHQUFTLEtBQWIsQ0FSeUM7QUFBQSxRQVV6QyxJQUFJUixZQUFKO0FBQUEsWUFBa0JNLE1BQUEsR0FBUzFSLElBQUEsQ0FBSzRGLElBQUwsQ0FBVTVGLElBQUEsQ0FBS3FMLEdBQUwsQ0FBU2xFLEVBQUEsR0FBS0UsRUFBZCxFQUFrQixDQUFsQixDQUF1QnJILEdBQUFBLElBQUEsQ0FBS3FMLEdBQUwsQ0FBU2pFLEVBQUEsR0FBS0UsRUFBZCxFQUFrQixDQUFsQixDQUFqQyxDQUFULENBVnVCO0FBQUEsUUFZekMsSUFBSTZKLENBQUEsR0FBSVYsRUFBUixFQUFZO0FBQUE7QUFFUixZQUFJdkQsSUFBQUEsQ0FBQSxHQUFJdUQsRUFBUixFQUFZO0FBQUEsZ0JBQ1J6SyxDQUFBLEdBQUl1TCxTQUFBLENBQVVGLEtBQVYsRUFBaUJsSyxFQUFqQixFQUFxQkMsRUFBckIsRUFBeUJDLEVBQXpCLEVBQTZCQyxFQUE3QixFQUFpQ21KLEVBQWpDLENBQUosQ0FEUTtBQUFBLGdCQUVSLElBQUlXLFlBQUo7QUFBQSxvQkFBa0JDLEtBQUEsQ0FBTWhCLEtBQU4sR0FBY3hRLEdBQUEsR0FBTTZSLE1BQUEsR0FBUzFMLENBQTdCLENBRlY7QUFBQSxhQUZKO0FBQUEsU0FBWixNQU1PLElBQUltTCxDQUFBLEdBQUlULEVBQVIsRUFBWTtBQUFBO0FBRWYsWUFBSXhELElBQUFBLENBQUEsR0FBSXdELEVBQVIsRUFBWTtBQUFBLGdCQUNSMUssQ0FBQSxHQUFJdUwsU0FBQSxDQUFVRixLQUFWLEVBQWlCbEssRUFBakIsRUFBcUJDLEVBQXJCLEVBQXlCQyxFQUF6QixFQUE2QkMsRUFBN0IsRUFBaUNvSixFQUFqQyxDQUFKLENBRFE7QUFBQSxnQkFFUixJQUFJVSxZQUFKO0FBQUEsb0JBQWtCQyxLQUFBLENBQU1oQixLQUFOLEdBQWN4USxHQUFBLEdBQU02UixNQUFBLEdBQVMxTCxDQUE3QixDQUZWO0FBQUEsYUFGRztBQUFBLFNBQVosTUFNQTtBQUFBLFlBQ0g2TCxRQUFBLENBQVNSLEtBQVQsRUFBZ0JsSyxFQUFoQixFQUFvQkMsRUFBcEIsRUFBd0J1SyxFQUF4QixDQURHLENBQUE7QUFBQSxTQXhCa0M7QUFBQSxRQTJCekMsSUFBSXpFLENBQUEsR0FBSXVELEVBQUosSUFBVVUsQ0FBQSxJQUFLVixFQUFuQixFQUF1QjtBQUFBO0FBRW5CLFlBQUF6SyxDQUFBLEdBQUl1TCxTQUFBLENBQVVGLEtBQVYsRUFBaUJsSyxFQUFqQixFQUFxQkMsRUFBckIsRUFBeUJDLEVBQXpCLEVBQTZCQyxFQUE3QixFQUFpQ21KLEVBQWpDLENBQUosQ0FGbUI7QUFBQSxZQUduQm1CLE1BQUEsR0FBUyxJQUFULENBSG1CO0FBQUEsU0EzQmtCO0FBQUEsUUFnQ3pDLElBQUkxRSxDQUFBLEdBQUl3RCxFQUFKLElBQVVTLENBQUEsSUFBS1QsRUFBbkIsRUFBdUI7QUFBQTtBQUVuQixZQUFBMUssQ0FBQSxHQUFJdUwsU0FBQSxDQUFVRixLQUFWLEVBQWlCbEssRUFBakIsRUFBcUJDLEVBQXJCLEVBQXlCQyxFQUF6QixFQUE2QkMsRUFBN0IsRUFBaUNvSixFQUFqQyxDQUFKLENBRm1CO0FBQUEsWUFHbkJrQixNQUFBLEdBQVMsSUFBVCxDQUhtQjtBQUFBLFNBaENrQjtBQUFBLFFBc0N6QyxJQUFJLENBQUMzQixTQUFELElBQWMyQixNQUFsQixFQUEwQjtBQUFBLFlBQ3RCLElBQUlSLFlBQUo7QUFBQSxnQkFBa0JDLEtBQUEsQ0FBTWYsR0FBTixHQUFZelEsR0FBQSxHQUFNNlIsTUFBQSxHQUFTMUwsQ0FBM0IsQ0FESTtBQUFBLFlBRXRCa0wsT0FBQSxDQUFRdFYsSUFBUixDQUFheVYsS0FBYixDQUZzQixDQUFBO0FBQUEsWUFHdEJBLEtBQUEsR0FBUUMsUUFBQSxDQUFTckMsSUFBVCxDQUFSLENBSHNCO0FBQUEsU0F0Q2U7QUFBQSxRQTRDekMsSUFBSW1DLFlBQUo7QUFBQSxZQUFrQnZSLEdBQUEsSUFBTzZSLE1BQVAsQ0E1Q3VCO0FBQUEsS0FQdUI7QUFBQTtBQXVEcEUsSUFBQSxJQUFJbEQsSUFBQSxHQUFPUyxJQUFBLENBQUt4VCxNQUFMLEdBQWMsQ0FBekIsQ0F2RG9FO0FBQUEsSUF3RHBFMEwsRUFBQSxHQUFLOEgsSUFBQSxDQUFLVCxJQUFMLENBQUwsQ0F4RG9FO0FBQUEsSUF5RHBFcEgsRUFBQSxHQUFLNkgsSUFBQSxDQUFLVCxJQUFBLEdBQU8sQ0FBWixDQUFMLENBekRvRTtBQUFBLElBMERwRW1ELEVBQUEsR0FBSzFDLElBQUEsQ0FBS1QsSUFBQSxHQUFPLENBQVosQ0FBTCxDQTFEb0U7QUFBQSxJQTJEcEUyQyxDQUFBLEdBQUl4SyxJQUFBLEtBQVMsQ0FBVCxHQUFhUSxFQUFiLEdBQWtCQyxFQUF0QixDQTNEb0U7QUFBQSxJQTREcEUsSUFBSStKLENBQUEsSUFBS1YsRUFBTCxJQUFXVSxDQUFBLElBQUtULEVBQXBCO0FBQUEsUUFBd0JtQixRQUFBLENBQVNSLEtBQVQsRUFBZ0JsSyxFQUFoQixFQUFvQkMsRUFBcEIsRUFBd0J1SyxFQUF4QixDQTVENEMsQ0FBQTtBQUFBO0FBK0RwRSxJQUFBbkQsSUFBQSxHQUFPNkMsS0FBQSxDQUFNNVYsTUFBTixHQUFlLENBQXRCLENBL0RvRTtBQUFBLElBZ0VwRSxJQUFJd1UsU0FBQSxJQUFhekIsSUFBQSxJQUFRLENBQXJCLEtBQTJCNkMsS0FBQSxDQUFNN0MsSUFBTixDQUFBLEtBQWdCNkMsS0FBQSxDQUFNLENBQU4sQ0FBaEIsSUFBNEJBLEtBQUEsQ0FBTTdDLElBQUEsR0FBTyxDQUFiLENBQW9CNkMsS0FBQUEsS0FBQSxDQUFNLENBQU4sQ0FBaEQsQ0FBL0IsRUFBMEY7QUFBQSxRQUN0RlEsUUFBQSxDQUFTUixLQUFULEVBQWdCQSxLQUFBLENBQU0sQ0FBTixDQUFoQixFQUEwQkEsS0FBQSxDQUFNLENBQU4sQ0FBMUIsRUFBb0NBLEtBQUEsQ0FBTSxDQUFOLENBQXBDLENBRHNGLENBQUE7QUFBQSxLQWhFdEI7QUFBQTtBQXFFcEUsSUFBSUEsSUFBQUEsS0FBQSxDQUFNNVYsTUFBVixFQUFrQjtBQUFBLFFBQ2R5VixPQUFBLENBQVF0VixJQUFSLENBQWF5VixLQUFiLENBRGMsQ0FBQTtBQUFBLEtBckVrRDtBQUFBLENBbkd4RTtBQTZLQSxTQUFTQyxRQUFULENBQWtCUSxJQUFsQixFQUF3QjtBQUFBLElBQ3BCLElBQUlULEtBQUEsR0FBUSxFQUFaLENBRG9CO0FBQUEsSUFFcEJBLEtBQUEsQ0FBTWpCLElBQU4sR0FBYTBCLElBQUEsQ0FBSzFCLElBQWxCLENBRm9CO0FBQUEsSUFHcEJpQixLQUFBLENBQU1oQixLQUFOLEdBQWN5QixJQUFBLENBQUt6QixLQUFuQixDQUhvQjtBQUFBLElBSXBCZ0IsS0FBQSxDQUFNZixHQUFOLEdBQVl3QixJQUFBLENBQUt4QixHQUFqQixDQUpvQjtBQUFBLElBS3BCLE9BQU9lLEtBQVAsQ0FMb0I7QUFBQSxDQTdLeEI7QUFxTEEsU0FBU0osU0FBVCxDQUFtQmhDLElBQW5CLEVBQXlCaUMsT0FBekIsRUFBa0NULEVBQWxDLEVBQXNDQyxFQUF0QyxFQUEwQy9KLElBQTFDLEVBQWdEc0osU0FBaEQsRUFBMkQ7QUFBQSxJQUN2RCxLQUFLLElBQUl6VSxDQUFBLEdBQUksQ0FBUixFQUFXQSxDQUFBLEdBQUl5VCxJQUFBLENBQUt4VCxNQUF6QixFQUFpQ0QsQ0FBQSxFQUFqQyxFQUFzQztBQUFBLFFBQ2xDd1YsUUFBQSxDQUFTL0IsSUFBQSxDQUFLelQsQ0FBTCxDQUFULEVBQWtCMFYsT0FBbEIsRUFBMkJULEVBQTNCLEVBQStCQyxFQUEvQixFQUFtQy9KLElBQW5DLEVBQXlDc0osU0FBekMsRUFBb0QsS0FBcEQsQ0FEa0MsQ0FBQTtBQUFBLEtBRGlCO0FBQUEsQ0FyTDNEO0FBMkxBLFNBQVM0QixRQUFULENBQWtCalAsR0FBbEIsRUFBdUJ4QixDQUF2QixFQUEwQkMsQ0FBMUIsRUFBNkJDLENBQTdCLEVBQWdDO0FBQUEsSUFDNUJzQixHQUFBLENBQUloSCxJQUFKLENBQVN3RixDQUFULENBRDRCLENBQUE7QUFBQSxJQUU1QndCLEdBQUEsQ0FBSWhILElBQUosQ0FBU3lGLENBQVQsQ0FGNEIsQ0FBQTtBQUFBLElBRzVCdUIsR0FBQSxDQUFJaEgsSUFBSixDQUFTMEYsQ0FBVCxDQUg0QixDQUFBO0FBQUEsQ0EzTGhDO0FBaU1BLFNBQVNrUSxVQUFULENBQW9CNU8sR0FBcEIsRUFBeUJ1RSxFQUF6QixFQUE2QkMsRUFBN0IsRUFBaUNDLEVBQWpDLEVBQXFDQyxFQUFyQyxFQUF5Q2xHLENBQXpDLEVBQTRDO0FBQUEsSUFDeEMsSUFBSTRFLENBQUEsR0FBSyxDQUFBNUUsQ0FBQSxHQUFJK0YsRUFBSixLQUFXRSxFQUFBLEdBQUtGLEVBQUwsQ0FBcEIsQ0FEd0M7QUFBQSxJQUV4Q3ZFLEdBQUEsQ0FBSWhILElBQUosQ0FBU3dGLENBQVQsQ0FGd0MsQ0FBQTtBQUFBLElBR3hDd0IsR0FBQSxDQUFJaEgsSUFBSixDQUFTd0wsRUFBQSxHQUFNLENBQUFFLEVBQUEsR0FBS0YsRUFBTCxJQUFXcEIsQ0FBMUIsQ0FId0MsQ0FBQTtBQUFBLElBSXhDcEQsR0FBQSxDQUFJaEgsSUFBSixDQUFTLENBQVQsQ0FKd0MsQ0FBQTtBQUFBLElBS3hDLE9BQU9vSyxDQUFQLENBTHdDO0FBQUEsQ0FqTTVDO0FBeU1BLFNBQVN5TCxVQUFULENBQW9CN08sR0FBcEIsRUFBeUJ1RSxFQUF6QixFQUE2QkMsRUFBN0IsRUFBaUNDLEVBQWpDLEVBQXFDQyxFQUFyQyxFQUF5Q2pHLENBQXpDLEVBQTRDO0FBQUEsSUFDeEMsSUFBSTJFLENBQUEsR0FBSyxDQUFBM0UsQ0FBQSxHQUFJK0YsRUFBSixLQUFXRSxFQUFBLEdBQUtGLEVBQUwsQ0FBcEIsQ0FEd0M7QUFBQSxJQUV4Q3hFLEdBQUEsQ0FBSWhILElBQUosQ0FBU3VMLEVBQUEsR0FBTSxDQUFBRSxFQUFBLEdBQUtGLEVBQUwsSUFBV25CLENBQTFCLENBRndDLENBQUE7QUFBQSxJQUd4Q3BELEdBQUEsQ0FBSWhILElBQUosQ0FBU3lGLENBQVQsQ0FId0MsQ0FBQTtBQUFBLElBSXhDdUIsR0FBQSxDQUFJaEgsSUFBSixDQUFTLENBQVQsQ0FKd0MsQ0FBQTtBQUFBLElBS3hDLE9BQU9vSyxDQUFQLENBTHdDO0FBQUE7O0FDck03QixTQUFTK0wsSUFBVCxDQUFjN1MsUUFBZCxFQUF3QjBDLE9BQXhCLEVBQWlDO0FBQUEsSUFDNUMsSUFBSW9RLE1BQUEsR0FBU3BRLE9BQUEsQ0FBUW9RLE1BQVIsR0FBaUJwUSxPQUFBLENBQVFsQixNQUF0QyxDQUQ0QztBQUFBLElBRTVDLElBQUl1UixNQUFBLEdBQVMvUyxRQUFiLENBRjRDO0FBQUEsSUFHNUMsSUFBSWdHLElBQUEsR0FBUXFMLElBQUEsQ0FBS3JSLFFBQUwsRUFBZSxDQUFmLEVBQWtCLENBQUMsQ0FBRCxHQUFLOFMsTUFBdkIsRUFBK0JBLE1BQS9CLEVBQTJDLENBQTNDLEVBQThDLENBQUMsQ0FBL0MsRUFBa0QsQ0FBbEQsRUFBcURwUSxPQUFyRCxDQUFaLENBSDRDO0FBQUE7QUFJNUMsSUFBSXVELElBQUFBLEtBQUEsR0FBUW9MLElBQUEsQ0FBS3JSLFFBQUwsRUFBZSxDQUFmLEVBQW1CLENBQUEsR0FBSThTLE1BQXZCLEVBQStCLElBQUlBLE1BQW5DLEVBQTJDLENBQTNDLEVBQThDLENBQUMsQ0FBL0MsRUFBa0QsQ0FBbEQsRUFBcURwUSxPQUFyRCxDQUFaLENBSjRDO0FBQUE7QUFNNUMsSUFBSXNELElBQUFBLElBQUEsSUFBUUMsS0FBWixFQUFtQjtBQUFBLFFBQ2Y4TSxNQUFBLEdBQVMxQixJQUFBLENBQUtyUixRQUFMLEVBQWUsQ0FBZixFQUFrQixDQUFDOFMsTUFBbkIsRUFBMkIsQ0FBQSxHQUFJQSxNQUEvQixFQUF1QyxDQUF2QyxFQUEwQyxDQUFDLENBQTNDLEVBQThDLENBQTlDLEVBQWlEcFEsT0FBakQsQ0FBNkQsSUFBQSxFQUF0RSxDQURlO0FBQUE7QUFHZixRQUFBLElBQUlzRCxJQUFKO0FBQUEsWUFBVStNLE1BQUEsR0FBU0Msa0JBQUEsQ0FBbUJoTixJQUFuQixFQUF5QixDQUF6QixDQUFBLENBQTRCaUYsTUFBNUIsQ0FBbUM4SCxNQUFuQyxDQUFULENBSEs7QUFBQTtBQUlmLFFBQUEsSUFBSTlNLEtBQUo7QUFBQSxZQUFXOE0sTUFBQSxHQUFTQSxNQUFBLENBQU85SCxNQUFQLENBQWMrSCxrQkFBQSxDQUFtQi9NLEtBQW5CLEVBQTBCLENBQUMsQ0FBM0IsQ0FBZCxDQUFULENBQUE7QUFKSSxLQU55QjtBQUFBLElBYTVDLE9BQU84TSxNQUFQLENBYjRDO0FBQUEsQ0FKaEQ7QUFvQkEsU0FBU0Msa0JBQVQsQ0FBNEJoVCxRQUE1QixFQUFzQ3dNLE1BQXRDLEVBQThDO0FBQUEsSUFDMUMsSUFBSXlHLFdBQUEsR0FBYyxFQUFsQixDQUQwQztBQUFBLElBRzFDLEtBQUssSUFBSTNXLENBQUEsR0FBSSxDQUFSLEVBQVdBLENBQUEsR0FBSTBELFFBQUEsQ0FBU3pELE1BQTdCLEVBQXFDRCxDQUFBLEVBQXJDLEVBQTBDO0FBQUEsUUFDdEMsSUFBSWdGLE9BQUEsR0FBVXRCLFFBQUEsQ0FBUzFELENBQVQsQ0FBZCxFQUNJeUQsSUFBQSxHQUFPdUIsT0FBQSxDQUFRdkIsSUFEbkIsQ0FEc0M7QUFBQSxRQUl0QyxJQUFJNlIsV0FBSixDQUpzQztBQUFBLFFBTXRDLElBQUk3UixJQUFBLEtBQVMsT0FBVCxJQUFvQkEsSUFBQSxLQUFTLFlBQTdCLElBQTZDQSxJQUFBLEtBQVMsWUFBMUQsRUFBd0U7QUFBQSxZQUNwRTZSLFdBQUEsR0FBY3NCLFdBQUEsQ0FBWTVSLE9BQUEsQ0FBUXBCLFFBQXBCLEVBQThCc00sTUFBOUIsQ0FBZCxDQURvRTtBQUFBLFNBQXhFLE1BR08sSUFBSXpNLElBQUEsS0FBUyxpQkFBVCxJQUE4QkEsSUFBQSxLQUFTLFNBQTNDLEVBQXNEO0FBQUEsWUFDekQ2UixXQUFBLEdBQWMsRUFBZCxDQUR5RDtBQUFBLFlBRXpELEtBQUssSUFBSWhSLENBQUEsR0FBSSxDQUFSLEVBQVdBLENBQUEsR0FBSVUsT0FBQSxDQUFRcEIsUUFBUixDQUFpQjNELE1BQXJDLEVBQTZDcUUsQ0FBQSxFQUE3QyxFQUFrRDtBQUFBLGdCQUM5Q2dSLFdBQUEsQ0FBWWxWLElBQVosQ0FBaUJ3VyxXQUFBLENBQVk1UixPQUFBLENBQVFwQixRQUFSLENBQWlCVSxDQUFqQixDQUFaLEVBQWlDNEwsTUFBakMsQ0FBakIsQ0FEOEMsQ0FBQTtBQUFBLGFBRk87QUFBQSxTQUF0RCxNQUtBLElBQUl6TSxJQUFBLEtBQVMsY0FBYixFQUE2QjtBQUFBLFlBQ2hDNlIsV0FBQSxHQUFjLEVBQWQsQ0FEZ0M7QUFBQSxZQUVoQyxLQUFLaFIsQ0FBQSxHQUFJLENBQVQsRUFBWUEsQ0FBQSxHQUFJVSxPQUFBLENBQVFwQixRQUFSLENBQWlCM0QsTUFBakMsRUFBeUNxRSxDQUFBLEVBQXpDLEVBQThDO0FBQUEsZ0JBQzFDLElBQUl1UyxVQUFBLEdBQWEsRUFBakIsQ0FEMEM7QUFBQSxnQkFFMUMsS0FBSyxJQUFJblgsQ0FBQSxHQUFJLENBQVIsRUFBV0EsQ0FBQSxHQUFJc0YsT0FBQSxDQUFRcEIsUUFBUixDQUFpQlUsQ0FBakIsQ0FBQSxDQUFvQnJFLE1BQXhDLEVBQWdEUCxDQUFBLEVBQWhELEVBQXFEO0FBQUEsb0JBQ2pEbVgsVUFBQSxDQUFXelcsSUFBWCxDQUFnQndXLFdBQUEsQ0FBWTVSLE9BQUEsQ0FBUXBCLFFBQVIsQ0FBaUJVLENBQWpCLENBQW9CNUUsQ0FBQUEsQ0FBcEIsQ0FBWixFQUFvQ3dRLE1BQXBDLENBQWhCLENBRGlELENBQUE7QUFBQSxpQkFGWDtBQUFBLGdCQUsxQ29GLFdBQUEsQ0FBWWxWLElBQVosQ0FBaUJ5VyxVQUFqQixDQUwwQyxDQUFBO0FBQUEsYUFGZDtBQUFBLFNBZEU7QUFBQSxRQXlCdENGLFdBQUEsQ0FBWXZXLElBQVosQ0FBaUJvVCxhQUFBLENBQWN4TyxPQUFBLENBQVE5RSxFQUF0QixFQUEwQnVELElBQTFCLEVBQWdDNlIsV0FBaEMsRUFBNkN0USxPQUFBLENBQVFLLElBQXJELENBQWpCLENBekJzQyxDQUFBO0FBQUEsS0FIQTtBQUFBLElBK0IxQyxPQUFPc1IsV0FBUCxDQS9CMEM7QUFBQSxDQXBCOUM7QUFzREEsU0FBU0MsV0FBVCxDQUFxQnpLLE1BQXJCLEVBQTZCK0QsTUFBN0IsRUFBcUM7QUFBQSxJQUNqQyxJQUFJNEcsU0FBQSxHQUFZLEVBQWhCLENBRGlDO0FBQUEsSUFFakNBLFNBQUEsQ0FBVWxDLElBQVYsR0FBaUJ6SSxNQUFBLENBQU95SSxJQUF4QixDQUZpQztBQUFBLElBSWpDLElBQUl6SSxNQUFBLENBQU8wSSxLQUFQLEtBQWlCL1YsU0FBckIsRUFBZ0M7QUFBQSxRQUM1QmdZLFNBQUEsQ0FBVWpDLEtBQVYsR0FBa0IxSSxNQUFBLENBQU8wSSxLQUF6QixDQUQ0QjtBQUFBLFFBRTVCaUMsU0FBQSxDQUFVaEMsR0FBVixHQUFnQjNJLE1BQUEsQ0FBTzJJLEdBQXZCLENBRjRCO0FBQUEsS0FKQztBQUFBLElBU2pDLEtBQUssSUFBSTlVLENBQUEsR0FBSSxDQUFSLEVBQVdBLENBQUEsR0FBSW1NLE1BQUEsQ0FBT2xNLE1BQTNCLEVBQW1DRCxDQUFBLElBQUssQ0FBeEMsRUFBMkM7QUFBQSxRQUN2QzhXLFNBQUEsQ0FBVTFXLElBQVYsQ0FBZStMLE1BQUEsQ0FBT25NLENBQVAsQ0FBWWtRLEdBQUFBLE1BQTNCLEVBQW1DL0QsTUFBQSxDQUFPbk0sQ0FBQSxHQUFJLENBQVgsQ0FBbkMsRUFBa0RtTSxNQUFBLENBQU9uTSxDQUFBLEdBQUksQ0FBWCxDQUFsRCxDQUR1QyxDQUFBO0FBQUEsS0FUVjtBQUFBLElBWWpDLE9BQU84VyxTQUFQLENBWmlDO0FBQUE7O0FDbkRyQztBQUFBO0FBQWUsU0FBU0MsYUFBVCxDQUF1QjVQLElBQXZCLEVBQTZCakMsTUFBN0IsRUFBcUM7QUFBQSxJQUNoRCxJQUFJaUMsSUFBQSxDQUFLNlAsV0FBVDtBQUFBLFFBQXNCLE9BQU83UCxJQUFQLENBRDBCO0FBQUEsSUFHaEQsSUFBSW1KLEVBQUEsR0FBSyxDQUFLbkosSUFBQUEsSUFBQSxDQUFLckIsQ0FBbkIsRUFDSW1SLEVBQUEsR0FBSzlQLElBQUEsQ0FBS3ZCLENBRGQsRUFFSXNSLEVBQUEsR0FBSy9QLElBQUEsQ0FBS3RCLENBRmQsRUFHSTdGLENBSEosRUFHT3NFLENBSFAsRUFHVTVFLENBSFYsQ0FIZ0Q7QUFBQSxJQVFoRCxLQUFLTSxDQUFBLEdBQUksQ0FBVCxFQUFZQSxDQUFBLEdBQUltSCxJQUFBLENBQUt6RCxRQUFMLENBQWN6RCxNQUE5QixFQUFzQ0QsQ0FBQSxFQUF0QyxFQUEyQztBQUFBLFFBQ3ZDLElBQUlnRixPQUFBLEdBQVVtQyxJQUFBLENBQUt6RCxRQUFMLENBQWMxRCxDQUFkLENBQWQsRUFDSXlULElBQUEsR0FBT3pPLE9BQUEsQ0FBUXBCLFFBRG5CLEVBRUlILElBQUEsR0FBT3VCLE9BQUEsQ0FBUXZCLElBRm5CLENBRHVDO0FBQUEsUUFLdkN1QixPQUFBLENBQVFwQixRQUFSLEdBQW1CLEVBQW5CLENBTHVDO0FBQUEsUUFPdkMsSUFBSUgsSUFBQSxLQUFTLENBQWIsRUFBZ0I7QUFBQSxZQUNaLEtBQUthLENBQUEsR0FBSSxDQUFULEVBQVlBLENBQUEsR0FBSW1QLElBQUEsQ0FBS3hULE1BQXJCLEVBQTZCcUUsQ0FBQSxJQUFLLENBQWxDLEVBQXFDO0FBQUEsZ0JBQ2pDVSxPQUFBLENBQVFwQixRQUFSLENBQWlCeEQsSUFBakIsQ0FBc0IrVyxjQUFBLENBQWUxRCxJQUFBLENBQUtuUCxDQUFMLENBQWYsRUFBd0JtUCxJQUFBLENBQUtuUCxDQUFBLEdBQUksQ0FBVCxDQUF4QixFQUFxQ1ksTUFBckMsRUFBNkNvTCxFQUE3QyxFQUFpRDJHLEVBQWpELEVBQXFEQyxFQUFyRCxDQUF0QixDQURpQyxDQUFBO0FBQUEsYUFEekI7QUFBQSxTQUFoQixNQUlPO0FBQUEsWUFDSCxLQUFLNVMsQ0FBQSxHQUFJLENBQVQsRUFBWUEsQ0FBQSxHQUFJbVAsSUFBQSxDQUFLeFQsTUFBckIsRUFBNkJxRSxDQUFBLEVBQTdCLEVBQWtDO0FBQUEsZ0JBQzlCLElBQUlMLElBQUEsR0FBTyxFQUFYLENBRDhCO0FBQUEsZ0JBRTlCLEtBQUt2RSxDQUFBLEdBQUksQ0FBVCxFQUFZQSxDQUFBLEdBQUkrVCxJQUFBLENBQUtuUCxDQUFMLEVBQVFyRSxNQUF4QixFQUFnQ1AsQ0FBQSxJQUFLLENBQXJDLEVBQXdDO0FBQUEsb0JBQ3BDdUUsSUFBQSxDQUFLN0QsSUFBTCxDQUFVK1csY0FBQSxDQUFlMUQsSUFBQSxDQUFLblAsQ0FBTCxDQUFBLENBQVE1RSxDQUFSLENBQWYsRUFBMkIrVCxJQUFBLENBQUtuUCxDQUFMLENBQUEsQ0FBUTVFLENBQUEsR0FBSSxDQUFaLENBQTNCLEVBQTJDd0YsTUFBM0MsRUFBbURvTCxFQUFuRCxFQUF1RDJHLEVBQXZELEVBQTJEQyxFQUEzRCxDQUFWLENBRG9DLENBQUE7QUFBQSxpQkFGVjtBQUFBLGdCQUs5QmxTLE9BQUEsQ0FBUXBCLFFBQVIsQ0FBaUJ4RCxJQUFqQixDQUFzQjZELElBQXRCLENBTDhCLENBQUE7QUFBQSxhQUQvQjtBQUFBLFNBWGdDO0FBQUEsS0FSSztBQUFBLElBOEJoRGtELElBQUEsQ0FBSzZQLFdBQUwsR0FBbUIsSUFBbkIsQ0E5QmdEO0FBQUEsSUFnQ2hELE9BQU83UCxJQUFQLENBaENnRDtBQUFBLENBSHBEO0FBc0NBLFNBQVNnUSxjQUFULENBQXdCdlIsQ0FBeEIsRUFBMkJDLENBQTNCLEVBQThCWCxNQUE5QixFQUFzQ29MLEVBQXRDLEVBQTBDMkcsRUFBMUMsRUFBOENDLEVBQTlDLEVBQWtEO0FBQUEsSUFDOUMsT0FBTztBQUFBLFFBQ0gxUyxJQUFBLENBQUs4TSxLQUFMLENBQVdwTSxNQUFBLElBQVVVLENBQUEsR0FBSTBLLEVBQUosR0FBUzJHLEVBQVQsQ0FBckIsQ0FERztBQUFBLFFBRUh6UyxJQUFBLENBQUs4TSxLQUFMLENBQVdwTSxNQUFBLElBQVVXLENBQUEsR0FBSXlLLEVBQUosR0FBUzRHLEVBQVQsQ0FBckIsQ0FGRztBQUFBLEtBQVAsQ0FEOEM7QUFBQTs7QUNyQ25DLFNBQVNFLFVBQVQsQ0FBb0IxVCxRQUFwQixFQUE4Qm9DLENBQTlCLEVBQWlDbVIsRUFBakMsRUFBcUNDLEVBQXJDLEVBQXlDOVEsT0FBekMsRUFBa0Q7QUFBQSxJQUM3RCxJQUFJNE4sU0FBQSxHQUFZbE8sQ0FBQSxLQUFNTSxPQUFBLENBQVF5RyxPQUFkLEdBQXdCLENBQXhCLEdBQTRCekcsT0FBQSxDQUFRNE4sU0FBUixJQUFzQixDQUFBLENBQUEsSUFBS2xPLENBQUwsSUFBVU0sT0FBQSxDQUFRbEIsTUFBbkIsQ0FBakUsQ0FENkQ7QUFBQSxJQUU3RCxJQUFJaUMsSUFBQSxHQUFPO0FBQUEsUUFDUHpELFFBQUEsRUFBVSxFQURIO0FBQUEsUUFFUHVMLFNBQUEsRUFBVyxDQUZKO0FBQUEsUUFHUG9JLGFBQUEsRUFBZSxDQUhSO0FBQUEsUUFJUEMsV0FBQSxFQUFhLENBSk47QUFBQSxRQUtQL1YsTUFBQSxFQUFRLElBTEQ7QUFBQSxRQU1QcUUsQ0FBQSxFQUFHcVIsRUFOSTtBQUFBLFFBT1BwUixDQUFBLEVBQUdxUixFQVBJO0FBQUEsUUFRUHBSLENBQUEsRUFBR0EsQ0FSSTtBQUFBLFFBU1BrUixXQUFBLEVBQWEsS0FUTjtBQUFBLFFBVVBsTSxJQUFBLEVBQU0sQ0FWQztBQUFBLFFBV1BDLElBQUEsRUFBTSxDQVhDO0FBQUEsUUFZUEMsSUFBQSxFQUFNLENBQUMsQ0FaQTtBQUFBLFFBYVBDLElBQUEsRUFBTSxDQWJDO0FBQUEsS0FBWCxDQUY2RDtBQUFBLElBaUI3RCxLQUFLLElBQUlqTCxDQUFBLEdBQUksQ0FBUixFQUFXQSxDQUFBLEdBQUkwRCxRQUFBLENBQVN6RCxNQUE3QixFQUFxQ0QsQ0FBQSxFQUFyQyxFQUEwQztBQUFBLFFBQ3RDbUgsSUFBQSxDQUFLbVEsV0FBTCxFQURzQyxDQUFBO0FBQUEsUUFFdENDLFVBQUEsQ0FBV3BRLElBQVgsRUFBaUJ6RCxRQUFBLENBQVMxRCxDQUFULENBQWpCLEVBQThCZ1UsU0FBOUIsRUFBeUM1TixPQUF6QyxDQUZzQyxDQUFBO0FBQUEsUUFJdEMsSUFBSTBFLElBQUEsR0FBT3BILFFBQUEsQ0FBUzFELENBQVQsQ0FBQSxDQUFZOEssSUFBdkIsQ0FKc0M7QUFBQSxRQUt0QyxJQUFJQyxJQUFBLEdBQU9ySCxRQUFBLENBQVMxRCxDQUFULENBQUEsQ0FBWStLLElBQXZCLENBTHNDO0FBQUEsUUFNdEMsSUFBSUMsSUFBQSxHQUFPdEgsUUFBQSxDQUFTMUQsQ0FBVCxDQUFBLENBQVlnTCxJQUF2QixDQU5zQztBQUFBLFFBT3RDLElBQUlDLElBQUEsR0FBT3ZILFFBQUEsQ0FBUzFELENBQVQsQ0FBQSxDQUFZaUwsSUFBdkIsQ0FQc0M7QUFBQSxRQVN0QyxJQUFJSCxJQUFBLEdBQU8zRCxJQUFBLENBQUsyRCxJQUFoQjtBQUFBLFlBQXNCM0QsSUFBQSxDQUFLMkQsSUFBTCxHQUFZQSxJQUFaLENBVGdCO0FBQUEsUUFVdEMsSUFBSUMsSUFBQSxHQUFPNUQsSUFBQSxDQUFLNEQsSUFBaEI7QUFBQSxZQUFzQjVELElBQUEsQ0FBSzRELElBQUwsR0FBWUEsSUFBWixDQVZnQjtBQUFBLFFBV3RDLElBQUlDLElBQUEsR0FBTzdELElBQUEsQ0FBSzZELElBQWhCO0FBQUEsWUFBc0I3RCxJQUFBLENBQUs2RCxJQUFMLEdBQVlBLElBQVosQ0FYZ0I7QUFBQSxRQVl0QyxJQUFJQyxJQUFBLEdBQU85RCxJQUFBLENBQUs4RCxJQUFoQjtBQUFBLFlBQXNCOUQsSUFBQSxDQUFLOEQsSUFBTCxHQUFZQSxJQUFaLENBWmdCO0FBQUEsS0FqQm1CO0FBQUEsSUErQjdELE9BQU85RCxJQUFQLENBL0I2RDtBQUFBLENBRGpFO0FBbUNBLFNBQVNvUSxVQUFULENBQW9CcFEsSUFBcEIsRUFBMEJuQyxPQUExQixFQUFtQ2dQLFNBQW5DLEVBQThDNU4sT0FBOUMsRUFBdUQ7QUFBQSxJQUVuRCxJQUFJcU4sSUFBQSxHQUFPek8sT0FBQSxDQUFRcEIsUUFBbkIsRUFDSUgsSUFBQSxHQUFPdUIsT0FBQSxDQUFRdkIsSUFEbkIsRUFFSStULFVBQUEsR0FBYSxFQUZqQixDQUZtRDtBQUFBLElBTW5ELElBQUkvVCxJQUFBLEtBQVMsT0FBVCxJQUFvQkEsSUFBQSxLQUFTLFlBQWpDLEVBQStDO0FBQUEsUUFDM0MsS0FBSyxJQUFJekQsQ0FBQSxHQUFJLENBQVIsRUFBV0EsQ0FBQSxHQUFJeVQsSUFBQSxDQUFLeFQsTUFBekIsRUFBaUNELENBQUEsSUFBSyxDQUF0QyxFQUF5QztBQUFBLFlBQ3JDd1gsVUFBQSxDQUFXcFgsSUFBWCxDQUFnQnFULElBQUEsQ0FBS3pULENBQUwsQ0FBaEIsQ0FEcUMsQ0FBQTtBQUFBLFlBRXJDd1gsVUFBQSxDQUFXcFgsSUFBWCxDQUFnQnFULElBQUEsQ0FBS3pULENBQUEsR0FBSSxDQUFULENBQWhCLENBRnFDLENBQUE7QUFBQSxZQUdyQ21ILElBQUEsQ0FBSzhILFNBQUwsRUFIcUMsQ0FBQTtBQUFBLFlBSXJDOUgsSUFBQSxDQUFLa1EsYUFBTCxFQUpxQyxDQUFBO0FBQUEsU0FERTtBQUFBLEtBQS9DLE1BUU8sSUFBSTVULElBQUEsS0FBUyxZQUFiLEVBQTJCO0FBQUEsUUFDOUJnVSxPQUFBLENBQVFELFVBQVIsRUFBb0IvRCxJQUFwQixFQUEwQnRNLElBQTFCLEVBQWdDNk0sU0FBaEMsRUFBMkMsS0FBM0MsRUFBa0QsS0FBbEQsQ0FEOEIsQ0FBQTtBQUFBLEtBQTNCLE1BR0EsSUFBSXZRLElBQUEsS0FBUyxpQkFBVCxJQUE4QkEsSUFBQSxLQUFTLFNBQTNDLEVBQXNEO0FBQUEsUUFDekQsS0FBS3pELENBQUEsR0FBSSxDQUFULEVBQVlBLENBQUEsR0FBSXlULElBQUEsQ0FBS3hULE1BQXJCLEVBQTZCRCxDQUFBLEVBQTdCLEVBQWtDO0FBQUEsWUFDOUJ5WCxPQUFBLENBQVFELFVBQVIsRUFBb0IvRCxJQUFBLENBQUt6VCxDQUFMLENBQXBCLEVBQTZCbUgsSUFBN0IsRUFBbUM2TSxTQUFuQyxFQUE4Q3ZRLElBQUEsS0FBUyxTQUF2RCxFQUFrRXpELENBQUEsS0FBTSxDQUF4RSxDQUQ4QixDQUFBO0FBQUEsU0FEdUI7QUFBQSxLQUF0RCxNQUtBLElBQUl5RCxJQUFBLEtBQVMsY0FBYixFQUE2QjtBQUFBLFFBRWhDLEtBQUssSUFBSS9ELENBQUEsR0FBSSxDQUFSLEVBQVdBLENBQUEsR0FBSStULElBQUEsQ0FBS3hULE1BQXpCLEVBQWlDUCxDQUFBLEVBQWpDLEVBQXNDO0FBQUEsWUFDbEMsSUFBSTRVLE9BQUEsR0FBVWIsSUFBQSxDQUFLL1QsQ0FBTCxDQUFkLENBRGtDO0FBQUEsWUFFbEMsS0FBS00sQ0FBQSxHQUFJLENBQVQsRUFBWUEsQ0FBQSxHQUFJc1UsT0FBQSxDQUFRclUsTUFBeEIsRUFBZ0NELENBQUEsRUFBaEMsRUFBcUM7QUFBQSxnQkFDakN5WCxPQUFBLENBQVFELFVBQVIsRUFBb0JsRCxPQUFBLENBQVF0VSxDQUFSLENBQXBCLEVBQWdDbUgsSUFBaEMsRUFBc0M2TSxTQUF0QyxFQUFpRCxJQUFqRCxFQUF1RGhVLENBQUEsS0FBTSxDQUE3RCxDQURpQyxDQUFBO0FBQUEsYUFGSDtBQUFBLFNBRk47QUFBQSxLQXRCZTtBQUFBLElBZ0NuRCxJQUFJd1gsVUFBQSxDQUFXdlgsTUFBZixFQUF1QjtBQUFBLFFBQ25CLElBQUlvRixJQUFBLEdBQU9MLE9BQUEsQ0FBUUssSUFBUixJQUFnQixJQUEzQixDQURtQjtBQUFBLFFBRW5CLElBQUk1QixJQUFBLEtBQVMsWUFBVCxJQUF5QjJDLE9BQUEsQ0FBUWdPLFdBQXJDLEVBQWtEO0FBQUEsWUFDOUMvTyxJQUFBLEdBQU8sRUFBUCxDQUQ4QztBQUFBLFlBRTlDLEtBQVNqRyxJQUFBQSxHQUFULElBQWdCNEYsT0FBQSxDQUFRSyxJQUF4QjtBQUFBLGdCQUE4QkEsSUFBQSxDQUFLakcsR0FBTCxDQUFZNEYsR0FBQUEsT0FBQSxDQUFRSyxJQUFSLENBQWFqRyxHQUFiLENBQVosQ0FGZ0I7QUFBQSxZQUc5Q2lHLElBQUEsQ0FBSyxtQkFBTCxDQUE0Qm9PLEdBQUFBLElBQUEsQ0FBS29CLEtBQUwsR0FBYXBCLElBQUEsQ0FBS21CLElBQTlDLENBSDhDO0FBQUEsWUFJOUN2UCxJQUFBLENBQUssaUJBQUwsQ0FBMEJvTyxHQUFBQSxJQUFBLENBQUtxQixHQUFMLEdBQVdyQixJQUFBLENBQUttQixJQUExQyxDQUo4QztBQUFBLFNBRi9CO0FBQUEsUUFRbkIsSUFBSThDLFdBQUEsR0FBYztBQUFBLFlBQ2Q5VCxRQUFBLEVBQVU0VCxVQURJO0FBQUEsWUFFZC9ULElBQUEsRUFBTUEsSUFBQSxLQUFTLFNBQVQsSUFBc0JBLElBQUEsS0FBUyxjQUEvQixHQUFnRCxDQUFoRCxHQUNGQSxJQUFBLEtBQVMsWUFBVCxJQUF5QkEsSUFBQSxLQUFTLGlCQUFsQyxHQUFzRCxDQUF0RCxHQUEwRCxDQUhoRDtBQUFBLFlBSWQ0QixJQUFBLEVBQU1BLElBSlE7QUFBQSxTQUFsQixDQVJtQjtBQUFBLFFBY25CLElBQUlMLE9BQUEsQ0FBUTlFLEVBQVIsS0FBZSxJQUFuQixFQUF5QjtBQUFBLFlBQ3JCd1gsV0FBQSxDQUFZeFgsRUFBWixHQUFpQjhFLE9BQUEsQ0FBUTlFLEVBQXpCLENBRHFCO0FBQUEsU0FkTjtBQUFBLFFBaUJuQmlILElBQUEsQ0FBS3pELFFBQUwsQ0FBY3RELElBQWQsQ0FBbUJzWCxXQUFuQixDQWpCbUIsQ0FBQTtBQUFBLEtBaEM0QjtBQUFBLENBbkN2RDtBQXdGQSxTQUFTRCxPQUFULENBQWlCcFgsTUFBakIsRUFBeUJvVCxJQUF6QixFQUErQnRNLElBQS9CLEVBQXFDNk0sU0FBckMsRUFBZ0RTLFNBQWhELEVBQTJEa0QsT0FBM0QsRUFBb0U7QUFBQSxJQUNoRSxJQUFJMUUsV0FBQSxHQUFjZSxTQUFBLEdBQVlBLFNBQTlCLENBRGdFO0FBQUEsSUFHaEUsSUFBSUEsU0FBQSxHQUFZLENBQVosSUFBa0JQLElBQUEsQ0FBS21CLElBQUwsSUFBYUgsU0FBQSxHQUFZeEIsV0FBWixHQUEwQmUsU0FBMUIsQ0FBbkMsRUFBMEU7QUFBQSxRQUN0RTdNLElBQUEsQ0FBSzhILFNBQUwsSUFBa0J3RSxJQUFBLENBQUt4VCxNQUFMLEdBQWMsQ0FBaEMsQ0FEc0U7QUFBQSxRQUV0RSxPQUZzRTtBQUFBLEtBSFY7QUFBQSxJQVFoRSxJQUFJZ0UsSUFBQSxHQUFPLEVBQVgsQ0FSZ0U7QUFBQSxJQVVoRSxLQUFLLElBQUlqRSxDQUFBLEdBQUksQ0FBUixFQUFXQSxDQUFBLEdBQUl5VCxJQUFBLENBQUt4VCxNQUF6QixFQUFpQ0QsQ0FBQSxJQUFLLENBQXRDLEVBQXlDO0FBQUEsUUFDckMsSUFBSWdVLFNBQUEsS0FBYyxDQUFkLElBQW1CUCxJQUFBLENBQUt6VCxDQUFBLEdBQUksQ0FBVCxDQUFjaVQsR0FBQUEsV0FBckMsRUFBa0Q7QUFBQSxZQUM5QzlMLElBQUEsQ0FBS2tRLGFBQUwsRUFEOEMsQ0FBQTtBQUFBLFlBRTlDcFQsSUFBQSxDQUFLN0QsSUFBTCxDQUFVcVQsSUFBQSxDQUFLelQsQ0FBTCxDQUFWLENBRjhDLENBQUE7QUFBQSxZQUc5Q2lFLElBQUEsQ0FBSzdELElBQUwsQ0FBVXFULElBQUEsQ0FBS3pULENBQUEsR0FBSSxDQUFULENBQVYsQ0FIOEMsQ0FBQTtBQUFBLFNBRGI7QUFBQSxRQU1yQ21ILElBQUEsQ0FBSzhILFNBQUwsRUFOcUMsQ0FBQTtBQUFBLEtBVnVCO0FBQUEsSUFtQmhFLElBQUl3RixTQUFKO0FBQUEsUUFBZW5SLE1BQUEsQ0FBT1csSUFBUCxFQUFhMFQsT0FBYixDQW5CaUQsQ0FBQTtBQUFBLElBcUJoRXRYLE1BQUEsQ0FBT0QsSUFBUCxDQUFZNkQsSUFBWixDQXJCZ0UsQ0FBQTtBQUFBLENBeEZwRTtBQWdIQSxTQUFTWCxNQUFULENBQWdCVyxJQUFoQixFQUFzQjJULFNBQXRCLEVBQWlDO0FBQUEsSUFDN0IsSUFBSXpULElBQUEsR0FBTyxDQUFYLENBRDZCO0FBQUEsSUFFN0IsS0FBSyxJQUFJbkUsQ0FBQSxHQUFJLENBQVIsRUFBV3FFLEdBQUEsR0FBTUosSUFBQSxDQUFLaEUsTUFBdEIsRUFBOEJxRSxDQUFBLEdBQUlELEdBQUEsR0FBTSxDQUF4QyxFQUEyQ3JFLENBQUEsR0FBSXFFLEdBQXBELEVBQXlEQyxDQUFBLEdBQUl0RSxDQUFKLEVBQU9BLENBQUEsSUFBSyxDQUFyRSxFQUF3RTtBQUFBLFFBQ3BFbUUsSUFBQSxJQUFTLENBQUFGLElBQUEsQ0FBS2pFLENBQUwsQ0FBQSxHQUFVaUUsSUFBQSxDQUFLSyxDQUFMLENBQVYsS0FBc0JMLElBQUEsQ0FBS2pFLENBQUEsR0FBSSxDQUFULENBQWNpRSxHQUFBQSxJQUFBLENBQUtLLENBQUEsR0FBSSxDQUFULENBQWQsQ0FBL0IsQ0FEb0U7QUFBQSxLQUYzQztBQUFBLElBSzdCLElBQUlILElBQUEsR0FBTyxDQUFQLEtBQWF5VCxTQUFqQixFQUE0QjtBQUFBLFFBQ3hCLEtBQUs1WCxDQUFBLEdBQUksQ0FBSixFQUFPcUUsR0FBQSxHQUFNSixJQUFBLENBQUtoRSxNQUF2QixFQUErQkQsQ0FBQSxHQUFJcUUsR0FBQSxHQUFNLENBQXpDLEVBQTRDckUsQ0FBQSxJQUFLLENBQWpELEVBQW9EO0FBQUEsWUFDaEQsSUFBSTRGLENBQUEsR0FBSTNCLElBQUEsQ0FBS2pFLENBQUwsQ0FBUixDQURnRDtBQUFBLFlBRWhELElBQUk2RixDQUFBLEdBQUk1QixJQUFBLENBQUtqRSxDQUFBLEdBQUksQ0FBVCxDQUFSLENBRmdEO0FBQUEsWUFHaERpRSxJQUFBLENBQUtqRSxDQUFMLENBQUEsR0FBVWlFLElBQUEsQ0FBS0ksR0FBQSxHQUFNLENBQU4sR0FBVXJFLENBQWYsQ0FBVixDQUhnRDtBQUFBLFlBSWhEaUUsSUFBQSxDQUFLakUsQ0FBQSxHQUFJLENBQVQsQ0FBY2lFLEdBQUFBLElBQUEsQ0FBS0ksR0FBQSxHQUFNLENBQU4sR0FBVXJFLENBQWYsQ0FBZCxDQUpnRDtBQUFBLFlBS2hEaUUsSUFBQSxDQUFLSSxHQUFBLEdBQU0sQ0FBTixHQUFVckUsQ0FBZixDQUFvQjRGLEdBQUFBLENBQXBCLENBTGdEO0FBQUEsWUFNaEQzQixJQUFBLENBQUtJLEdBQUEsR0FBTSxDQUFOLEdBQVVyRSxDQUFmLENBQW9CNkYsR0FBQUEsQ0FBcEIsQ0FOZ0Q7QUFBQSxTQUQ1QjtBQUFBLEtBTEM7QUFBQTs7QUN6R2pDO0FBQWUsU0FBU2dTLFNBQVQsQ0FBbUJoRSxJQUFuQixFQUF5QnpOLE9BQXpCLEVBQWtDO0FBQUEsSUFDN0MsT0FBTyxJQUFJMFIsU0FBSixDQUFjakUsSUFBZCxFQUFvQnpOLE9BQXBCLENBQVAsQ0FENkM7QUFBQSxDQVBqRDtBQVdBLFNBQVMwUixTQUFULENBQW1CakUsSUFBbkIsRUFBeUJ6TixPQUF6QixFQUFrQztBQUFBLElBQzlCQSxPQUFBLEdBQVUsSUFBQSxDQUFLQSxPQUFMLEdBQWVrSCxNQUFBLENBQU9qTyxNQUFBLENBQU9rTyxNQUFQLENBQWMsSUFBS25ILENBQUFBLE9BQW5CLENBQVAsRUFBb0NBLE9BQXBDLENBQXpCLENBRDhCO0FBQUEsSUFHOUIsSUFBSTJSLEtBQUEsR0FBUTNSLE9BQUEsQ0FBUTJSLEtBQXBCLENBSDhCO0FBQUEsSUFLOUIsSUFBSUEsS0FBSjtBQUFBLFFBQVdySyxPQUFBLENBQVFDLElBQVIsQ0FBYSxpQkFBYixDQUxtQixDQUFBO0FBQUEsSUFPOUIsSUFBSXZILE9BQUEsQ0FBUXlHLE9BQVIsR0FBa0IsQ0FBbEIsSUFBdUJ6RyxPQUFBLENBQVF5RyxPQUFSLEdBQWtCLEVBQTdDO0FBQUEsUUFBaUQsTUFBTSxJQUFJOEMsS0FBSixDQUFVLHFDQUFWLENBQU4sQ0FQbkI7QUFBQSxJQVE5QixJQUFJdkosT0FBQSxDQUFRNk4sU0FBUixJQUFxQjdOLE9BQUEsQ0FBUTRHLFVBQWpDO0FBQUEsUUFBNkMsTUFBTSxJQUFJMkMsS0FBSixDQUFVLG1EQUFWLENBQU4sQ0FSZjtBQUFBLElBVTlCLElBQUlqTSxRQUFBLEdBQVdrUSxPQUFBLENBQVFDLElBQVIsRUFBY3pOLE9BQWQsQ0FBZixDQVY4QjtBQUFBLElBWTlCLElBQUs0UixDQUFBQSxLQUFMLEdBQWEsRUFBYixDQVo4QjtBQUFBLElBYTlCLElBQUtDLENBQUFBLFVBQUwsR0FBa0IsRUFBbEIsQ0FiOEI7QUFBQSxJQWU5QixJQUFJRixLQUFKLEVBQVc7QUFBQSxRQUNQckssT0FBQSxDQUFRSyxPQUFSLENBQWdCLGlCQUFoQixDQURPLENBQUE7QUFBQSxRQUVQTCxPQUFBLENBQVExRCxHQUFSLENBQVksbUNBQVosRUFBaUQ1RCxPQUFBLENBQVE4UixZQUF6RCxFQUF1RTlSLE9BQUEsQ0FBUStSLGNBQS9FLENBRk8sQ0FBQTtBQUFBLFFBR1B6SyxPQUFBLENBQVFDLElBQVIsQ0FBYSxnQkFBYixDQUhPLENBQUE7QUFBQSxRQUlQLElBQUt5SyxDQUFBQSxLQUFMLEdBQWEsRUFBYixDQUpPO0FBQUEsUUFLUCxJQUFLQyxDQUFBQSxLQUFMLEdBQWEsQ0FBYixDQUxPO0FBQUEsS0FmbUI7QUFBQSxJQXVCOUIzVSxRQUFBLEdBQVc2UyxJQUFBLENBQUs3UyxRQUFMLEVBQWUwQyxPQUFmLENBQVgsQ0F2QjhCO0FBQUE7QUEwQjlCLElBQUkxQyxJQUFBQSxRQUFBLENBQVN6RCxNQUFiO0FBQUEsUUFBcUIsSUFBQSxDQUFLcVksU0FBTCxDQUFlNVUsUUFBZixFQUF5QixDQUF6QixFQUE0QixDQUE1QixFQUErQixDQUEvQixDQTFCUyxDQUFBO0FBQUEsSUE0QjlCLElBQUlxVSxLQUFKLEVBQVc7QUFBQSxRQUNQLElBQUlyVSxRQUFBLENBQVN6RCxNQUFiO0FBQUEsWUFBcUJ5TixPQUFBLENBQVExRCxHQUFSLENBQVksMEJBQVosRUFBd0MsS0FBS2dPLEtBQUwsQ0FBVyxDQUFYLENBQWNWLENBQUFBLFdBQXRELEVBQW1FLElBQUtVLENBQUFBLEtBQUwsQ0FBVyxDQUFYLENBQUEsQ0FBYy9JLFNBQWpGLENBRGQsQ0FBQTtBQUFBLFFBRVB2QixPQUFBLENBQVFLLE9BQVIsQ0FBZ0IsZ0JBQWhCLENBRk8sQ0FBQTtBQUFBLFFBR1BMLE9BQUEsQ0FBUTFELEdBQVIsQ0FBWSxrQkFBWixFQUFnQyxJQUFLcU8sQ0FBQUEsS0FBckMsRUFBNEN0WixJQUFBLENBQUtILFNBQUwsQ0FBZSxJQUFLd1osQ0FBQUEsS0FBcEIsQ0FBNUMsQ0FITyxDQUFBO0FBQUEsS0E1Qm1CO0FBQUEsQ0FYbEM7QUE4Q0FOLFNBQUEsQ0FBVWhULFNBQVYsQ0FBb0JzQixPQUFwQixHQUE4QjtBQUFBLElBQzFCeUcsT0FBQSxFQUFTLEVBRGlCO0FBQUE7QUFFMUIsSUFBQXFMLFlBQUEsRUFBYyxDQUZZO0FBQUE7QUFHMUIsSUFBQUMsY0FBQSxFQUFnQixNQUhVO0FBQUE7QUFJMUIsSUFBQW5FLFNBQUEsRUFBVyxDQUplO0FBQUE7QUFLMUIsSUFBQTlPLE1BQUEsRUFBUSxJQUxrQjtBQUFBO0FBTTFCLElBQUFzUixNQUFBLEVBQVEsRUFOa0I7QUFBQTtBQU8xQixJQUFBcEMsV0FBQSxFQUFhLEtBUGE7QUFBQTtBQVExQixJQUFBSCxTQUFBLEVBQVcsSUFSZTtBQUFBO0FBUzFCLElBQUFqSCxVQUFBLEVBQVksS0FUYztBQUFBO0FBVTFCLElBQUErSyxLQUFBLEVBQU8sQ0FBQTtBQVZtQixDQUE5QixDQTlDQTtBQTJEQUQsU0FBQSxDQUFVaFQsU0FBVixDQUFvQndULFNBQXBCLEdBQWdDLFVBQVU1VSxRQUFWLEVBQW9Cb0MsQ0FBcEIsRUFBdUJGLENBQXZCLEVBQTBCQyxDQUExQixFQUE2QjBTLEVBQTdCLEVBQWlDQyxFQUFqQyxFQUFxQ0MsRUFBckMsRUFBeUM7QUFBQSxJQUVyRSxJQUFJdk4sS0FBQSxHQUFRO0FBQUEsWUFBQ3hILFFBQUQ7QUFBQSxZQUFXb0MsQ0FBWDtBQUFBLFlBQWNGLENBQWQ7QUFBQSxZQUFpQkMsQ0FBakI7QUFBQSxTQUFaLEVBQ0lPLE9BQUEsR0FBVSxJQUFLQSxDQUFBQSxPQURuQixFQUVJMlIsS0FBQSxHQUFRM1IsT0FBQSxDQUFRMlIsS0FGcEIsQ0FGcUU7QUFBQTtBQU9yRSxJQUFPN00sT0FBQUEsS0FBQSxDQUFNakwsTUFBYixFQUFxQjtBQUFBLFFBQ2pCNEYsQ0FBQSxHQUFJcUYsS0FBQSxDQUFNRSxHQUFOLEVBQUosQ0FEaUI7QUFBQSxRQUVqQnhGLENBQUEsR0FBSXNGLEtBQUEsQ0FBTUUsR0FBTixFQUFKLENBRmlCO0FBQUEsUUFHakJ0RixDQUFBLEdBQUlvRixLQUFBLENBQU1FLEdBQU4sRUFBSixDQUhpQjtBQUFBLFFBSWpCMUgsUUFBQSxHQUFXd0gsS0FBQSxDQUFNRSxHQUFOLEVBQVgsQ0FKaUI7QUFBQSxRQU1qQixJQUFJa0YsRUFBQSxHQUFLLENBQUt4SyxJQUFBQSxDQUFkLEVBQ0k1RixFQUFBLEdBQUt3WSxJQUFBLENBQUs1UyxDQUFMLEVBQVFGLENBQVIsRUFBV0MsQ0FBWCxDQURULEVBRUlzQixJQUFBLEdBQU8sS0FBSzZRLEtBQUwsQ0FBVzlYLEVBQVgsQ0FGWCxDQU5pQjtBQUFBLFFBVWpCLElBQUksQ0FBQ2lILElBQUwsRUFBVztBQUFBLFlBQ1AsSUFBSTRRLEtBQUEsR0FBUSxDQUFaO0FBQUEsZ0JBQWVySyxPQUFBLENBQVFDLElBQVIsQ0FBYSxVQUFiLENBRFIsQ0FBQTtBQUFBLFlBR1B4RyxJQUFBLEdBQU8sSUFBQSxDQUFLNlEsS0FBTCxDQUFXOVgsRUFBWCxJQUFpQmtYLFVBQUEsQ0FBVzFULFFBQVgsRUFBcUJvQyxDQUFyQixFQUF3QkYsQ0FBeEIsRUFBMkJDLENBQTNCLEVBQThCTyxPQUE5QixDQUF4QixDQUhPO0FBQUEsWUFJUCxJQUFLNlIsQ0FBQUEsVUFBTCxDQUFnQjdYLElBQWhCLENBQXFCO0FBQUEsZ0JBQUMwRixDQUFBLEVBQUdBLENBQUo7QUFBQSxnQkFBT0YsQ0FBQSxFQUFHQSxDQUFWO0FBQUEsZ0JBQWFDLENBQUEsRUFBR0EsQ0FBaEI7QUFBQSxhQUFyQixDQUpPLENBQUE7QUFBQSxZQU1QLElBQUlrUyxLQUFKLEVBQVc7QUFBQSxnQkFDUCxJQUFJQSxLQUFBLEdBQVEsQ0FBWixFQUFlO0FBQUEsb0JBQ1hySyxPQUFBLENBQVExRCxHQUFSLENBQVksMkRBQVosRUFDSWxFLENBREosRUFDT0YsQ0FEUCxFQUNVQyxDQURWLEVBQ2FzQixJQUFBLENBQUttUSxXQURsQixFQUMrQm5RLElBQUEsQ0FBSzhILFNBRHBDLEVBQytDOUgsSUFBQSxDQUFLa1EsYUFEcEQsQ0FEVyxDQUFBO0FBQUEsb0JBR1gzSixPQUFBLENBQVFLLE9BQVIsQ0FBZ0IsVUFBaEIsQ0FIVyxDQUFBO0FBQUEsaUJBRFI7QUFBQSxnQkFNUCxJQUFJM08sR0FBQSxHQUFNLEdBQUEsR0FBTTBHLENBQWhCLENBTk87QUFBQSxnQkFPUCxJQUFLc1MsQ0FBQUEsS0FBTCxDQUFXaFosR0FBWCxJQUFtQixDQUFLZ1osSUFBQUEsQ0FBQUEsS0FBTCxDQUFXaFosR0FBWCxDQUFtQixJQUFBLENBQW5CLElBQXdCLENBQTNDLENBUE87QUFBQSxnQkFRUCxJQUFBLENBQUtpWixLQUFMLEVBUk8sQ0FBQTtBQUFBLGFBTko7QUFBQSxTQVZNO0FBQUE7QUE2QmpCLFFBQUFsUixJQUFBLENBQUs1RixNQUFMLEdBQWNtQyxRQUFkLENBN0JpQjtBQUFBO0FBZ0NqQixRQUFJLElBQUEsQ0FBQzZVLEVBQUwsRUFBUztBQUFBO0FBRUwsWUFBQSxJQUFJelMsQ0FBQSxLQUFNTSxPQUFBLENBQVE4UixZQUFkLElBQThCL1EsSUFBQSxDQUFLOEgsU0FBTCxJQUFrQjdJLE9BQUEsQ0FBUStSLGNBQTVEO0FBQUEsZ0JBQTRFLFNBQUE7QUFGdkUsU0FBVCxNQUtPO0FBQUE7QUFFSCxZQUFJclMsSUFBQUEsQ0FBQSxLQUFNTSxPQUFBLENBQVF5RyxPQUFkLElBQXlCL0csQ0FBQSxLQUFNeVMsRUFBbkM7QUFBQSxnQkFBdUMsU0FGcEM7QUFBQTtBQUtILFlBQUEsSUFBSWhVLENBQUEsR0FBSSxDQUFBLElBQU1nVSxFQUFBLEdBQUt6UyxDQUFuQixDQUxHO0FBQUEsWUFNSCxJQUFJRixDQUFBLEtBQU1wQixJQUFBLENBQUs4RixLQUFMLENBQVdrTyxFQUFBLEdBQUtqVSxDQUFoQixDQUFOLElBQTRCc0IsQ0FBQSxLQUFNckIsSUFBQSxDQUFLOEYsS0FBTCxDQUFXbU8sRUFBQSxHQUFLbFUsQ0FBaEIsQ0FBdEM7QUFBQSxnQkFBMEQsU0FOdkQ7QUFBQSxTQXJDVTtBQUFBO0FBK0NqQixRQUFBNEMsSUFBQSxDQUFLNUYsTUFBTCxHQUFjLElBQWQsQ0EvQ2lCO0FBQUEsUUFpRGpCLElBQUltQyxRQUFBLENBQVN6RCxNQUFULEtBQW9CLENBQXhCO0FBQUEsWUFBMkIsU0FqRFY7QUFBQSxRQW1EakIsSUFBSThYLEtBQUEsR0FBUSxDQUFaO0FBQUEsWUFBZXJLLE9BQUEsQ0FBUUMsSUFBUixDQUFhLFVBQWIsQ0FuREUsQ0FBQTtBQUFBO0FBc0RqQixRQUFBLElBQUlzSCxFQUFBLEdBQUssR0FBQSxHQUFNN08sT0FBQSxDQUFRb1EsTUFBZCxHQUF1QnBRLE9BQUEsQ0FBUWxCLE1BQXhDLEVBQ0lnUSxFQUFBLEdBQUssR0FBTUQsR0FBQUEsRUFEZixFQUVJMEQsRUFBQSxHQUFLLEdBQU0xRCxHQUFBQSxFQUZmLEVBR0kyRCxFQUFBLEdBQUssQ0FBSTNELEdBQUFBLEVBSGIsRUFJSTRELEVBSkosRUFJUUMsRUFKUixFQUlZQyxFQUpaLEVBSWdCQyxFQUpoQixFQUlvQnRQLElBSnBCLEVBSTBCQyxLQUoxQixDQXREaUI7QUFBQSxRQTREakJrUCxFQUFBLEdBQUtDLEVBQUEsR0FBS0MsRUFBQSxHQUFLQyxFQUFBLEdBQUssSUFBcEIsQ0E1RGlCO0FBQUEsUUE4RGpCdFAsSUFBQSxHQUFRcUwsSUFBQSxDQUFLclIsUUFBTCxFQUFlNE0sRUFBZixFQUFtQjFLLENBQUEsR0FBSXFQLEVBQXZCLEVBQTJCclAsQ0FBQSxHQUFJK1MsRUFBL0IsRUFBbUMsQ0FBbkMsRUFBc0N4UixJQUFBLENBQUsyRCxJQUEzQyxFQUFpRDNELElBQUEsQ0FBSzZELElBQXRELEVBQTRENUUsT0FBNUQsQ0FBUixDQTlEaUI7QUFBQSxRQStEakJ1RCxLQUFBLEdBQVFvTCxJQUFBLENBQUtyUixRQUFMLEVBQWU0TSxFQUFmLEVBQW1CMUssQ0FBQSxHQUFJc1AsRUFBdkIsRUFBMkJ0UCxDQUFBLEdBQUlnVCxFQUEvQixFQUFtQyxDQUFuQyxFQUFzQ3pSLElBQUEsQ0FBSzJELElBQTNDLEVBQWlEM0QsSUFBQSxDQUFLNkQsSUFBdEQsRUFBNEQ1RSxPQUE1RCxDQUFSLENBL0RpQjtBQUFBLFFBZ0VqQjFDLFFBQUEsR0FBVyxJQUFYLENBaEVpQjtBQUFBLFFBa0VqQixJQUFJZ0csSUFBSixFQUFVO0FBQUEsWUFDTm1QLEVBQUEsR0FBSzlELElBQUEsQ0FBS3JMLElBQUwsRUFBVzRHLEVBQVgsRUFBZXpLLENBQUEsR0FBSW9QLEVBQW5CLEVBQXVCcFAsQ0FBQSxHQUFJOFMsRUFBM0IsRUFBK0IsQ0FBL0IsRUFBa0N4UixJQUFBLENBQUs0RCxJQUF2QyxFQUE2QzVELElBQUEsQ0FBSzhELElBQWxELEVBQXdEN0UsT0FBeEQsQ0FBTCxDQURNO0FBQUEsWUFFTjBTLEVBQUEsR0FBSy9ELElBQUEsQ0FBS3JMLElBQUwsRUFBVzRHLEVBQVgsRUFBZXpLLENBQUEsR0FBSXFQLEVBQW5CLEVBQXVCclAsQ0FBQSxHQUFJK1MsRUFBM0IsRUFBK0IsQ0FBL0IsRUFBa0N6UixJQUFBLENBQUs0RCxJQUF2QyxFQUE2QzVELElBQUEsQ0FBSzhELElBQWxELEVBQXdEN0UsT0FBeEQsQ0FBTCxDQUZNO0FBQUEsWUFHTnNELElBQUEsR0FBTyxJQUFQLENBSE07QUFBQSxTQWxFTztBQUFBLFFBd0VqQixJQUFJQyxLQUFKLEVBQVc7QUFBQSxZQUNQb1AsRUFBQSxHQUFLaEUsSUFBQSxDQUFLcEwsS0FBTCxFQUFZMkcsRUFBWixFQUFnQnpLLENBQUEsR0FBSW9QLEVBQXBCLEVBQXdCcFAsQ0FBQSxHQUFJOFMsRUFBNUIsRUFBZ0MsQ0FBaEMsRUFBbUN4UixJQUFBLENBQUs0RCxJQUF4QyxFQUE4QzVELElBQUEsQ0FBSzhELElBQW5ELEVBQXlEN0UsT0FBekQsQ0FBTCxDQURPO0FBQUEsWUFFUDRTLEVBQUEsR0FBS2pFLElBQUEsQ0FBS3BMLEtBQUwsRUFBWTJHLEVBQVosRUFBZ0J6SyxDQUFBLEdBQUlxUCxFQUFwQixFQUF3QnJQLENBQUEsR0FBSStTLEVBQTVCLEVBQWdDLENBQWhDLEVBQW1DelIsSUFBQSxDQUFLNEQsSUFBeEMsRUFBOEM1RCxJQUFBLENBQUs4RCxJQUFuRCxFQUF5RDdFLE9BQXpELENBQUwsQ0FGTztBQUFBLFlBR1B1RCxLQUFBLEdBQVEsSUFBUixDQUhPO0FBQUEsU0F4RU07QUFBQSxRQThFakIsSUFBSW9PLEtBQUEsR0FBUSxDQUFaO0FBQUEsWUFBZXJLLE9BQUEsQ0FBUUssT0FBUixDQUFnQixVQUFoQixDQTlFRSxDQUFBO0FBQUEsUUFnRmpCN0MsS0FBQSxDQUFNOUssSUFBTixDQUFXeVksRUFBQSxJQUFNLEVBQWpCLEVBQXFCL1MsQ0FBQSxHQUFJLENBQXpCLEVBQTRCRixDQUFBLEdBQUksQ0FBaEMsRUFBdUNDLENBQUEsR0FBSSxDQUEzQyxDQWhGaUIsQ0FBQTtBQUFBLFFBaUZqQnFGLEtBQUEsQ0FBTTlLLElBQU4sQ0FBVzBZLEVBQUEsSUFBTSxFQUFqQixFQUFxQmhULENBQUEsR0FBSSxDQUF6QixFQUE0QkYsQ0FBQSxHQUFJLENBQWhDLEVBQXVDQyxDQUFBLEdBQUksQ0FBSixHQUFRLENBQS9DLENBakZpQixDQUFBO0FBQUEsUUFrRmpCcUYsS0FBQSxDQUFNOUssSUFBTixDQUFXMlksRUFBQSxJQUFNLEVBQWpCLEVBQXFCalQsQ0FBQSxHQUFJLENBQXpCLEVBQTRCRixDQUFBLEdBQUksQ0FBSixHQUFRLENBQXBDLEVBQXVDQyxDQUFBLEdBQUksQ0FBM0MsQ0FsRmlCLENBQUE7QUFBQSxRQW1GakJxRixLQUFBLENBQU05SyxJQUFOLENBQVc0WSxFQUFBLElBQU0sRUFBakIsRUFBcUJsVCxDQUFBLEdBQUksQ0FBekIsRUFBNEJGLENBQUEsR0FBSSxDQUFKLEdBQVEsQ0FBcEMsRUFBdUNDLENBQUEsR0FBSSxDQUFKLEdBQVEsQ0FBL0MsQ0FuRmlCLENBQUE7QUFBQSxLQVBnRDtBQUFBLENBQXpFLENBM0RBO0FBeUpBaVMsU0FBQSxDQUFVaFQsU0FBVixDQUFvQnVMLE9BQXBCLEdBQThCLFVBQVV2SyxDQUFWLEVBQWFGLENBQWIsRUFBZ0JDLENBQWhCLEVBQW1CO0FBQUEsSUFDN0MsSUFBSU8sT0FBQSxHQUFVLElBQUtBLENBQUFBLE9BQW5CLEVBQ0lsQixNQUFBLEdBQVNrQixPQUFBLENBQVFsQixNQURyQixFQUVJNlMsS0FBQSxHQUFRM1IsT0FBQSxDQUFRMlIsS0FGcEIsQ0FENkM7QUFBQSxJQUs3QyxJQUFJalMsQ0FBQSxHQUFJLENBQUosSUFBU0EsQ0FBQSxHQUFJLEVBQWpCO0FBQUEsUUFBcUIsT0FBTyxJQUFQLENBTHdCO0FBQUEsSUFPN0MsSUFBSXdLLEVBQUEsR0FBSyxDQUFBLElBQUt4SyxDQUFkLENBUDZDO0FBQUEsSUFRN0NGLENBQUEsR0FBSyxDQUFDQSxDQUFBLEdBQUkwSyxFQUFMLEdBQVdBLEVBQVgsSUFBaUJBLEVBQXRCLENBUjZDO0FBQUE7QUFVN0MsSUFBSXBRLElBQUFBLEVBQUEsR0FBS3dZLElBQUEsQ0FBSzVTLENBQUwsRUFBUUYsQ0FBUixFQUFXQyxDQUFYLENBQVQsQ0FWNkM7QUFBQSxJQVc3QyxJQUFJLElBQUEsQ0FBS21TLEtBQUwsQ0FBVzlYLEVBQVgsQ0FBSjtBQUFBLFFBQW9CLE9BQU8rWSxhQUFBLENBQVUsSUFBS2pCLENBQUFBLEtBQUwsQ0FBVzlYLEVBQVgsQ0FBVixFQUEwQmdGLE1BQTFCLENBQVAsQ0FYeUI7QUFBQSxJQWE3QyxJQUFJNlMsS0FBQSxHQUFRLENBQVo7QUFBQSxRQUFlckssT0FBQSxDQUFRMUQsR0FBUixDQUFZLDRCQUFaLEVBQTBDbEUsQ0FBMUMsRUFBNkNGLENBQTdDLEVBQWdEQyxDQUFoRCxDQWI4QixDQUFBO0FBQUEsSUFlN0MsSUFBSXFULEVBQUEsR0FBS3BULENBQVQsRUFDSTRPLEVBQUEsR0FBSzlPLENBRFQsRUFFSStPLEVBQUEsR0FBSzlPLENBRlQsRUFHSXNULE1BSEosQ0FmNkM7QUFBQSxJQW9CN0MsT0FBTyxDQUFDQSxNQUFELElBQVdELEVBQUEsR0FBSyxDQUF2QixFQUEwQjtBQUFBLFFBQ3RCQSxFQUFBLEVBRHNCLENBQUE7QUFBQSxRQUV0QnhFLEVBQUEsR0FBS2xRLElBQUEsQ0FBSzhGLEtBQUwsQ0FBV29LLEVBQUEsR0FBSyxDQUFoQixDQUFMLENBRnNCO0FBQUEsUUFHdEJDLEVBQUEsR0FBS25RLElBQUEsQ0FBSzhGLEtBQUwsQ0FBV3FLLEVBQUEsR0FBSyxDQUFoQixDQUFMLENBSHNCO0FBQUEsUUFJdEJ3RSxNQUFBLEdBQVMsSUFBS25CLENBQUFBLEtBQUwsQ0FBV1UsSUFBQSxDQUFLUSxFQUFMLEVBQVN4RSxFQUFULEVBQWFDLEVBQWIsQ0FBWCxDQUFULENBSnNCO0FBQUEsS0FwQm1CO0FBQUEsSUEyQjdDLElBQUksQ0FBQ3dFLE1BQUQsSUFBVyxDQUFDQSxNQUFBLENBQU81WCxNQUF2QjtBQUFBLFFBQStCLE9BQU8sSUFBUCxDQTNCYztBQUFBO0FBOEI3QyxJQUFJd1csSUFBQUEsS0FBQSxHQUFRLENBQVo7QUFBQSxRQUFlckssT0FBQSxDQUFRMUQsR0FBUixDQUFZLDZCQUFaLEVBQTJDa1AsRUFBM0MsRUFBK0N4RSxFQUEvQyxFQUFtREMsRUFBbkQsQ0E5QjhCLENBQUE7QUFBQSxJQWdDN0MsSUFBSW9ELEtBQUEsR0FBUSxDQUFaO0FBQUEsUUFBZXJLLE9BQUEsQ0FBUUMsSUFBUixDQUFhLGVBQWIsQ0FoQzhCLENBQUE7QUFBQSxJQWlDN0MsS0FBSzJLLFNBQUwsQ0FBZWEsTUFBQSxDQUFPNVgsTUFBdEIsRUFBOEIyWCxFQUE5QixFQUFrQ3hFLEVBQWxDLEVBQXNDQyxFQUF0QyxFQUEwQzdPLENBQTFDLEVBQTZDRixDQUE3QyxFQUFnREMsQ0FBaEQsQ0FqQzZDLENBQUE7QUFBQSxJQWtDN0MsSUFBSWtTLEtBQUEsR0FBUSxDQUFaO0FBQUEsUUFBZXJLLE9BQUEsQ0FBUUssT0FBUixDQUFnQixlQUFoQixDQWxDOEIsQ0FBQTtBQUFBLElBb0M3QyxPQUFPLElBQUtpSyxDQUFBQSxLQUFMLENBQVc5WCxFQUFYLENBQUEsR0FBaUIrWSxhQUFBLENBQVUsSUFBQSxDQUFLakIsS0FBTCxDQUFXOVgsRUFBWCxDQUFWLEVBQTBCZ0YsTUFBMUIsQ0FBakIsR0FBcUQsSUFBNUQsQ0FwQzZDO0FBQUEsQ0FBakQsQ0F6SkE7QUFnTUEsU0FBU3dULElBQVQsQ0FBYzVTLENBQWQsRUFBaUJGLENBQWpCLEVBQW9CQyxDQUFwQixFQUF1QjtBQUFBLElBQ25CLE9BQVMsQ0FBQyxDQUFLQyxDQUFBQSxJQUFBQSxDQUFMLElBQVVELENBQVgsR0FBZUQsQ0FBZixJQUFvQixFQUF0QixHQUE0QkUsQ0FBbkMsQ0FEbUI7QUFBQSxDQWhNdkI7QUFvTUEsU0FBU3dILE1BQVQsQ0FBZ0JzRixJQUFoQixFQUFzQkMsR0FBdEIsRUFBMkI7QUFBQSxJQUN2QixLQUFBLElBQVM3UyxDQUFULElBQWM2UyxHQUFkO0FBQUEsUUFBbUJELElBQUEsQ0FBSzVTLENBQUwsQ0FBQSxHQUFVNlMsR0FBQSxDQUFJN1MsQ0FBSixDQUFWLENBREk7QUFBQSxJQUV2QixPQUFPNFMsSUFBUCxDQUZ1QjtBQUFBOztBQ2xNM0I7QUFtREEsU0FBU3dHLGVBQVQsQ0FBeUJ0WCxNQUF6QixFQUEwREMsUUFBMUQsRUFBNEY7QUFBQSxJQUN4RixNQUFNc1gsU0FBQSxHQUFZdlgsTUFBQSxDQUFPd1gsTUFBUCxDQUFjRCxTQUFoQyxDQUR3RjtBQUFBLElBR3hGLElBQUksQ0FBQyxJQUFLRSxDQUFBQSxhQUFWLEVBQXlCO0FBQUEsUUFDckIsT0FBT3hYLFFBQUEsQ0FBUyxJQUFULEVBQWUsSUFBZixDQUFQLENBQUE7QUFEcUIsS0FIK0Q7QUFBQSxJQU94RixNQUFNeVgsV0FBQSxHQUFjLElBQUEsQ0FBS0QsYUFBTCxDQUFtQmxKLE9BQW5CLENBQTJCZ0osU0FBQSxDQUFVdlQsQ0FBckMsRUFBd0N1VCxTQUFBLENBQVV6VCxDQUFsRCxFQUFxRHlULFNBQUEsQ0FBVXhULENBQS9ELENBQXBCLENBUHdGO0FBQUEsSUFReEYsSUFBSSxDQUFDMlQsV0FBTCxFQUFrQjtBQUFBLFFBQ2QsT0FBT3pYLFFBQUEsQ0FBUyxJQUFULEVBQWUsSUFBZixDQUFQLENBQUE7QUFEYyxLQVJzRTtBQUFBLElBWXhGLE1BQU0wWCxjQUFBLEdBQWlCLElBQUl6VCxnQkFBSixDQUFtQndULFdBQUEsQ0FBWTlWLFFBQS9CLENBQXZCLENBWndGO0FBQUE7QUFpQnhGO0FBQUE7QUFBQSxJQUFBLElBQUkrRCxHQUFBLEdBQU1pUyxhQUFBLENBQU1ELGNBQU4sQ0FBVixDQWpCd0Y7QUFBQSxJQWtCeEYsSUFBSWhTLEdBQUEsQ0FBSWtTLFVBQUosS0FBbUIsQ0FBbkIsSUFBd0JsUyxHQUFBLENBQUltUyxVQUFKLEtBQW1CblMsR0FBQSxDQUFJK08sTUFBSixDQUFXb0QsVUFBMUQsRUFBc0U7QUFBQTtBQUVsRSxRQUFBblMsR0FBQSxHQUFNLElBQUlvUyxVQUFKLENBQWVwUyxHQUFmLENBQU4sQ0FGa0U7QUFBQSxLQWxCa0I7QUFBQSxJQXVCeEYxRixRQUFBLENBQVMsSUFBVCxFQUFlO0FBQUEsUUFDWCtYLFVBQUEsRUFBWUwsY0FERDtBQUFBLFFBRVhNLE9BQUEsRUFBU3RTLEdBQUEsQ0FBSStPLE1BRkY7QUFBQSxLQUFmLENBdkJ3RixDQUFBO0FBQUEsQ0FyRDVGO0FBNEZBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUEsTUFBTXdELG1CQUFOLFNBQWtDQyxxQ0FBbEMsQ0FBeUQ7QUFBQTtBQVVyRDtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUEsSUFBQTFaLFdBQUEsQ0FBWTJaLEtBQVosRUFBMEJDLFVBQTFCLEVBQXVEQyxlQUF2RCxFQUF1RkMsY0FBdkYsRUFBZ0hDLFdBQWhILEVBQTJJO0FBQUEsUUFDdkksS0FBQSxDQUFNSixLQUFOLEVBQWFDLFVBQWIsRUFBeUJDLGVBQXpCLEVBQTBDQyxjQUExQyxFQUEwRGpCLGVBQTFELENBRHVJLENBQUE7QUFBQSxRQUV2SSxJQUFJa0IsV0FBSixFQUFpQjtBQUFBLFlBQ2IsSUFBS0EsQ0FBQUEsV0FBTCxHQUFtQkEsV0FBbkIsQ0FEYTtBQUFBLFNBRnNIO0FBQUEsS0FWdEY7QUFBQTtBQWtDckQ7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQSxJQUFBQyxRQUFBLENBQVN6WSxNQUFULEVBQXdDQyxRQUF4QyxFQUFnSTtBQUFBLFFBQzVILE1BQU15WSxZQUFBLEdBQWUxWSxNQUFBLElBQVVBLE1BQUEsQ0FBTzJZLE9BQXRDLENBRDRIO0FBQUEsUUFFNUgsTUFBTUMsSUFBQSxHQUFPRixZQUFBLElBQWdCQSxZQUFBLENBQWFHLHFCQUExQyxDQUY0SDtBQUFBLFFBSTVILElBQUEsQ0FBS0wsV0FBTCxDQUFpQnhZLE1BQWpCLEVBQXlCLENBQUNzQyxHQUFELEVBQWN5UCxJQUFkLEtBQWdDO0FBQUEsWUFDckQsSUFBSXpQLEdBQUEsSUFBTyxDQUFDeVAsSUFBWixFQUFrQjtBQUFBLGdCQUNkLE9BQU85UixRQUFBLENBQVNxQyxHQUFULENBQVAsQ0FEYztBQUFBLGFBQWxCLE1BRU8sSUFBSSxPQUFPeVAsSUFBUCxLQUFnQixRQUFwQixFQUE4QjtBQUFBLGdCQUNqQyxPQUFPOVIsUUFBQSxDQUFTLElBQUk0TixLQUFKLENBQVUsQ0FBQyxxQkFBRCxHQUF3QjdOLE1BQUEsQ0FBT1AsTUFBL0IsRUFBc0MsZ0NBQXRDLENBQVYsQ0FBVCxDQUFQLENBRGlDO0FBQUEsYUFBOUIsTUFFQTtBQUFBLGdCQUNIK0IsYUFBQSxDQUFPdVEsSUFBUCxFQUFhLElBQWIsQ0FERyxDQUFBO0FBQUEsZ0JBR0gsSUFBSTtBQUFBLG9CQUNBLElBQUkvUixNQUFBLENBQU84WSxNQUFYLEVBQW1CO0FBQUEsd0JBQ2YsTUFBTUMsUUFBQSxHQUFXQywrQkFBQSxDQUFpQmhaLE1BQUEsQ0FBTzhZLE1BQXhCLEVBQWdDO0FBQUEsNEJBQUNuWCxJQUFBLEVBQU0sU0FBUDtBQUFBLDRCQUFrQixpQkFBaUIsYUFBbkM7QUFBQSw0QkFBa0RzWCxXQUFBLEVBQWEsS0FBL0Q7QUFBQSw0QkFBc0VDLFVBQUEsRUFBWSxLQUFsRjtBQUFBLHlCQUFoQyxDQUFqQixDQURlO0FBQUEsd0JBRWYsSUFBSUgsUUFBQSxDQUFTeGEsTUFBVCxLQUFvQixPQUF4QjtBQUFBLDRCQUNJLE1BQU0sSUFBSXNQLEtBQUosQ0FBVWtMLFFBQUEsQ0FBU3hTLEtBQVQsQ0FBZWpILEdBQWYsQ0FBbUJnRCxHQUFBLElBQU8sQ0FBR0EsR0FBQUEsR0FBQSxDQUFJaEYsR0FBUCxFQUFXLEVBQVgsR0FBZWdGLEdBQUEsQ0FBSTZXLE9BQW5CLEVBQTFCLENBQUEsQ0FBQSxDQUF3REMsSUFBeEQsQ0FBNkQsSUFBN0QsQ0FBVixDQUFOLENBSFc7QUFBQSx3QkFLZixNQUFNeFgsUUFBQSxHQUFXbVEsSUFBQSxDQUFLblEsUUFBTCxDQUFja1gsTUFBZCxDQUFxQjVWLE9BQUEsSUFBVzZWLFFBQUEsQ0FBU3hTLEtBQVQsQ0FBZThTLFFBQWYsQ0FBd0IsRUFBQy9NLElBQUEsRUFBTSxDQUFQLEVBQXhCLEVBQW1DcEosT0FBbkMsQ0FBaEMsQ0FBakIsQ0FMZTtBQUFBLHdCQU1mNk8sSUFBQSxHQUFPO0FBQUEsNEJBQUNwUSxJQUFBLEVBQU0sbUJBQVA7QUFBQSw0QkFBNEJDLFFBQTVCO0FBQUEseUJBQVAsQ0FOZTtBQUFBLHFCQURuQjtBQUFBLG9CQVVBLElBQUs2VixDQUFBQSxhQUFMLEdBQXFCelgsTUFBQSxDQUFPaVAsT0FBUCxHQUNqQixJQUFJMUQsWUFBSixDQUFpQitOLHNCQUFBLENBQXVCdFosTUFBdkIsQ0FBakIsQ0FBQSxDQUFpRDJMLElBQWpELENBQXNEb0csSUFBQSxDQUFLblEsUUFBM0QsQ0FEaUIsR0FFakJtVSxTQUFBLENBQVVoRSxJQUFWLEVBQWdCL1IsTUFBQSxDQUFPdVosZ0JBQXZCLENBRkosQ0FWQTtBQUFBLGlCQUFKLENBYUUsT0FBT2pYLEdBQVAsRUFBWTtBQUFBLG9CQUNWLE9BQU9yQyxRQUFBLENBQVNxQyxHQUFULENBQVAsQ0FEVTtBQUFBLGlCQWhCWDtBQUFBLGdCQW9CSCxJQUFLa1gsQ0FBQUEsTUFBTCxHQUFjLEVBQWQsQ0FwQkc7QUFBQSxnQkFzQkgsTUFBTWpiLE1BQUEsR0FBUyxFQUFmLENBdEJHO0FBQUEsZ0JBdUJILElBQUlxYSxJQUFKLEVBQVU7QUFBQSxvQkFDTixNQUFNYSxrQkFBQSxHQUFxQkMsd0NBQUEsQ0FBMEJoQixZQUExQixDQUEzQixDQURNO0FBQUE7QUFJTjtBQUFBLG9CQUFBLElBQUllLGtCQUFKLEVBQXdCO0FBQUEsd0JBQ3BCbGIsTUFBQSxDQUFPb2IsY0FBUCxHQUF3QixFQUF4QixDQURvQjtBQUFBLHdCQUVwQnBiLE1BQUEsQ0FBT29iLGNBQVAsQ0FBc0IzWixNQUFBLENBQU9QLE1BQTdCLENBQXVDeEMsR0FBQUEsSUFBQSxDQUFLMmMsS0FBTCxDQUFXM2MsSUFBQSxDQUFLSCxTQUFMLENBQWUyYyxrQkFBZixDQUFYLENBQXZDLENBRm9CO0FBQUEscUJBSmxCO0FBQUEsaUJBdkJQO0FBQUEsZ0JBZ0NIeFosUUFBQSxDQUFTLElBQVQsRUFBZTFCLE1BQWYsQ0FoQ0csQ0FBQTtBQUFBLGFBTDhDO0FBQUEsU0FBekQsQ0FKNEgsQ0FBQTtBQUFBLEtBbEMzRTtBQUFBO0FBMEZyRDtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQSxJQUFBc2IsVUFBQSxDQUFXN1osTUFBWCxFQUF5Q0MsUUFBekMsRUFBNkU7QUFBQSxRQUN6RSxNQUFNdVosTUFBQSxHQUFTLElBQUtBLENBQUFBLE1BQXBCLEVBQ0l0WixHQUFBLEdBQU1GLE1BQUEsQ0FBT0UsR0FEakIsQ0FEeUU7QUFBQSxRQUl6RSxJQUFJc1osTUFBQSxJQUFVQSxNQUFBLENBQU90WixHQUFQLENBQWQsRUFBMkI7QUFBQSxZQUN2QixPQUFPLE1BQU0yWixVQUFOLENBQWlCN1osTUFBakIsRUFBeUJDLFFBQXpCLENBQVAsQ0FEdUI7QUFBQSxTQUEzQixNQUVPO0FBQUEsWUFDSCxPQUFPLEtBQUtGLFFBQUwsQ0FBY0MsTUFBZCxFQUFzQkMsUUFBdEIsQ0FBUCxDQURHO0FBQUEsU0FOa0U7QUFBQSxLQTFGeEI7QUFBQTtBQWlIckQ7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBLElBQUF1WSxXQUFBLENBQVl4WSxNQUFaLEVBQTJDQyxRQUEzQyxFQUFxRjtBQUFBO0FBS2pGO0FBQUE7QUFBQTtBQUFBLFFBQUlELElBQUFBLE1BQUEsQ0FBTzJZLE9BQVgsRUFBb0I7QUFBQSxZQUNoQm1CLHNCQUFBLENBQVE5WixNQUFBLENBQU8yWSxPQUFmLEVBQXdCMVksUUFBeEIsQ0FEZ0IsQ0FBQTtBQUFBLFNBQXBCLE1BRU8sSUFBSSxPQUFPRCxNQUFBLENBQU8rUixJQUFkLEtBQXVCLFFBQTNCLEVBQXFDO0FBQUEsWUFDeEMsSUFBSTtBQUFBLGdCQUNBLE9BQU85UixRQUFBLENBQVMsSUFBVCxFQUFlaEQsSUFBQSxDQUFLMmMsS0FBTCxDQUFXNVosTUFBQSxDQUFPK1IsSUFBbEIsQ0FBZixDQUFQLENBREE7QUFBQSxhQUFKLENBRUUsT0FBT2dJLENBQVAsRUFBVTtBQUFBLGdCQUNSLE9BQU85WixRQUFBLENBQVMsSUFBSTROLEtBQUosQ0FBVSxDQUFDLHFCQUFELEdBQXdCN04sTUFBQSxDQUFPUCxNQUEvQixFQUFzQyxnQ0FBdEMsQ0FBVixDQUFULENBQVAsQ0FEUTtBQUFBLGFBSDRCO0FBQUEsU0FBckMsTUFNQTtBQUFBLFlBQ0gsT0FBT1EsUUFBQSxDQUFTLElBQUk0TixLQUFKLENBQVUsQ0FBQyxxQkFBRCxHQUF3QjdOLE1BQUEsQ0FBT1AsTUFBL0IsRUFBc0MsZ0NBQXRDLENBQVYsQ0FBVCxDQUFQLENBREc7QUFBQSxTQWIwRTtBQUFBLEtBakhoQztBQUFBLElBbUlyRG1QLHVCQUFBLENBQXdCNU8sTUFBeEIsRUFBcURDLFFBQXJELEVBQWlGO0FBQUEsUUFDN0UsSUFBSTtBQUFBLFlBQ0FBLFFBQUEsQ0FBUyxJQUFULEVBQWUsSUFBS3dYLENBQUFBLGFBQUwsQ0FBbUI3SSx1QkFBbkIsQ0FBMkM1TyxNQUFBLENBQU91TixTQUFsRCxDQUFmLENBREEsQ0FBQTtBQUFBLFNBQUosQ0FFRSxPQUFPd00sQ0FBUCxFQUFVO0FBQUEsWUFDUjlaLFFBQUEsQ0FBUzhaLENBQVQsQ0FEUSxDQUFBO0FBQUEsU0FIaUU7QUFBQSxLQW5JNUI7QUFBQSxJQTJJckRDLGtCQUFBLENBQW1CaGEsTUFBbkIsRUFBZ0RDLFFBQWhELEVBQTJGO0FBQUEsUUFDdkYsSUFBSTtBQUFBLFlBQ0FBLFFBQUEsQ0FBUyxJQUFULEVBQWUsSUFBS3dYLENBQUFBLGFBQUwsQ0FBbUJuSyxXQUFuQixDQUErQnROLE1BQUEsQ0FBT3VOLFNBQXRDLENBQWYsQ0FEQSxDQUFBO0FBQUEsU0FBSixDQUVFLE9BQU93TSxDQUFQLEVBQVU7QUFBQSxZQUNSOVosUUFBQSxDQUFTOFosQ0FBVCxDQURRLENBQUE7QUFBQSxTQUgyRTtBQUFBLEtBM0l0QztBQUFBLElBbUpyREUsZ0JBQUEsQ0FBaUJqYSxNQUFqQixFQUE2RUMsUUFBN0UsRUFBd0g7QUFBQSxRQUNwSCxJQUFJO0FBQUEsWUFDQUEsUUFBQSxDQUFTLElBQVQsRUFBZSxLQUFLd1gsYUFBTCxDQUFtQnZKLFNBQW5CLENBQTZCbE8sTUFBQSxDQUFPdU4sU0FBcEMsRUFBK0N2TixNQUFBLENBQU9tTyxLQUF0RCxFQUE2RG5PLE1BQUEsQ0FBT29PLE1BQXBFLENBQWYsQ0FEQSxDQUFBO0FBQUEsU0FBSixDQUVFLE9BQU8yTCxDQUFQLEVBQVU7QUFBQSxZQUNSOVosUUFBQSxDQUFTOFosQ0FBVCxDQURRLENBQUE7QUFBQSxTQUh3RztBQUFBLEtBbkpuRTtBQUFBLENBNUZ6RDtBQXdQQSxTQUFTVCxzQkFBVCxDQUFnQyxDQUFDWSxtQkFBRCxFQUFzQm5LLGlCQUF0QixDQUFoQyxFQUEwRTtBQUFBLElBQ3RFLElBQUksQ0FBQ0EsaUJBQUQsSUFBc0IsQ0FBQ21LLG1CQUEzQjtBQUFBLFFBQWdELE9BQU9BLG1CQUFQLENBRHNCO0FBQUEsSUFHdEUsTUFBTUMsY0FBQSxHQUFpQixFQUF2QixDQUhzRTtBQUFBLElBSXRFLE1BQU1DLGlCQUFBLEdBQW9CLEVBQTFCLENBSnNFO0FBQUEsSUFLdEUsTUFBTUMsT0FBQSxHQUFVO0FBQUEsUUFBQ0MsV0FBQSxFQUFhLElBQWQ7QUFBQSxRQUFvQmhPLElBQUEsRUFBTSxDQUExQjtBQUFBLEtBQWhCLENBTHNFO0FBQUEsSUFNdEUsTUFBTXBKLE9BQUEsR0FBVSxFQUFDSSxVQUFBLEVBQVksSUFBYixFQUFoQixDQU5zRTtBQUFBLElBT3RFLE1BQU1pWCxhQUFBLEdBQWdCaGQsTUFBQSxDQUFPQyxJQUFQLENBQVl1UyxpQkFBWixDQUF0QixDQVBzRTtBQUFBLElBU3RFLEtBQVd6UyxNQUFBQSxHQUFYLElBQWtCaWQsYUFBbEIsRUFBaUM7QUFBQSxRQUM3QixNQUFNLENBQUNDLFFBQUQsRUFBV0MsYUFBWCxJQUE0QjFLLGlCQUFBLENBQWtCelMsR0FBbEIsQ0FBbEMsQ0FENkI7QUFBQSxRQUc3QixNQUFNb2QsbUJBQUEsR0FBc0IxQiwrQkFBQSxDQUFpQnlCLGFBQWpCLENBQTVCLENBSDZCO0FBQUEsUUFJN0IsTUFBTUUsc0JBQUEsR0FBeUIzQiwrQkFBQSxDQUMzQixPQUFPd0IsUUFBUCxLQUFvQixRQUFwQixHQUErQjtBQUFBLFlBQUNBLFFBQUQ7QUFBQSxZQUFXLENBQUMsYUFBRCxDQUFYO0FBQUEsWUFBNEI7QUFBQSxnQkFBQyxLQUFEO0FBQUEsZ0JBQVFsZCxHQUFSO0FBQUEsYUFBNUI7QUFBQSxTQUEvQixHQUEyRWtkLFFBRGhELENBQS9CLENBSjZCO0FBQUEsUUFVN0JMLGNBQUEsQ0FBZTdjLEdBQWYsSUFBc0JvZCxtQkFBQSxDQUFvQm5VLEtBQTFDLENBVjZCO0FBQUEsUUFXN0I2VCxpQkFBQSxDQUFrQjljLEdBQWxCLElBQXlCcWQsc0JBQUEsQ0FBdUJwVSxLQUFoRCxDQVg2QjtBQUFBLEtBVHFDO0FBQUEsSUF1QnRFMlQsbUJBQUEsQ0FBb0I1YSxHQUFwQixHQUEyQnNiLGVBQUQsSUFBcUI7QUFBQSxRQUMzQzFYLE9BQUEsQ0FBUUksVUFBUixHQUFxQnNYLGVBQXJCLENBRDJDO0FBQUEsUUFFM0MsTUFBTXRYLFVBQUEsR0FBYSxFQUFuQixDQUYyQztBQUFBLFFBRzNDLEtBQVdoRyxNQUFBQSxHQUFYLElBQWtCaWQsYUFBbEIsRUFBaUM7QUFBQSxZQUM3QmpYLFVBQUEsQ0FBV2hHLEdBQVgsQ0FBQSxHQUFrQjZjLGNBQUEsQ0FBZTdjLEdBQWYsQ0FBb0IrYixDQUFBQSxRQUFwQixDQUE2QmdCLE9BQTdCLEVBQXNDblgsT0FBdEMsQ0FBbEIsQ0FENkI7QUFBQSxTQUhVO0FBQUEsUUFNM0MsT0FBT0ksVUFBUCxDQU4yQztBQUFBLEtBQS9DLENBdkJzRTtBQUFBLElBK0J0RTRXLG1CQUFBLENBQW9CL08sTUFBcEIsR0FBNkIsQ0FBQ21QLFdBQUQsRUFBY3ZLLGlCQUFkLEtBQW9DO0FBQUEsUUFDN0Q3TSxPQUFBLENBQVFJLFVBQVIsR0FBcUJ5TSxpQkFBckIsQ0FENkQ7QUFBQSxRQUU3RCxLQUFXelMsTUFBQUEsR0FBWCxJQUFrQmlkLGFBQWxCLEVBQWlDO0FBQUEsWUFDN0JGLE9BQUEsQ0FBUUMsV0FBUixHQUFzQkEsV0FBQSxDQUFZaGQsR0FBWixDQUF0QixDQUQ2QjtBQUFBLFlBRTdCZ2QsV0FBQSxDQUFZaGQsR0FBWixDQUFBLEdBQW1COGMsaUJBQUEsQ0FBa0I5YyxHQUFsQixDQUF1QitiLENBQUFBLFFBQXZCLENBQWdDZ0IsT0FBaEMsRUFBeUNuWCxPQUF6QyxDQUFuQixDQUY2QjtBQUFBLFNBRjRCO0FBQUEsS0FBakUsQ0EvQnNFO0FBQUEsSUF1Q3RFLE9BQU9nWCxtQkFBUCxDQXZDc0U7QUFBQTs7QUN0UDFFO0FBZ0NBO0FBQUE7QUFBQTtBQUFlLE1BQU1XLE1BQU4sQ0FBYTtBQUFBLElBY3hCcGMsV0FBQSxDQUFZcWMsSUFBWixFQUE4QztBQUFBLFFBRTFDLElBQUtBLENBQUFBLElBQUwsR0FBWUEsSUFBWixDQUYwQztBQUFBLFFBRzFDLElBQUEsQ0FBSzFDLEtBQUwsR0FBYSxJQUFJMkMsb0JBQUosQ0FBVUQsSUFBVixFQUFnQixJQUFoQixDQUFiLENBSDBDO0FBQUEsUUFLMUMsSUFBS0UsQ0FBQUEsWUFBTCxHQUFvQixFQUFwQixDQUwwQztBQUFBLFFBTTFDLElBQUsxQyxDQUFBQSxlQUFMLEdBQXVCLEVBQXZCLENBTjBDO0FBQUEsUUFPMUMsSUFBS0MsQ0FBQUEsY0FBTCxHQUFzQixFQUF0QixDQVAwQztBQUFBLFFBUzFDLElBQUswQyxDQUFBQSxXQUFMLEdBQW1CLEVBQW5CLENBVDBDO0FBQUEsUUFVMUMsSUFBQSxDQUFLQyxpQkFBTCxHQUF5QkMsNEJBQUEsQ0FBYyxFQUFDaFgsSUFBQSxFQUFNLFVBQVAsRUFBZCxDQUF6QixDQVYwQztBQUFBLFFBWTFDLElBQUEsQ0FBS2lYLGlCQUFMLEdBQXlCO0FBQUEsWUFDckJDLE1BQUEsRUFBUWxELHFDQURhO0FBQUEsWUFFckJsRyxPQUFBLEVBQVNpRyxtQkFGWTtBQUFBLFNBQXpCLENBWjBDO0FBQUE7QUFrQjFDLFFBQUtvRCxJQUFBQSxDQUFBQSxhQUFMLEdBQXFCLEVBQXJCLENBbEIwQztBQUFBLFFBbUIxQyxJQUFLQyxDQUFBQSxnQkFBTCxHQUF3QixFQUF4QixDQW5CMEM7QUFBQSxRQXFCMUMsSUFBQSxDQUFLVCxJQUFMLENBQVVVLG9CQUFWLEdBQWlDLENBQUNyWCxJQUFELEVBQWVzWCxZQUFmLEtBQXFEO0FBQUEsWUFDbEYsSUFBSSxJQUFLTCxDQUFBQSxpQkFBTCxDQUF1QmpYLElBQXZCLENBQUosRUFBa0M7QUFBQSxnQkFDOUIsTUFBTSxJQUFJMEosS0FBSixDQUFVLENBQUMseUJBQUQsR0FBNEIxSixJQUE1QixFQUFpQyxxQkFBakMsQ0FBVixDQUFOLENBRDhCO0FBQUEsYUFEZ0Q7QUFBQSxZQUlsRixLQUFLaVgsaUJBQUwsQ0FBdUJqWCxJQUF2QixDQUFBLEdBQStCc1gsWUFBL0IsQ0FKa0Y7QUFBQSxTQUF0RixDQXJCMEM7QUFBQTtBQTZCMUMsUUFBQSxJQUFBLENBQUtYLElBQUwsQ0FBVVkscUJBQVYsR0FBbUNDLGFBQUQsSUFBa0k7QUFBQSxZQUNoSyxJQUFJQyxxQkFBQSxDQUFvQkMsUUFBcEIsRUFBSixFQUFvQztBQUFBLGdCQUNoQyxNQUFNLElBQUloTyxLQUFKLENBQVUscUNBQVYsQ0FBTixDQURnQztBQUFBLGFBRDRIO0FBQUEsWUFJaEsrTixxQkFBQSxDQUFvQixvQkFBcEIsSUFBNENELGFBQUEsQ0FBY0csa0JBQTFELENBSmdLO0FBQUEsWUFLaEtGLHFCQUFBLENBQW9CLDBCQUFwQixJQUFrREQsYUFBQSxDQUFjSSx3QkFBaEUsQ0FMZ0s7QUFBQSxZQU1oS0gscUJBQUEsQ0FBb0IsZ0NBQXBCLElBQXdERCxhQUFBLENBQWNLLDhCQUF0RSxDQU5nSztBQUFBLFNBQXBLLENBN0IwQztBQUFBLEtBZHRCO0FBQUEsSUFxRHhCQyxXQUFBLENBQVlDLEtBQVosRUFBMkJDLE1BQTNCLEVBQTBDbGMsUUFBMUMsRUFBd0U7QUFBQSxRQUNwRSxPQUFPLElBQUsrYSxDQUFBQSxZQUFMLENBQWtCa0IsS0FBbEIsQ0FBUCxDQURvRTtBQUFBLFFBRXBFLE9BQU8sSUFBSzVELENBQUFBLGVBQUwsQ0FBcUI0RCxLQUFyQixDQUFQLENBRm9FO0FBQUEsUUFHcEUsT0FBTyxJQUFLWixDQUFBQSxhQUFMLENBQW1CWSxLQUFuQixDQUFQLENBSG9FO0FBQUEsUUFJcEUsT0FBTyxJQUFLWCxDQUFBQSxnQkFBTCxDQUFzQlcsS0FBdEIsQ0FBUCxDQUpvRTtBQUFBLFFBS3BFamMsUUFBQSxFQUxvRSxDQUFBO0FBQUEsS0FyRGhEO0FBQUEsSUE2RHhCbWMsWUFBQSxDQUFhQyxLQUFiLEVBQTRCRixNQUE1QixFQUEyQ2xjLFFBQTNDLEVBQXlFO0FBQUE7QUFFckUsUUFBQUEsUUFBQSxFQUZxRSxDQUFBO0FBQUEsS0E3RGpEO0FBQUEsSUFrRXhCcWMsV0FBQSxDQUFZRCxLQUFaLEVBQTJCRSxRQUEzQixFQUE2QztBQUFBLFFBQ3pDLElBQUtBLENBQUFBLFFBQUwsR0FBZ0JBLFFBQWhCLENBRHlDO0FBQUEsS0FsRXJCO0FBQUEsSUFzRXhCQyxZQUFBLENBQWFOLEtBQWIsRUFBNEJPLElBQTVCLEVBQTJDO0FBQUEsUUFDdkMsS0FBS2xFLGNBQUwsQ0FBb0IyRCxLQUFwQixDQUFBLEdBQTZCTyxJQUE3QixDQUR1QztBQUFBLFFBRXZDLEtBQUEsTUFBV0MsWUFBWCxJQUEyQixJQUFBLENBQUtwQixhQUFMLENBQW1CWSxLQUFuQixDQUEzQixFQUFzRDtBQUFBLFlBQ2xELE1BQU1TLEVBQUEsR0FBSyxJQUFLckIsQ0FBQUEsYUFBTCxDQUFtQlksS0FBbkIsQ0FBQSxDQUEwQlEsWUFBMUIsQ0FBWCxDQURrRDtBQUFBLFlBRWxELEtBQVdqZCxNQUFBQSxNQUFYLElBQXFCa2QsRUFBckIsRUFBeUI7QUFBQSxnQkFDckIsSUFBSUEsRUFBQSxDQUFHbGQsTUFBSCxDQUFBLFlBQXNCMFkscUNBQTFCLEVBQWtEO0FBQUEsb0JBQzlDd0UsRUFBQSxDQUFHbGQsTUFBSCxFQUFXOFksY0FBWCxHQUE0QmtFLElBQTVCLENBRDhDO0FBQUEsb0JBRTlDRSxFQUFBLENBQUdsZCxNQUFILENBQVdtZCxDQUFBQSxJQUFYLENBQWdCLElBQUlDLG9CQUFKLENBQVUsZ0JBQVYsQ0FBaEIsQ0FGOEMsQ0FBQTtBQUFBLGlCQUQ3QjtBQUFBLGFBRnlCO0FBQUEsU0FGZjtBQUFBLEtBdEVuQjtBQUFBLElBbUZ4QkMsU0FBQSxDQUFVWixLQUFWLEVBQXlCYSxNQUF6QixFQUFnRDljLFFBQWhELEVBQThFO0FBQUEsUUFDMUUsS0FBS3FZLGVBQUwsQ0FBcUI0RCxLQUFyQixDQUFBLEdBQThCYSxNQUE5QixDQUQwRTtBQUFBLFFBRTFFLEtBQUEsTUFBV0wsWUFBWCxJQUEyQixJQUFBLENBQUtwQixhQUFMLENBQW1CWSxLQUFuQixDQUEzQixFQUFzRDtBQUFBLFlBQ2xELE1BQU1TLEVBQUEsR0FBSyxJQUFLckIsQ0FBQUEsYUFBTCxDQUFtQlksS0FBbkIsQ0FBQSxDQUEwQlEsWUFBMUIsQ0FBWCxDQURrRDtBQUFBLFlBRWxELEtBQVdqZCxNQUFBQSxNQUFYLElBQXFCa2QsRUFBckIsRUFBeUI7QUFBQSxnQkFDckJBLEVBQUEsQ0FBR2xkLE1BQUgsRUFBVzZZLGVBQVgsR0FBNkJ5RSxNQUE3QixDQURxQjtBQUFBLGFBRnlCO0FBQUEsU0FGb0I7QUFBQSxRQVExRTljLFFBQUEsRUFSMEUsQ0FBQTtBQUFBLEtBbkZ0RDtBQUFBLElBOEZ4QitjLGFBQUEsQ0FBY2QsS0FBZCxFQUE2QmUsTUFBN0IsRUFBOENoZCxRQUE5QyxFQUE0RTtBQUFBLFFBQ3hFLElBQUtpZCxDQUFBQSxPQUFMLEdBQWVELE1BQWYsQ0FEd0U7QUFBQSxRQUV4RWhkLFFBQUEsRUFGd0UsQ0FBQTtBQUFBLEtBOUZwRDtBQUFBLElBbUd4QmtkLGFBQUEsQ0FBY2pCLEtBQWQsRUFBNkJrQixNQUE3QixFQUE4RDtBQUFBLFFBQzFELElBQUEsQ0FBS25DLFdBQUwsQ0FBaUJpQixLQUFqQixJQUEwQmYsNEJBQUEsQ0FBY2lDLE1BQWQsQ0FBMUIsQ0FEMEQ7QUFBQSxLQW5HdEM7QUFBQSxJQXVHeEJDLFNBQUEsQ0FBVW5CLEtBQVYsRUFBeUJuZSxNQUF6QixFQUE0RGtDLFFBQTVELEVBQTBGO0FBQUEsUUFDdEYsSUFBQSxDQUFLcWQsYUFBTCxDQUFtQnBCLEtBQW5CLEVBQTBCdGQsT0FBMUIsQ0FBa0NiLE1BQWxDLENBRHNGLENBQUE7QUFBQSxRQUV0RmtDLFFBQUEsRUFGc0YsQ0FBQTtBQUFBLEtBdkdsRTtBQUFBLElBNEd4QnNkLFlBQUEsQ0FBYXJCLEtBQWIsRUFBNEJsYyxNQUE1QixFQUFvR0MsUUFBcEcsRUFBa0k7QUFBQSxRQUM5SCxJQUFLcWQsQ0FBQUEsYUFBTCxDQUFtQnBCLEtBQW5CLEVBQTBCbmQsTUFBMUIsQ0FBaUNpQixNQUFBLENBQU9qQyxNQUF4QyxFQUFnRGlDLE1BQUEsQ0FBT2hCLFVBQXZELENBRDhILENBQUE7QUFBQSxRQUU5SGlCLFFBQUEsRUFGOEgsQ0FBQTtBQUFBLEtBNUcxRztBQUFBLElBaUh4QkYsUUFBQSxDQUFTbWMsS0FBVCxFQUF3QmxjLE1BQXhCLEVBQXVFQyxRQUF2RSxFQUFxRztBQUFBLFFBRWpHLE1BQU1pSyxDQUFBLEdBQUksS0FBSzhTLGFBQUwsR0FBcUJ4UixxQkFBQSxDQUFPLEVBQUN3UixhQUFBLEVBQWUsSUFBQSxDQUFLRSxPQUFyQixFQUFQLEVBQXNDbGQsTUFBdEMsQ0FBckIsR0FBcUVBLE1BQS9FLENBRmlHO0FBQUEsUUFHakdrSyxDQUFBLENBQUVzVCxVQUFGLEdBQWUsSUFBQSxDQUFLdkMsV0FBTCxDQUFpQmlCLEtBQWpCLENBQUEsSUFBMkIsSUFBS2hCLENBQUFBLGlCQUEvQyxDQUhpRztBQUFBLFFBSWpHLEtBQUt1QyxlQUFMLENBQXFCdkIsS0FBckIsRUFBNEJsYyxNQUFBLENBQU8yQixJQUFuQyxFQUF5QzNCLE1BQUEsQ0FBT1AsTUFBaEQsQ0FBd0RNLENBQUFBLFFBQXhELENBQWlFbUssQ0FBakUsRUFBb0VqSyxRQUFwRSxDQUppRyxDQUFBO0FBQUEsS0FqSDdFO0FBQUEsSUF3SHhCeWQsV0FBQSxDQUFZeEIsS0FBWixFQUEyQmxjLE1BQTNCLEVBQTREQyxRQUE1RCxFQUE2RjtBQUFBLFFBQ3pGLE1BQU1pSyxDQUFBLEdBQUksS0FBSzhTLGFBQUwsR0FBcUJ4UixxQkFBQSxDQUFPLEVBQUNsTCxhQUFBLEVBQWUsSUFBQSxDQUFLNGMsT0FBckIsRUFBUCxFQUFzQ2xkLE1BQXRDLENBQXJCLEdBQXFFQSxNQUEvRSxDQUR5RjtBQUFBLFFBRXpGLElBQUsyZCxDQUFBQSxrQkFBTCxDQUF3QnpCLEtBQXhCLEVBQStCbGMsTUFBQSxDQUFPUCxNQUF0QyxDQUFBLENBQThDTSxRQUE5QyxDQUF1RG1LLENBQXZELEVBQTBEakssUUFBMUQsQ0FGeUYsQ0FBQTtBQUFBLEtBeEhyRTtBQUFBLElBNkh4QjRaLFVBQUEsQ0FBV3FDLEtBQVgsRUFBMEJsYyxNQUExQixFQUF5RUMsUUFBekUsRUFBdUc7QUFBQSxRQUVuRyxNQUFNaUssQ0FBQSxHQUFJLEtBQUs4UyxhQUFMLEdBQXFCeFIscUJBQUEsQ0FBTyxFQUFDd1IsYUFBQSxFQUFlLElBQUEsQ0FBS0UsT0FBckIsRUFBUCxFQUFzQ2xkLE1BQXRDLENBQXJCLEdBQXFFQSxNQUEvRSxDQUZtRztBQUFBLFFBR25Ha0ssQ0FBQSxDQUFFc1QsVUFBRixHQUFlLElBQUEsQ0FBS3ZDLFdBQUwsQ0FBaUJpQixLQUFqQixDQUFBLElBQTJCLElBQUtoQixDQUFBQSxpQkFBL0MsQ0FIbUc7QUFBQSxRQUluRyxLQUFLdUMsZUFBTCxDQUFxQnZCLEtBQXJCLEVBQTRCbGMsTUFBQSxDQUFPMkIsSUFBbkMsRUFBeUMzQixNQUFBLENBQU9QLE1BQWhELENBQXdEb2EsQ0FBQUEsVUFBeEQsQ0FBbUUzUCxDQUFuRSxFQUFzRWpLLFFBQXRFLENBSm1HLENBQUE7QUFBQSxLQTdIL0U7QUFBQSxJQW9JeEIyZCxTQUFBLENBQVUxQixLQUFWLEVBQXlCbGMsTUFBekIsRUFBa0VDLFFBQWxFLEVBQWdHO0FBQUEsUUFFNUYsS0FBS3dkLGVBQUwsQ0FBcUJ2QixLQUFyQixFQUE0QmxjLE1BQUEsQ0FBTzJCLElBQW5DLEVBQXlDM0IsTUFBQSxDQUFPUCxNQUFoRCxDQUF3RG1lLENBQUFBLFNBQXhELENBQWtFNWQsTUFBbEUsRUFBMEVDLFFBQTFFLENBRjRGLENBQUE7QUFBQSxLQXBJeEU7QUFBQSxJQXlJeEI0ZCxVQUFBLENBQVczQixLQUFYLEVBQTBCbGMsTUFBMUIsRUFBbUVDLFFBQW5FLEVBQWlHO0FBQUEsUUFFN0YsS0FBS3dkLGVBQUwsQ0FBcUJ2QixLQUFyQixFQUE0QmxjLE1BQUEsQ0FBTzJCLElBQW5DLEVBQXlDM0IsTUFBQSxDQUFPUCxNQUFoRCxDQUF3RG9lLENBQUFBLFVBQXhELENBQW1FN2QsTUFBbkUsRUFBMkVDLFFBQTNFLENBRjZGLENBQUE7QUFBQSxLQXpJekU7QUFBQSxJQThJeEI2ZCxZQUFBLENBQWE1QixLQUFiLEVBQTRCbGMsTUFBNUIsRUFBdUVDLFFBQXZFLEVBQXFHO0FBQUEsUUFJakcsSUFBSSxDQUFDLElBQUtxYixDQUFBQSxhQUFMLENBQW1CWSxLQUFuQixDQUFELElBQ0EsQ0FBQyxJQUFLWixDQUFBQSxhQUFMLENBQW1CWSxLQUFuQixFQUEwQmxjLE1BQUEsQ0FBTzJCLElBQWpDLENBREQsSUFFQSxDQUFDLElBQUsyWixDQUFBQSxhQUFMLENBQW1CWSxLQUFuQixDQUFBLENBQTBCbGMsTUFBQSxDQUFPMkIsSUFBakMsQ0FBdUMzQixDQUFBQSxNQUFBLENBQU9QLE1BQTlDLENBRkwsRUFFNEQ7QUFBQSxZQUN4RCxPQUR3RDtBQUFBLFNBTnFDO0FBQUEsUUFVakcsTUFBTXNlLE1BQUEsR0FBUyxJQUFLekMsQ0FBQUEsYUFBTCxDQUFtQlksS0FBbkIsQ0FBQSxDQUEwQmxjLE1BQUEsQ0FBTzJCLElBQWpDLENBQXVDM0IsQ0FBQUEsTUFBQSxDQUFPUCxNQUE5QyxDQUFmLENBVmlHO0FBQUEsUUFXakcsT0FBTyxJQUFBLENBQUs2YixhQUFMLENBQW1CWSxLQUFuQixDQUEwQmxjLENBQUFBLE1BQUEsQ0FBTzJCLElBQWpDLENBQXVDM0IsQ0FBQUEsTUFBQSxDQUFPUCxNQUE5QyxDQUFQLENBWGlHO0FBQUEsUUFhakcsSUFBSXNlLE1BQUEsQ0FBT0QsWUFBUCxLQUF3QjlnQixTQUE1QixFQUF1QztBQUFBLFlBQ25DK2dCLE1BQUEsQ0FBT0QsWUFBUCxDQUFvQjlkLE1BQXBCLEVBQTRCQyxRQUE1QixDQURtQyxDQUFBO0FBQUEsU0FBdkMsTUFFTztBQUFBLFlBQ0hBLFFBQUEsRUFERyxDQUFBO0FBQUEsU0FmMEY7QUFBQSxLQTlJN0U7QUFBQTtBQXdLeEI7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBLElBQUErZCxnQkFBQSxDQUFpQjFlLEdBQWpCLEVBQThCVSxNQUE5QixFQUF1REMsUUFBdkQsRUFBaUY7QUFBQSxRQUM3RSxJQUFJO0FBQUEsWUFDQSxJQUFBLENBQUs2YSxJQUFMLENBQVVtRCxhQUFWLENBQXdCamUsTUFBQSxDQUFPa2UsR0FBL0IsQ0FEQSxDQUFBO0FBQUEsWUFFQWplLFFBQUEsRUFGQSxDQUFBO0FBQUEsU0FBSixDQUdFLE9BQU84WixDQUFQLEVBQVU7QUFBQSxZQUNSOVosUUFBQSxDQUFTOFosQ0FBQSxDQUFFb0UsUUFBRixFQUFULENBRFEsQ0FBQTtBQUFBLFNBSmlFO0FBQUEsS0F4S3pEO0FBQUEsSUFpTHhCQyxrQkFBQSxDQUFtQjllLEdBQW5CLEVBQWdDK2UsS0FBaEMsRUFBb0RwZSxRQUFwRCxFQUFpRjtBQUFBLFFBQzdFLElBQUk7QUFBQSxZQUNBMmIscUJBQUEsQ0FBb0IwQyxRQUFwQixDQUE2QkQsS0FBN0IsQ0FEQSxDQUFBO0FBQUEsWUFFQSxNQUFNRSxTQUFBLEdBQVkzQyxxQkFBQSxDQUFvQjRDLFlBQXBCLEVBQWxCLENBRkE7QUFBQSxZQUdBLElBQ0k1QyxxQkFBQSxDQUFvQjZDLFFBQXBCLEVBQ0EsSUFBQSxDQUFDN0MscUJBQUEsQ0FBb0JDLFFBQXBCLEVBREQsSUFFQTBDLFNBQUEsSUFBYSxJQUFBO0FBSGpCLEVBSUU7QUFBQSxnQkFDRSxLQUFLekQsSUFBTCxDQUFVbUQsYUFBVixDQUF3Qk0sU0FBeEIsQ0FERixDQUFBO0FBQUEsZ0JBRUUsTUFBTUcsUUFBQSxHQUFXOUMscUJBQUEsQ0FBb0JDLFFBQXBCLEVBQWpCLENBRkY7QUFBQSxnQkFHRSxNQUFNOEMsS0FBQSxHQUFRRCxRQUFBLEdBQVcxaEIsU0FBWCxHQUF1QixJQUFJNlEsS0FBSixDQUFVLENBQUMsOENBQUQsR0FBaUQwUSxTQUFqRCxFQUFBLENBQVYsQ0FBckMsQ0FIRjtBQUFBLGdCQUlFdGUsUUFBQSxDQUFTMGUsS0FBVCxFQUFnQkQsUUFBaEIsQ0FKRixDQUFBO0FBQUEsYUFQRjtBQUFBLFNBQUosQ0FhRSxPQUFPM0UsQ0FBUCxFQUFVO0FBQUEsWUFDUjlaLFFBQUEsQ0FBUzhaLENBQUEsQ0FBRW9FLFFBQUYsRUFBVCxDQURRLENBQUE7QUFBQSxTQWRpRTtBQUFBLEtBakx6RDtBQUFBLElBb014QlMsa0JBQUEsQ0FBbUIxQyxLQUFuQixFQUFpRDtBQUFBLFFBQzdDLElBQUk1RCxlQUFBLEdBQWtCLElBQUEsQ0FBS0EsZUFBTCxDQUFxQjRELEtBQXJCLENBQXRCLENBRDZDO0FBQUEsUUFHN0MsSUFBSSxDQUFDNUQsZUFBTCxFQUFzQjtBQUFBLFlBQ2xCQSxlQUFBLEdBQWtCLEVBQWxCLENBRGtCO0FBQUEsU0FIdUI7QUFBQSxRQU83QyxPQUFPQSxlQUFQLENBUDZDO0FBQUEsS0FwTXpCO0FBQUEsSUE4TXhCZ0YsYUFBQSxDQUFjcEIsS0FBZCxFQUE4QztBQUFBLFFBQzFDLElBQUlsQixZQUFBLEdBQWUsSUFBQSxDQUFLQSxZQUFMLENBQWtCa0IsS0FBbEIsQ0FBbkIsQ0FEMEM7QUFBQSxRQUUxQyxJQUFJLENBQUNsQixZQUFMLEVBQW1CO0FBQUEsWUFDZkEsWUFBQSxHQUFlLElBQUEsQ0FBS0EsWUFBTCxDQUFrQmtCLEtBQWxCLENBQTJCLEdBQUEsSUFBSTFkLGVBQUosRUFBMUMsQ0FEZTtBQUFBLFNBRnVCO0FBQUEsUUFLMUMsT0FBT3djLFlBQVAsQ0FMMEM7QUFBQSxLQTlNdEI7QUFBQSxJQXNOeEJ5QyxlQUFBLENBQWdCdkIsS0FBaEIsRUFBK0J2YSxJQUEvQixFQUE2Q2xDLE1BQTdDLEVBQTJFO0FBQUEsUUFDdkUsSUFBSSxDQUFDLElBQUEsQ0FBSzZiLGFBQUwsQ0FBbUJZLEtBQW5CLENBQUw7QUFBQSxZQUNJLEtBQUtaLGFBQUwsQ0FBbUJZLEtBQW5CLENBQUEsR0FBNEIsRUFBNUIsQ0FGbUU7QUFBQSxRQUd2RSxJQUFJLENBQUMsSUFBS1osQ0FBQUEsYUFBTCxDQUFtQlksS0FBbkIsQ0FBQSxDQUEwQnZhLElBQTFCLENBQUw7QUFBQSxZQUNJLElBQUEsQ0FBSzJaLGFBQUwsQ0FBbUJZLEtBQW5CLEVBQTBCdmEsSUFBMUIsQ0FBQSxHQUFrQyxFQUFsQyxDQUptRTtBQUFBLFFBTXZFLElBQUksQ0FBQyxJQUFBLENBQUsyWixhQUFMLENBQW1CWSxLQUFuQixDQUFBLENBQTBCdmEsSUFBMUIsQ0FBQSxDQUFnQ2xDLE1BQWhDLENBQUwsRUFBOEM7QUFBQTtBQUcxQztBQUFBLFlBQUEsTUFBTTJZLEtBQUEsR0FBUTtBQUFBLGdCQUNWeUcsSUFBQSxFQUFNLENBQUNsZCxJQUFELEVBQU9vUSxJQUFQLEVBQWE5UixRQUFiLEVBQXVCNmUsQ0FBdkIsRUFBMEJDLFNBQTFCLEVBQXFDQyxRQUFyQyxLQUFrRDtBQUFBLG9CQUNwRCxLQUFLNUcsS0FBTCxDQUFXeUcsSUFBWCxDQUFnQmxkLElBQWhCLEVBQXNCb1EsSUFBdEIsRUFBNEI5UixRQUE1QixFQUFzQ2ljLEtBQXRDLEVBQTZDNkMsU0FBN0MsRUFBd0RDLFFBQXhELENBRG9ELENBQUE7QUFBQSxpQkFEOUM7QUFBQSxnQkFJVkMsU0FBQSxFQUFXLElBQUs3RyxDQUFBQSxLQUFMLENBQVc2RyxTQUpaO0FBQUEsYUFBZCxDQUgwQztBQUFBLFlBUzFDLEtBQUszRCxhQUFMLENBQW1CWSxLQUFuQixDQUFBLENBQTBCdmEsSUFBMUIsQ0FBZ0NsQyxDQUFBQSxNQUFoQyxDQUEwQyxHQUFBLElBQUssS0FBSzJiLGlCQUFMLENBQXVCelosSUFBdkIsQ0FBTCxDQUF5Q3lXLEtBQXpDLEVBQXNELElBQUtrRixDQUFBQSxhQUFMLENBQW1CcEIsS0FBbkIsQ0FBdEQsRUFBaUYsSUFBQSxDQUFLMEMsa0JBQUwsQ0FBd0IxQyxLQUF4QixDQUFqRixFQUFpSCxLQUFLM0QsY0FBTCxDQUFvQjJELEtBQXBCLENBQWpILENBQTFDLENBVDBDO0FBQUEsU0FOeUI7QUFBQSxRQWtCdkUsT0FBTyxJQUFLWixDQUFBQSxhQUFMLENBQW1CWSxLQUFuQixFQUEwQnZhLElBQTFCLENBQUEsQ0FBZ0NsQyxNQUFoQyxDQUFQLENBbEJ1RTtBQUFBLEtBdE5uRDtBQUFBLElBMk94QmtlLGtCQUFBLENBQW1CekIsS0FBbkIsRUFBa0N6YyxNQUFsQyxFQUE2RTtBQUFBLFFBQ3pFLElBQUksQ0FBQyxJQUFBLENBQUs4YixnQkFBTCxDQUFzQlcsS0FBdEIsQ0FBTDtBQUFBLFlBQ0ksS0FBS1gsZ0JBQUwsQ0FBc0JXLEtBQXRCLENBQUEsR0FBK0IsRUFBL0IsQ0FGcUU7QUFBQSxRQUl6RSxJQUFJLENBQUMsSUFBS1gsQ0FBQUEsZ0JBQUwsQ0FBc0JXLEtBQXRCLENBQUEsQ0FBNkJ6YyxNQUE3QixDQUFMLEVBQTJDO0FBQUEsWUFDdkMsSUFBQSxDQUFLOGIsZ0JBQUwsQ0FBc0JXLEtBQXRCLENBQUEsQ0FBNkJ6YyxNQUE3QixDQUF1QyxHQUFBLElBQUlLLHlCQUFKLEVBQXZDLENBRHVDO0FBQUEsU0FKOEI7QUFBQSxRQVF6RSxPQUFPLEtBQUt5YixnQkFBTCxDQUFzQlcsS0FBdEIsQ0FBNkJ6YyxDQUFBQSxNQUE3QixDQUFQLENBUnlFO0FBQUEsS0EzT3JEO0FBQUEsSUFzUHhCeWYscUJBQUEsQ0FBc0JoRCxLQUF0QixFQUFxQy9OLEtBQXJDLEVBQW9EO0FBQUEsUUFDaEQrUSxvQ0FBQSxDQUFzQi9RLEtBQXRCLENBRGdELENBQUE7QUFBQSxLQXRQNUI7QUFBQSxJQTBQeEJnUiwyQkFBQSxDQUE0QmpELEtBQTVCLEVBQTJDbGMsTUFBM0MsRUFBd0RDLFFBQXhELEVBQTRHO0FBQUEsUUFDeEdBLFFBQUEsQ0FBU2pELFNBQVQsRUFBb0IsS0FBQSxDQUFwQixDQUR3RyxDQUFBO0FBQUEsS0ExUHBGO0FBQUEsQ0FsQzVCO0FBa1NBO0FBQUEsSUFBSSxPQUFPb2lCLGlCQUFQLEtBQTZCLFdBQTdCLElBQ0EsT0FBT3RFLElBQVAsS0FBZ0IsV0FEaEIsSUFFQUEsSUFBQSxZQUFnQnNFLGlCQUZwQixFQUV1QztBQUFBLElBQ25DdEUsSUFBQSxDQUFLaUQsTUFBTCxHQUFjLElBQUlsRCxNQUFKLENBQVdDLElBQVgsQ0FBZCxDQURtQztBQUFBOzs7Ozs7OzsifQ==
