class MouseDragHelper {

    /**
     * @param {!Document} document The document to attach to.
     * @param {!EventTarget} target The element where dragging is started.
     */
    constructor(document, target) {
        this.document = document;
        this.target = target;
        this.targetRect = target.getBoundingClientRect(); //FIXME: Should we update this in onMouseDown?
        this.listeners = [];

        target.addEventListener("mousedown", event => {this.onMouseDown(event)});
        document.addEventListener("mousemove", event => {this.onMouseMove(event)});
        document.addEventListener("mouseup", event => {this.onMouseUp(event)});
    }

    onMouseDown(event) {
        this.listeners.forEach(element => {
            element.isActive = event.buttons === element.buttons;
            if (element.isDragging) {
                element.isDragging = false;
                element.onStopDragging && element.onStopDragging(this.getRelativeCoordinates(event));
            }
            if (element.isActive) {
                element.mouseDownLocation = this.getRelativeCoordinates(event);
                element.onMouseDown && element.onMouseDown(element.mouseDownLocation);
            }
        });
    }

    onMouseMove(event) {
        this.listeners.forEach(element => {
            if (element.isActive) {
                if (!element.isDragging) {
                    element.isDragging = true;
                    element.onStartDragging && element.onStartDragging(element.mouseDownLocation);
                }
                element.onDrag && element.onDrag(this.getRelativeCoordinates(event));
            }
        });
    }

    onMouseUp(event) {
        this.listeners.forEach(element => {
            if (event.buttons !== element.buttons) element.isActive = false;
            if (element.isDragging) {
                element.isDragging = false;
                element.onStopDragging && element.onStopDragging(this.getRelativeCoordinates(event));
            }
        });
    }

    getRelativeCoordinates(event) {
        return [Math.floor(event.clientX - this.targetRect.left) + 0.5, Math.floor(event.clientY - this.targetRect.top) + 0.5];
    }

    addMouseDragListener(listener) {
        if (typeof listener.buttons === "undefined") {
            listener.buttons = 1;
        }
        this.listeners.push(listener);
    }

}
