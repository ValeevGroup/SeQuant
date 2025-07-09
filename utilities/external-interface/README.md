# External interface

The external interface can be thought of as the glue between SeQuant and any arbitrary other program. At the moment, the intended workflow is to
derive the underlying diagrams in an external program and then use SeQuant for post-processing and code generation.

The details of the processing is defined in a _driver_, which is just a [JSON](https://www.json.org/json-en.html) file. See the
[examples](./examples/) for more details.
