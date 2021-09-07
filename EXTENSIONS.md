# homelette Extensions

Extensions are homology modelling building blocks (model generation, model evaluation) that are developed by users and expand the homelette interface. homelette can and should be extended by custom Routines and Evaluations. We strongly encourage users to share extensions they themselves found useful with the community.

## Using extensions

Extensions are placed in the extension folder in the homelette package. The extension folder on your device can be found in the following way:

```python
import homelette.extension as ext
print(ext.__file__)
```

After an extension has been placed in the extension folder, it can be used as such:

```python
import homelette.extension.your_extension as ext_1
```

## Submitting an extension

Please contact us with a Pull Request on GitHub or via Email (philipp.junk@ucdconnect.ie) if you want to share your extension! Please make sure your extension is sufficiently annotated for others to use, in particular mentioning dependencies or other requirements.

## More information

More informations on extensions can be found in Tutorial 4 and in our online [documentation](https://homelette.readthedocs.io).
