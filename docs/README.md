# How to compile documentation
```shell script
cd docs/api
pip install -r requirements.txt
pip install -r ../../requirements.txt
make html
```
And then go to _build/html and open `index.html` in your browser. Please do not commit compiled doc files.