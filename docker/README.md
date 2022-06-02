To run the check in the main package folder:

```bash
docker run -v$(pwd):/pkg/ -w/pkg --rm -i uscbiostats/fdrci:latest make check
```

To update the version just type

```bash
make push
```
