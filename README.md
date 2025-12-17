# valores-iniciais
Biblioteca para geração de valores iniciais. É utilizada no [gravidade-fortran](https://github.com/potalej/gravidade-fortran).

# Binding para Python.

Através do [F2Py](https://numpy.org/doc/stable/f2py/), é possível compilar esta biblioteca e utilizá-la no python. Para isso, siga os passos:

1. Instale o NumPy, o F2PY é um de seus módulos:

```bash
python -m pip install numpy
```

2. Tenha o OpenBLAS/LAPACK instalado. Isso varia de sistema para sistema.

3. Compile:
```bash
cd f2py
cmake -B build
cd build
make
```
> **Dica**: se estiver usando o conda, talvez deva rodar `cmake -B build -DPython3_EXECUTABLE=$CONDA_PREFIX/bin/python` em vez de `cmake -B build`.

A biblioteca compilada será copiada automaticamente para o diretório "python/valores_iniciais". O script "python/exemplo.py" contém um exemplo simples de uso.

Por hora, para usar, copie essa pasta "valores_iniciais" para onde quiser e faça como no exemplo. Mais para frente transformo em uma biblioteca de verdade.