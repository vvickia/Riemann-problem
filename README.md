# Задача Римана о распаде произвольного разрыва

Программа реализует точное решение задачи распада произвольного газодинамического разрыва в идеальном газе с показателем адиабаты $\gamma = 5/3$.
Первоначально разрыв находится в начале координат.

Программа протестирована для трех начальных условий:
1. $\{\rho, v, p\}_L = \{1, 0, 3\}$,
   $\{\rho, v, p\}_R = \{1, 0, 1\}$;
2. $\{\rho, v, p\}_L = \{1, 1, 3\}$,
   $\{\rho, v, p\}_R = \{1, -1, 1\}$;
3. $\{\rho, v, p\}_L = \{1, -0.1, 1\}$,
   $\{\rho, v, p\}_R = \{1, 0.2, 1\}$,

здесь индексы $L$ и $R$ обозначают характеристики газа (плотность, скорость, давление) слева и справа от разрыва соответственно.

Структура течения проанализирована в моменты времени $t = 0.18$ в тесте 1 и $t = 0.1$ в тестах 2, 3.

## Запуск:

Из терминала:  
```bash
g++ main.cpp de_allocate.cpp iof.cpp configurations.cpp equations.cpp newton.cpp riemann.cpp -o main
./main
```

Чтобы построить графики (так же из терминала):
```bash
python3 graphs/graphs.py
```
