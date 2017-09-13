# CE3202-Project1-PolynomialSolve

En la comunidad de C++, la biblioteca boost ofrece gran cantidad de estructuras de datos
y algoritmos. De allı́, en esta tarea partirá de la clase boost::math::tools::polynomial,
que ofrece una forma sencilla de representar polinomios.
Dicha propuesta utiliza el método externo a la clase evaluate polynomial para evaluar
un valor dado en el polinomio. También ofrece una funcionalidad básica para realizar la
división polinomial.
En este proyecto usted ofrecerá un método para encontrar todas las raı́ces, reales y com-
plejas de polinomios que pueden tener como elementos objetos de los tipos float, double,
std::complex<float> o std::complex<double>, para lo cual deberá entonces implementar todo su código por medio de plantillas (C++ templates).

Pasos para ejecutar en Eclipse:

1. Ir a Project - Properties - C/C++ Build - Settings.
2. Escribir en command:
	g++ -std=c++11
3. Abrir terminal e instalar la biblioteca Boost:
	sudo apt-get install libboost-all-dev con Ubuntu 17

