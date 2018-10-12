% Algoritmo para hacer un ajuste lineal
% x = vector con valores x0, x1
% y = vector con valores y0, y1
% y_ast = es el valor constante
% x_ast = el valor a encontrar
% Dado una constante y*, y dos puntos (x0, y0) y (x1, y1)
% x* = (y*-y0)*(x1-x0)/(y1-y0) + x0

function x_ast = inter_extra_polation(x,y,y_ast)
% close all, hold on
% x = [1 3]; y = [2 4]; y_ast = 5;
x_ast = x(1) + (y_ast-y(1))*(x(2)-x(1))/(y(2)-y(1));
%%%%%% Graficar los
% plot(x,y, 'b-')
% plot(x_ast,y_ast, 'r*')
% plot([min([x,x_ast,0]) x_ast],[y_ast, max([0 y_ast])],'-b')
% plot([x_ast x_ast], [min([0, y,y_ast])  y_ast]  ,':r')
end