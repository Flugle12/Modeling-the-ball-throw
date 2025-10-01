import math
import matplotlib.pyplot as plt

# Входные данные
m = 0.6  # кг
R_m = 0.12  # м
x0, y0 = 0, 0  # м
V0 = 7.8  # м/с
alpha = math.radians(45)  # градусы в радианы
x_k, y_k = 6.25, 0  # м
R_k = 0.225  # м
g = 9.8  # м/с²


# Функции движения
def x(t):
    return x0 + V0 * math.cos(alpha) * t


def y(t):
    return y0 + V0 * math.sin(alpha) * t - 0.5 * g * t ** 2


# Определение времени полета до высоты корзины
try:
    discriminant = (V0 * math.sin(alpha)) ** 2 - 2 * g * (y0 - y_k)
    if discriminant < 0:
        t_k = 2 * V0 * math.sin(alpha) / g  # полное время полета до земли
        print("Мяч не достигает высоты корзины, используем полное время полета")
    else:
        t_k = (V0 * math.sin(alpha) + math.sqrt(discriminant)) / g
except:
    t_k = 2 * V0 * math.sin(alpha) / g  # полное время полета до земли


# Поиск оптимального dt с проверкой попадания в корзину
def simulate(dt):
    t = 0
    while t <= t_k:
        x_curr, y_curr = x(t), y(t)
        # Проверяем, находится ли мяч в корзине
        distance_to_basket = math.sqrt((x_curr - x_k) ** 2 + (y_curr - y_k) ** 2)
        if distance_to_basket <= (R_k - R_m):
            return True
        t += dt
    return False


# Находим минимальный dt, при котором мяч попадает в корзину
dt = 0.1
eps = 0.000001  # Требуемая точность
hit_in_basket = False

try:
    # Сначала убедимся, что мяч вообще попадает в корзину
    if not simulate(dt):
        # Если не попадаем с dt=0.1, увеличиваем точность
        while not simulate(dt) and dt > eps:
            dt /= 2
            if dt < eps:
                hit_in_basket = False
                dt = 0.001  # используем маленький шаг для графика
                break
    else:
        # Если попадаем, ищем минимальный dt, при котором все еще попадаем
        hit_in_basket = True
        prev_dt = dt
        new_dt = dt / 2
        while abs(new_dt - dt) < eps:
            # Пробуем уменьшить dt в 2 раза
            new_dt = dt / 2

            # Проверяем, достигли ли требуемой точности
            # if abs(new_dt - dt) < eps:
            #     break

            if simulate(new_dt):
                # Если с новым dt все еще попадаем, продолжаем уменьшать
                prev_dt = dt
                dt = new_dt
            else:
                # Если с новым dt уже не попадаем, используем предыдущий
                dt = prev_dt
                break
except:
    hit_in_basket = False
    #dt = 0.001  # используем маленький шаг для графика

print(f"Оптимальный шаг dt: {dt:.6f}")
print(f"Попадание в корзину: {hit_in_basket}")

# Построение графика в любом случае
t_max = 2 * V0 * math.sin(alpha) / g  # полное время полета до земли
dt_plot = dt
t_values = []
x_values = []
y_values = []

t = 0
while t <= t_max:
    x_curr, y_curr = x(t), y(t)
    t_values.append(t)
    x_values.append(x_curr)
    y_values.append(y_curr)

    # Останавливаемся, если мяч упал на землю
    if y_curr < 0:
        break

    t += dt_plot

plt.figure(figsize=(10, 6))
plt.plot(x_values, y_values, label="Траектория мяча")
plt.scatter([x_k], [y_k], color="red", label="Центр корзины", zorder=5)

# Рисуем корзину
basket_circle = plt.Circle((x_k, y_k), R_k, color='red', fill=False, linestyle='--', label='Границы корзины')
plt.gca().add_patch(basket_circle)

# Рисуем мяч в точке попадания только если попадает
if hit_in_basket:
    t_hit = 0
    hit_found = False
    while t_hit <= t_k and not hit_found:
        x_curr, y_curr = x(t_hit), y(t_hit)
        distance_to_basket = math.sqrt((x_curr - x_k) ** 2 + (y_curr - y_k) ** 2)
        if distance_to_basket <= (R_k - R_m):
            plt.scatter([x_curr], [y_curr], color='green', s=100, zorder=10, label='Положение мяча при попадании')
            hit_found = True
            print(f"Мяч попадает в корзину при t = {t_hit:.4f} с, x = {x_curr:.3f} м")
        t_hit += dt_plot

# Определяем границы графика
x_max = max(max(x_values) if x_values else 0, x_k + R_k + 1)
y_max = max(max(y_values) if y_values else 0, y_k + R_k + 0.5)

plt.xlim(0, x_max)
plt.ylim(0, y_max)
plt.xlabel("Расстояние (м)")
plt.ylabel("Высота (м)")
plt.title("Траектория мяча и корзина")
plt.legend()
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

# Дополнительная информация
print(f"Время полета до корзины: {t_k:.3f} с")
print(f"Полное время полета: {2 * V0 * math.sin(alpha) / g:.3f} с")
if x_values and y_values:
    print(f"Максимальная высота: {max(y_values):.3f} м")
    print(f"Дальность полета: {x_values[-1]:.3f} м")
print(f"Координаты корзины: ({x_k:.3f}, {y_k:.3f})")
print(f"Радиус корзины: {R_k} м, радиус мяча: {R_m} м")
print(f"Допустимое расстояние для попадания: {R_k - R_m:.3f} м")
