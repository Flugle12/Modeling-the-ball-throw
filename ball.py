import math
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('TkAgg')
# Входные данные
m = 0.6  # кг
R_m = 0.12  # м
x0, y0 = 0, 0  # м
V0 = 8.1  # м/с
alpha_deg = 45  # угол в градусах
alpha = math.radians(alpha_deg)  # угол в радианах
x_k, y_k = 6.25, 0  # м
R_k = 0.225  # м
g = 9.8  # м/с²


# Функции движения (аналитическое решение)
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


# Численное решение методом Эйлера для проверки попадания
def numerical_check_hit(dt):
    """
    Численная проверка попадания в корзину методом Эйлера

    Возвращает:
    success - флаг попадания в корзину
    t_hit - время попадания
    x_hit - координата x в момент попадания
    y_hit - координата y в момент попадания
    min_distance - минимальное расстояние до корзины
    """
    # Начальные условия
    t = 0
    x_curr = x0
    y_curr = y0
    Vx = V0 * math.cos(alpha)
    Vy = V0 * math.sin(alpha)

    min_distance = float('inf')
    closest_point = (0, 0, 0)  # (t, x, y)

    # Интегрирование методом Эйлера
    while y_curr >= y_k and t < 10:  # ограничение по времени
        # Сохраняем предыдущие значения
        x_prev = x_curr
        y_prev = y_curr
        t_prev = t

        # Шаг интегрирования
        x_curr = x_curr + Vx * dt
        y_curr = y_curr + Vy * dt
        Vy = Vy - g * dt
        t = t + dt

        # Вычисляем расстояние до корзины
        distance_to_basket = math.sqrt((x_curr - x_k) ** 2 + (y_curr - y_k) ** 2)

        # Запоминаем точку минимального расстояния
        if distance_to_basket < min_distance:
            min_distance = distance_to_basket
            closest_point = (t, x_curr, y_curr)

        # Проверяем пересечение уровня корзины и попадание
        if y_curr < y_k and y_prev >= y_k:
            # Линейная интерполяция для точного определения момента
            fraction = (y_k - y_prev) / (y_curr - y_prev)
            t_interp = t_prev + fraction * dt
            x_interp = x_prev + fraction * (x_curr - x_prev)
            y_interp = y_k

            # Проверяем, находится ли мяч в корзине
            distance_interp = math.sqrt((x_interp - x_k) ** 2 + (y_interp - y_k) ** 2)
            if distance_interp <= (R_k - R_m):
                return True, t_interp, x_interp, y_interp, min_distance

    # Если не достигли уровня корзины или не попали
    return False, closest_point[0], closest_point[1], closest_point[2], min_distance


# Аналитическая проверка попадания в корзину
def analytical_check_hit(dt):
    """
    Аналитическая проверка попадания в корзину

    Возвращает:
    success - флаг попадания в корзину
    t_hit - время попадания
    x_hit - координата x в момент попадания
    y_hit - координата y в момент попадания
    min_distance - минимальное расстояние до корзины
    """
    t = 0
    min_distance = float('inf')
    closest_point = (0, 0, 0)  # (t, x, y)

    while t <= t_k:
        x_curr, y_curr = x(t), y(t)
        # Вычисляем расстояние до корзины
        distance_to_basket = math.sqrt((x_curr - x_k) ** 2 + (y_curr - y_k) ** 2)

        # Запоминаем точку минимального расстояния
        if distance_to_basket < min_distance:
            min_distance = distance_to_basket
            closest_point = (t, x_curr, y_curr)

        # Проверяем, находится ли мяч в корзине
        if distance_to_basket <= (R_k - R_m):
            return True, t, x_curr, y_curr, min_distance

        t += dt

    return False, closest_point[0], closest_point[1], closest_point[2], min_distance


# Автоматический подбор оптимального шага dt для заданной точности (численный метод)
def find_optimal_dt_numerical(precision=0.0001, initial_dt=0.1, max_iterations=50):
    """
    Автоматический подбор оптимального шага dt для заданной точности с использованием численного метода

    Возвращает:
    optimal_dt - оптимальный шаг
    hit_in_basket - флаг попадания
    hit_time - время попадания
    hit_x - координата x попадания
    hit_y - координата y попадания
    min_distance - минимальное расстояние до корзины
    iterations - количество итераций
    """
    dt = initial_dt
    min_dt = 0.0001  # минимальный шаг
    iterations = 0

    print("\nПоиск оптимального шага dt (численный метод):")
    print("Итерация | Шаг dt (с)   | Попадание | Расстояние (м) | Погрешность (м)")
    print("---------|--------------|-----------|----------------|----------------")

    # Используем численное решение с очень малым шагом как "точное"
    hit_ref, t_ref, x_ref, y_ref, min_dist_ref = numerical_check_hit(min_dt)

    best_dt = dt
    best_hit = False
    best_t = 0
    best_x = 0
    best_y = 0
    best_min_dist = float('inf')

    while (x(dt) - x(dt/2)) > min_dt:
        iterations += 1

        # Выполняем проверку с текущим шагом (численный метод)
        hit, t_hit, x_hit, y_hit, min_dist = numerical_check_hit(dt)

        # Для оценки погрешности используем численное решение с малым шагом
        error = 0
        if hit and hit_ref:
            error = math.sqrt((x_hit - x_ref) ** 2 + (y_hit - y_ref) ** 2)
        elif not hit and not hit_ref:
            error = abs(min_dist - min_dist_ref)
        else:
            # Если тип результата разный, считаем погрешность максимальной
            error = float('inf')

        status = "Попадание" if hit else "Промах"
        print(f"{iterations:8d} | {dt:<12.6f} | {status:<9} | {min_dist:<14.6f} | {error:<14.6f}")

        # Сохраняем лучший результат
        if error < best_min_dist:
            best_dt = dt
            best_hit = hit
            best_t = t_hit
            best_x = x_hit
            best_y = y_hit
            best_min_dist = min_dist

        # Проверяем достижение требуемой точности
        # if error <= precision:
        #     print("\nТребуемая точность достигнута!")
        #     return dt, hit, t_hit, x_hit, y_hit, min_dist, iterations

        # Уменьшаем шаг для следующей итерации
        dt /= 2

    print("\nДостигнута максимальная точность с текущими параметрами")
    return best_dt, best_hit, best_t, best_x, best_y, best_min_dist, iterations


# Основная программа
if __name__ == "__main__":
    # Поиск оптимального dt численным методом
    dt, hit_in_basket, hit_time, hit_x, hit_y, min_distance, iterations = find_optimal_dt_numerical()

    print(f"\nРезультаты численного подбора dt:")
    print(f"Оптимальный шаг dt: {dt:.6f}")
    print(f"Количество итераций: {iterations}")
    print(f"Попадание в корзину: {hit_in_basket}")

    if not hit_in_basket:
        print(f"Минимальное расстояние до корзины: {min_distance:.4f} м")
        print(f"Допустимое расстояние для попадания: {R_k - R_m:.3f} м")
    else:
        print(f"Мяч попадает в корзину при t = {hit_time:.4f} с, x = {hit_x:.3f} м, y = {hit_y:.3f} м")

    # Проверяем результат аналитическим методом с найденным dt
    print("\nПроверка аналитическим методом с найденным dt:")
    hit_analytical, t_analytical, x_analytical, y_analytical, min_dist_analytical = analytical_check_hit(dt)

    status_analytical = "Попадание" if hit_analytical else "Промах"
    print(f"Аналитический метод: {status_analytical}")
    if hit_analytical:
        print(f"Время: {t_analytical:.4f} с, x = {x_analytical:.3f} м, y = {y_analytical:.3f} м")
    else:
        print(f"Минимальное расстояние: {min_dist_analytical:.4f} м")

    # Построение графика с использованием аналитического решения
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
    plt.plot(x_values, y_values, label="Траектория мяча (аналитическая)")
    plt.scatter([x_k], [y_k], color="red", label="Центр корзины", zorder=5)

    # Рисуем корзину
    basket_circle = plt.Circle((x_k, y_k), R_k, color='red', fill=False, linestyle='--', label='Границы корзины')
    plt.gca().add_patch(basket_circle)

    # Рисуем мяч в точке попадания только если попадает
    if hit_in_basket:
        plt.scatter([hit_x], [hit_y], color='green', s=100, zorder=10, label='Положение мяча при попадании (численное)')

    # Также показываем точку попадания по аналитическому методу
    if hit_analytical:
        plt.scatter([x_analytical], [y_analytical], color='blue', s=100, zorder=10,
                    label='Положение мяча при попадании (аналитическое)')
    else:
        # Показываем точку минимального расстояния
        hit, t_min, x_min, y_min, dist_min = analytical_check_hit(dt)
        plt.scatter([x_min], [y_min], color='orange', s=100, zorder=10,
                    label='Точка минимального расстояния (аналитическая)')

    # Определяем границы графика
    x_max_val = max(max(x_values) if x_values else 0, x_k + R_k + 1)
    y_max_val = max(max(y_values) if y_values else 0, y_k + R_k + 0.5)

    plt.xlim(0, x_max_val)
    plt.ylim(0, y_max_val)
    plt.xlabel("Расстояние (м)")
    plt.ylabel("Высота (м)")
    plt.title("Траектория мяча и корзина (аналитическое решение)")
    plt.legend()
    plt.grid(True)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()

    # Дополнительная информация
    print(f"\nДополнительная информация:")
    print(f"Время полета до корзины: {t_k:.3f} с")
    print(f"Полное время полета: {2 * V0 * math.sin(alpha) / g:.3f} с")
    if x_values and y_values:
        print(f"Максимальная высота: {max(y_values):.3f} м")
        print(f"Дальность полета: {x_values[-1]:.3f} м")
    print(f"Координаты корзины: ({x_k:.3f}, {y_k:.3f})")
    print(f"Радиус корзины: {R_k} м, радиус мяча: {R_m} м")
