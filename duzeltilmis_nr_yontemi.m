tolerans = input('Hata toleransını girin (örn: 1e-6): ');

f = @(x) x.^3 + x.^2 - 8*x - 12;
df = @(x) 3*x.^2 + 2*x - 8;
ddf = @(x) 6*x + 2;

% düzeltilmiş newton-raphson
function [root, iterations] = modified_newton(f, df, ddf, x0, tolerans, max_iter)
    iterations = [x0];
    for k = 1:max_iter
        x = iterations(end);
        fx = f(x);
        dfx = df(x);
        ddfx = ddf(x);
        pay = fx * dfx;
        payda = dfx^2 - fx * ddfx;
        if abs(payda) < 1e-12
            error('Payda sıfır!');
        end
        x_new = x - pay / payda;
        iterations = [iterations; x_new];
        if abs(f(x_new)) < tolerans
            break;
        end
    end
    root = x_new;
end

% standart newton-raphson
function [root, iterations] = standard_newton(f, df, x0, tolerans, max_iter)
    iterations = [x0];
    for k = 1:max_iter
        x = iterations(end);
        fx = f(x);
        dfx = df(x);
        if dfx == 0
            error('Türev sıfır!');
        end
        x_new = x - fx / dfx;
        iterations = [iterations; x_new];
        if abs(f(x_new)) < tolerans
            break;
        end
    end
    root = x_new;
end

% katlı kök (x0=0):
[root_mult, iter_mult] = modified_newton(f, df, ddf, 0, tolerans, 100);

% tek kök (x0=4):
[root_single, iter_single] = standard_newton(f, df, 4, tolerans, 100);

fprintf('Katlı kök: x = %.8f\n', root_mult);
fprintf('Tek kök: x = %.8f\n', root_single);

% grafikler:
figure;
plot(0:length(iter_mult)-1, iter_mult, '-o');
xlabel('İterasyon');
ylabel('x');
title('Katlı Köke Yakınsama');
grid on;

figure;
plot(0:length(iter_single)-1, iter_single, '-o');
xlabel('İterasyon');
ylabel('x');
title('Tek Köke Yakınsama');
grid on;
