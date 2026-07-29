function newton1d_interactive(f, dfdx, x0_initial, N, tol)

    % Create figure
    fig = figure('Name','Newton 1D Interactive','NumberTitle','off');

    ax = axes('Parent',fig,'Position',[0.08 0.20 0.88 0.75]);
    
    % Slider for x0
    slider = uicontrol('Style','slider',...
        'Min',-5,'Max',5,...
        'Value',x0_initial,...
        'Units','normalized',...
        'Position',[0.2 0.02 0.6 0.05],...
        'Callback',@updatePlot);

    % Label to display x0
    label = uicontrol('Style','text',...
        'Units','normalized',...
        'Position',[0.42 0.08 0.2 0.04],...
        'String',['x0 = ' num2str(x0_initial)]);

    % Initial plot
    updatePlot(slider)

    function updatePlot(src,~)
        x0 = src.Value;
        label.String = ['x0 = ' num2str(x0,3)];

        cla(ax);
        hold(ax, 'on');

        % Run Newton
        [~,xs,ys,~,~] = newton1d(f, dfdx, x0, N, tol);

        % Plot
        plot_newton1d(f,xs,ys);

        hold off
    end

end


function [x,xs, ys, ds, converged] = newton1d(f, dfdx, x0, N, tol)
    x = x0;
    xs = zeros(1,N+1);
    ys = zeros(1,N+1);
    ds = zeros(1,N+1);

    n = 0;
    converged = false;

    z = 1e-12;

    % Invariant: 
    % - n = number of elements in arr = index of last computed x
    % - x = last computed x_n, x not in arr
    % - y and dfdx of last computed x is not computed before loop

    while n<N

        if(abs(dfdx(x))<z)
            error("dfdx = 0 at %g",x)
        end

        n = n+1;

        xs(n) = x;
        ys(n) = f(x);
        ds(n) = dfdx(x);

        x = xs(n) - ys(n)/ds(n);

        if(abs(x-xs(n))<tol*(1+abs(x)))
            converged = true;
            n=n+1;
            xs(n)=x;
            ys(n)=f(x);
            ds(n)=dfdx(x);
            break;
        end
    end

    
    xs = xs(1:n);
    ys = ys(1:n);
    ds = ds(1:n);

    plot_newton1d(f,xs,ys);

end


function plot_newton1d(f,xs,ys)
    %figure("Name","Newton1d","NumberTitle","off")
    %clf;

    ratio = 0.6;
    density = 200;

    xmin = min(xs) - ratio*abs(min(xs));
    xmax = max(xs) + ratio*abs(max(xs));
    xmin = min(xmin,-5);
    xmax = max(xmax,5);
    range = xmax-xmin;

    points = density * range;

    n = length(xs);

    %fig = figure;
    
    x = linspace(xmin,xmax,points);
    y = arrayfun(f,x);
    plot(x,y,"-b");
    title("Newton 1D");
    xlabel('x');
    ylabel('f(x)');
    grid on;
    hold on;
    
    plot(xs(n),ys(n),"or");
    ax = gca;            
    ax.XAxisLocation = 'origin'; 
    ax.YAxisLocation = 'origin'; 
    ax.Box = 'off';  

    for i = 1:(n-1)
        if i==1
            plot(xs(1),0,".r")
        else
            plot(xs(i),0,".k")
        end

        plot([xs(i) xs(i)],[0 ys(i)],"--k")
        plot(xs(i),ys(i),".k")
        plot([xs(i) xs(i+1)],[ys(i) 0], "-k")
    end

    grid off;
    hold off;
    
end