classdef SHO < handle
    properties (Access=private)
        m double; % mass
        b double; % damping constant
        k double; % spring constant

        y0 double; % initial position
        v0 double; % initial velocity

    end
    properties (Access=public)
        yh function_handle;
        yt function_handle;
    end

    properties (Access=public)
        r1 double;
        r2 double;
        C1 double;
        C2 double;
    end

    methods (Access=public)
        function obj = SHO(m,b,k,y0,v0)
            arguments
                m double; % mass
                b double; % damping constant
                k double; % spring constant
        
                y0 double = 0; % initial position
                v0 double = 0; % initial velocity
            end

            obj.m = m;
            obj.b = b;
            obj.k = k;
            obj.y0 = y0;
            obj.v0 = v0;

            obj.calcHomogenous();
            obj.yt = @(t) obj.yh(t);
        end

        function yt = getMotion(obj)
            yt = obj.yt;
        end

        function graph(obj,a,b)
            xx = linspace(a,b,1000);
            yy = obj.yt(xx);

            plot(xx,yy); axis padded;
            xlim([a b]);

            if obj.b^2 < 4*obj.m*obj.k  % under damped
                hold on;
                R = ((obj.C1+obj.C2)^2 + (i*(obj.C1-obj.C2))^2)^.5;
                yy = R * exp(real(obj.r1)*xx);
                plot(xx,yy, 'r--');
            end
        end
    end
    methods (Access=private)
        function calcHomogenous(obj)


            obj.r1 = (-obj.b + sqrt(obj.b^2 - 4*obj.m*obj.k)) / (2*obj.m);
            obj.r2 = (-obj.b - sqrt(obj.b^2 - 4*obj.m*obj.k)) / (2*obj.m);

            disp("r1: " + obj.r1)
            disp("r2: " + obj.r2)

            if obj.r1 ~= obj.r2
                % yh = C1*e^(r1*t) + C2*e^(r2*t)
                % y(0) = C1+C2
                % y'(0) = r1*C1 + r2*C2

                A = [1 1; obj.r1 obj.r2];
                ic = [obj.y0; obj.v0];
                C = A\ic;
                obj.C1 = C(1); obj.C2 = C(2);
                obj.yh = @(t) obj.C1*exp(obj.r1*t) + obj.C2*exp(obj.r2*t);
                
            else
                % yh = C1*e^(r*t) + C2*t*e^(r*t)
                % y(0) = C1
                % y'(0) = r1*C1 + C2

                obj.C1 = obj.y0;
                obj.C2 = obj.v0 - obj.r1*obj.C1;
                obj.yh = @(t) obj.C1*exp(obj.r1*t) + obj.C2*t.*exp(obj.r2*t);

            end
        end
    end
end
