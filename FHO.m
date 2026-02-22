classdef FHO < handle
    properties (Access=private)
        m double; % mass
        b double; % damping constant
        k double; % spring constant

        y0 double; % initial position
        v0 double; % initial velocity

        F function_handle;
    end
    properties (Access=public)
        ti double = 0
        tf double = 10
        dt double = 0.01;

        yt (1,:) double;
        T (1,:) double;
    end

    properties (Access=public)
        r1 double;
        r2 double;
        C1 double;
        C2 double;
    end

    methods (Access=public)
        function obj = FHO(m,b,k,F, y0,v0, ti,tf,dt)
            arguments
                m double; % mass
                b double; % damping constant
                k double; % spring constant
                F function_handle = @(t) 0;
        
                y0 double = 0; % initial position
                v0 double = 0; % initial velocity

                ti double = 0
                tf double = 10
                dt double = 0.01;
            end

            obj.m = m;
            obj.b = b;
            obj.k = k;
            obj.y0 = y0;
            obj.v0 = v0;
            obj.F = F;

            obj.ti = ti;
            obj.tf = tf;
            obj.dt = dt;

            obj.calcSolution();
        end

        function yt = getMotion(obj)
            yt = obj.yt;
        end

        function graph(obj)
            plot(obj.yt);
        end
    end
    methods (Access=private)
        function calcSolution(obj)
            % my'' + by' + ky - F = 0

            % y' = v
            % mv' + bv + ky - F(t) = 0
            % v' = (1/m)*(F(t)-bv-ky)
            
            obj.T = obj.ti:obj.dt:obj.tf;
            y = zeros(1,length(obj.T)); v = zeros(1,length(obj.T));
            y(1) = obj.y0; v(1) = obj.v0;
            for i = 1:length(y)-1
                y(i+1) = y(i) + obj.dt*v(i);
                v(i+1) = v(i) + obj.dt*(1/obj.m)*(obj.F(obj.T(i)) - obj.b*v(i) - obj.k*y(i)); %
            end

            obj.yt = y;
        end
    end
end
