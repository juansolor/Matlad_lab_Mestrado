classdef ESP
    properties
        parameters;
        states;        outputs;        inputs;       simulation;
        stationary; symbolic; LinearModel; envelope
    end
    methods
        function obj = ESP()
            % Valores padrões para o posso
            obj.parameters.g   = 9.81;   % Gravitational acceleration constant [m/s²]
            obj.parameters.Cc = 2e-5 ;   % Choke valve constant
            obj.parameters.A1 = 0.008107;% Cross-section area of pipe below ESP [m²]
            obj.parameters.A2 = 0.008107;% Cross-section area of pipe above ESP [m²]
            obj.parameters.D1 = 0.1016;  % Pipe diameter below ESP [m]
            obj.parameters.D2 = 0.1016;  % Pipe diameter above ESP [m]
            obj.parameters.h1 = 200;     % Heigth from reservoir to ESP [m]
            obj.parameters.hw = 1000;    % Total vertical distance in well [m]
            obj.parameters.L1 =  500;    % Length from reservoir to ESP [m]
            obj.parameters.L2 = 1200;    % Length from ESP to choke [m]
            obj.parameters.V1 = 4.054;   % Pipe volume below ESP [m3]
            obj.parameters.V2 = 9.729;   % Pipe volume above ESP [m3]
            obj.parameters.f0 = 60;      % ESP characteristics reference freq [Hz]
            obj.parameters.q0_dt = 0.00542847756833871; % Downtrhust flow at f0 m3/s
            obj.parameters.q0_ut = 0.0134581949664923;  % Uptrhust flow at f0 m3/s
            obj.parameters.Inp = 65;     % ESP motor nominal current [A]
            obj.parameters.Pnp = 1.625e5;% ESP motor nominal Power [W]
            obj.parameters.b1 = 1.5e9;   % Bulk modulus below ESP [Pa]
            obj.parameters.b2 = 1.5e9;   % Bulk modulus above ESP [Pa]
            obj.parameters.M  = 1.992e8; % Fluid inertia parameters [kg/m4]
            obj.parameters.rho = 950;    % Density of produced fluid [kg/m³]
            obj.parameters.PI = 2.32e-9; % Well productivy index [m3/s/Pa]
            obj.parameters.mu  = 0.025;  % Viscosity [Pa*s]
            obj.parameters.dfq_max = 0.5;   % máxima variação em f/s
            obj.parameters.dzc_max = 1;   % máxima variação em zc %/s
            obj.parameters.tp =[1/obj.parameters.dfq_max;1/obj.parameters.dzc_max];  % Actuator Response time
            obj.parameters.CH = -0.03*obj.parameters.mu + 1;
            obj.parameters.Cq = 2.7944*obj.parameters.mu^4 - 6.8104*obj.parameters.mu^3 + 6.0032*obj.parameters.mu^2 - 2.6266*obj.parameters.mu + 1;
            obj.parameters.Cp = -4.4376*obj.parameters.mu^4 + 11.091*obj.parameters.mu^3 -9.9306*obj.parameters.mu^2 + 3.9042*obj.parameters.mu + 1;
                    
            
            % Valores nominais para as entradas
            obj.inputs.fqref = 50; % Rotational frequency
            obj.inputs.zcref = 50; % Chocke opening
            obj.inputs.pr = 1.26e7; % Reservoir pressure
            obj.inputs.pm = 2e5;%101325;%1atm fez o foco do sinal ficar mais centrado no envelope   nominal: 2e6;    % Manifold pressure
            obj.inputs.u = [obj.inputs.fqref;obj.inputs.zcref;obj.inputs.pm;obj.inputs.pr];
            
            % Calculo dos valores nominais das saidas
            obj.states.x = [8311024.82175957;2990109.06207437;0.00995042241351780;50;50];
            obj = obj.EvaluateModel(obj.states.x,obj.inputs.u);
            
            % Estacionário
            stationary = obj.StateSolve(obj.states.x,obj.inputs.u);
            obj.stationary = stationary;
            % Construção do modelo simbolic
            obj = obj.SymbolicModel();
            
            % Construção do modelo linear para o ponto
            obj = obj.Linearization(obj.stationary.xss,obj.stationary.uss);
            
            % Envelope
            %obj = obj.Envelope();
            obj = obj.OperationalEnvelop();
            
            fprintf('Modelo BCS carregado para o ponto de operação\n')
            fprintf('zc = %.2f%% e fq = %.2fHz \n',obj.inputs.fqref,obj.inputs.zcref);
        end
        
        function obj = EvaluateModel(obj,x,u)
            % Estados
            pbh = x(1); pwh = x(2); q = x(3); fq = x(4); zc = x(5);
            % Entradas
            fqref = u(1); zcref = u(2); pm = u(3); pr = u(4);
            
            % SEA
            % Calculo do HEAD e delta de pressão
            q0 = q/obj.parameters.Cq*(obj.parameters.f0/fq);
            H0 = -1.2454e6*q0^2 + 7.4959e3*q0 + 9.5970e2;
            H = obj.parameters.CH*H0*(fq/obj.parameters.f0)^2; % Head
            Dp = obj.parameters.rho*obj.parameters.g*H;       % Delta de pressão
            % Calculo da Potencia e corrente da bomba
            P0 = -2.3599e9*q0^3 -1.8082e7*q0^2 +4.3346e6*q0 + 9.4355e4;
            P = obj.parameters.Cp*P0*(fq/obj.parameters.f0)^3; % Potencia
            I = obj.parameters.Inp*P/obj.parameters.Pnp;       % Corrente
            % Calculo da pressão de intaike
            F1 = 0.158*((obj.parameters.rho*obj.parameters.L1*q^2)/(obj.parameters.D1*obj.parameters.A1^2))*(obj.parameters.mu/(obj.parameters.rho*obj.parameters.D1*q))^(1/4);
            F2 = 0.158*((obj.parameters.rho*obj.parameters.L2*q^2)/(obj.parameters.D2*obj.parameters.A2^2))*(obj.parameters.mu/(obj.parameters.rho*obj.parameters.D2*q))^(1/4);
            pin = pbh - obj.parameters.rho*obj.parameters.g*obj.parameters.h1 - F1;
            % Vazao do rezervatorio vazao da chocke
            qr  = obj.parameters.PI*(pr - pbh);
            qc  = obj.parameters.Cc*(zc/100)*sign((pwh - pm))*sqrt(abs(pwh - pm));
                   
            % System of nonlinear ordinary diferenctial equations
            dpbhdt = obj.parameters.b1/obj.parameters.V1*(qr - q);
            dpwhdt = obj.parameters.b2/obj.parameters.V2*(q - qc);
            dqdt = 1/obj.parameters.M*(pbh - pwh - obj.parameters.rho*obj.parameters.g*obj.parameters.hw - F1 - F2 + Dp);
            dfqdt = (fqref - fq)/obj.parameters.tp(1);
            if abs(dfqdt) > obj.parameters.dfq_max
                dfqdt = obj.parameters.dfq_max*sign(dfqdt);
            end
            dzcdt = (zcref - zc)/obj.parameters.tp(2);
            if abs(dzcdt) > obj.parameters.dzc_max
                dzcdt = obj.parameters.dzc_max*sign(dzcdt);
            end 
            
            % States
            dxdt = [dpbhdt;dpwhdt;dqdt;dfqdt;dzcdt];
            obj.states = struct('pbh',pbh,'pwh',pwh,'q',q,'fq',fq,'zc',zc,'x',x,'dxdt',dxdt);
            
            % inputs
            obj.inputs = struct('fqref',fqref,'zcref',zcref,'pm',pm,'pr',pr,'u',u);
            
            % Outputs
            y = [pin,H,qc,qr,P,I,Dp]';
            obj.outputs = struct('pin',pin,'H',H,'qc',qc,'qr',qr,'P',P,'I',I,'Dp',Dp,'y',y);
        end
        
        function [y,x,sim,dxdt] = SimModel(obj,xk,uk_1,tsim)
            [~,xsim] = ode15s(@(t,x) EqEstado(obj,x,uk_1),tsim,xk);
            sim.x = xsim';
            sim.u = repmat(uk_1,1,size(xsim,1));
            sim.y = zeros(length(obj.outputs.y),size(xsim,1));
            for ix = 1:size(xsim,1)
                obj = obj.EvaluateModel(xsim(ix,:)',uk_1);
                sim.y(:,ix) = obj.outputs.y;
            end
            y = obj.outputs.y;
            x = obj.states.x;
            dxdt = EqEstado(obj,x,uk_1);
        end
        
        function stationary = StateSolve(obj,xk,uk_1)
            % obj = obj.SimModel(xk,uk_1,[0,100]);
            [xss,fval,ef] = fsolve(@(x) obj.EqEstado(x,uk_1),[xk(1:3);uk_1(1:2)],optimoptions(@fsolve,'Display','off'));%xk ou [xk(1:3);uk_1(1:2)] como condiçao inicial
            obj = obj.EvaluateModel(xss,uk_1);
            
            stationary = struct('xss',xss,'yss',obj.outputs.y,'uss',obj.inputs.u,'opt',struct('fval',fval,'ef',ef));
        end
        
        function obj = Linearization(obj,xss,uss)
            A = obj.LinearModel.A(xss,uss);
            B = obj.LinearModel.B(xss,uss);
            C = obj.LinearModel.C(xss,uss);
            obj.LinearModel.SS = ss(A,B,C,[]);
            obj.LinearModel.TF = tf(obj.LinearModel.SS);
            
            An = obj.LinearModel.Normalized.A(xss,uss);
            Bn = obj.LinearModel.Normalized.B(xss,uss);
            Cn = obj.LinearModel.Normalized.C(xss,uss);
            obj.LinearModel.Normalized.SS = ss(An,Bn,Cn,[]);
            obj.LinearModel.Normalized.TF = tf(obj.LinearModel.Normalized.SS);
        end
        
        % Private
        
        function dxdt = EqEstado(obj,xk,uk_1)
            obj = obj.EvaluateModel(xk,uk_1);
            dxdt = obj.states.dxdt;
        end
        
        function obj = SymbolicModel(obj)
            x = sym('x',[5,1]); u = sym('u',[4,1]);            
           
            % Estados
            pbh = x(1); pwh = x(2); q = x(3); fq = x(4); zc = x(5);
            % Entradas
            fqref = u(1); zcref = u(2); pm = u(3); pr = u(4);
            
            % SEA
            % Calculo do HEAD e delta de pressão
            q0 = q/obj.parameters.Cq*(obj.parameters.f0/fq);
            H0 = -1.2454e6*q0^2 + 7.4959e3*q0 + 9.5970e2;
            H = obj.parameters.CH*H0*(fq/obj.parameters.f0)^2; % Head
            Dp = obj.parameters.rho*obj.parameters.g*H;       % Delta de pressão
            
            % Calculo da Potencia e corrente da bomba
            P0 = -2.3599e9*q0^3 -1.8082e7*q0^2 +4.3346e6*q0 + 9.4355e4;
            P = obj.parameters.Cp*P0*(fq/obj.parameters.f0)^3; % Potencia
            I = obj.parameters.Inp*P/obj.parameters.Pnp;       % Corrente
            
            % Calculo da pressão de intaike
            F1 = 0.158*((obj.parameters.rho*obj.parameters.L1*q^2)/(obj.parameters.D1*obj.parameters.A1^2))*(obj.parameters.mu/(obj.parameters.rho*obj.parameters.D1*q))^(1/4);
            F2 = 0.158*((obj.parameters.rho*obj.parameters.L2*q^2)/(obj.parameters.D2*obj.parameters.A2^2))*(obj.parameters.mu/(obj.parameters.rho*obj.parameters.D2*q))^(1/4);
            pin = pbh - obj.parameters.rho*obj.parameters.g*obj.parameters.h1 - F1;
            
            % Vazao do rezervatorio vazao da chocke
            qr  = obj.parameters.PI*(pr - pbh);
            qc  = obj.parameters.Cc*(zc/100)*sign((pwh - pm))*sqrt(abs(pwh - pm));
                   
            % System of nonlinear ordinary diferenctial equations
            dpbhdt = obj.parameters.b1/obj.parameters.V1*(qr - q);
            dpwhdt = obj.parameters.b2/obj.parameters.V2*(q - qc);
            dqdt = 1/obj.parameters.M*(pbh - pwh - obj.parameters.rho*obj.parameters.g*obj.parameters.hw - F1 - F2 + Dp);
            dfqdt = (fqref - fq)/obj.parameters.tp(1);
            dzcdt = (zcref - zc)/obj.parameters.tp(2);
            
            % States
            obj.symbolic.x = x;
            
            % inputs
            obj.symbolic.u = u;
            
            % Nonlinear Model
            obj.symbolic.NonLinear.dxdt = [dpbhdt;dpwhdt;dqdt;dfqdt;dzcdt];
            obj.symbolic.NonLinear.y = [pin;H;qc;qr;P;I;Dp];
            
            % Linear model
            A = symfun(jacobian(obj.symbolic.NonLinear.dxdt,obj.symbolic.x),[x;u]);
            B = symfun(jacobian(obj.symbolic.NonLinear.dxdt,obj.symbolic.u),[x;u]);
            C = symfun(jacobian(obj.symbolic.NonLinear.y,obj.symbolic.x),[x;u]);
            
            obj.LinearModel.A = @(xss,uss) double(A(xss(1),xss(2),xss(3),xss(4),xss(5),uss(1),uss(2),uss(3),uss(4)));
            obj.LinearModel.B = @(xss,uss) double(B(xss(1),xss(2),xss(3),xss(4),xss(5),uss(1),uss(2),uss(3),uss(4)));
            obj.LinearModel.C = @(xss,uss) double(C(xss(1),xss(2),xss(3),xss(4),xss(5),uss(1),uss(2),uss(3),uss(4)));
            
            % Normalized -1 to 1
            An = sym('An',[length(x),length(x)]);
            Bn = sym('Bn',[length(x),length(u)]);
            Cn = sym('Cn',[length(obj.symbolic.NonLinear.y),length(x)]);
            for ix = 1:length(x)
                An(ix,:) = transpose(x)/x(ix);
                Bn(ix,:) = transpose(u)/x(ix);
            end
            for iy = 1:length(obj.symbolic.NonLinear.y)
                Cn(iy,:) = transpose(x)/obj.symbolic.NonLinear.y(iy);
            end
            
            An = symfun(A.*An,[x;u]); Bn = symfun(B.*Bn,[x;u]); Cn = symfun(C.*Cn,[x;u]);
            obj.LinearModel.Normalized.A = @(xss,uss) double(An(xss(1),xss(2),xss(3),xss(4),xss(5),uss(1),uss(2),uss(3),uss(4)));
            obj.LinearModel.Normalized.B = @(xss,uss) double(Bn(xss(1),xss(2),xss(3),xss(4),xss(5),uss(1),uss(2),uss(3),uss(4)));
            obj.LinearModel.Normalized.C = @(xss,uss) double(Cn(xss(1),xss(2),xss(3),xss(4),xss(5),uss(1),uss(2),uss(3),uss(4)));
            
        end
        
        function obj = Envelope_afinidade(obj)
                        
            H0_dt = -1.2454e6*obj.parameters.q0_dt.^2 + 7.4959e3*obj.parameters.q0_dt + 9.5970e2;
            H0_dt = obj.parameters.CH*H0_dt*(obj.parameters.f0/obj.parameters.f0).^2;
            H0_ut = -1.2454e6*obj.parameters.q0_ut.^2 + 7.4959e3*obj.parameters.q0_ut + 9.5970e2;
            H0_ut = obj.parameters.CH*H0_ut*(obj.parameters.f0/obj.parameters.f0).^2;
            
            % Variacao frequencia
            f = linspace(30,70,1000); % Hz
            H_ut = H0_ut*(f./obj.parameters.f0).^2;
            H_dt = H0_dt*(f./obj.parameters.f0).^2;
            % corrige lei da afinidade
            Qdt = obj.parameters.q0_dt.*f/obj.parameters.f0;
            Qut = obj.parameters.q0_ut.*f/obj.parameters.f0;
            % Variacao vazao
            flim = 35:1:65;
            qop = linspace(0,obj.parameters.q0_ut*flim(end)/obj.parameters.f0,1000); % m3/s
            Hop = zeros(length(flim),length(qop));
            for i = 1:length(flim)
                q0 = qop./obj.parameters.Cq*(obj.parameters.f0/flim(i));
                H0 = -1.2454e6*q0.^2 + 7.4959e3*q0 + 9.5970e2;
                Hop(i,:) = obj.parameters.CH*H0*(flim(i)/obj.parameters.f0).^2;
            end
            [ip(1,1),ip(1,2)] = polyxpoly(qop*3600,Hop(1,:),Qdt*3600,H_dt);
            [ip(2,1),ip(2,2)] = polyxpoly(Qdt*3600,H_dt,qop*3600,Hop(end,:));
            [ip(3,1),ip(3,2)] = polyxpoly(qop*3600,Hop(end,:),Qut*3600,H_ut);
            [ip(4,1),ip(4,2)] = polyxpoly(Qut*3600,H_ut,qop*3600,Hop(1,:));

            p_35hz = polyfit(qop*3600,Hop(1,:),3);
            H_35hz = @(qk) p_35hz*[cumprod(repmat(qk,length(p_35hz)-1,1),1,'reverse');ones(1,length(qk))];
            q_35hz = linspace(ip(1,1),ip(4,1),100);

            p_65hz = polyfit(qop*3600,Hop(end,:),3);
            H_65hz = @(qk) p_65hz*[cumprod(repmat(qk,length(p_65hz)-1,1),1,'reverse');ones(1,length(qk))];
            q_65hz = linspace(ip(2,1),ip(3,1),100);

            p_dt = polyfit(Qdt*3600,H_dt,2);
            H_dt = @(qk) p_dt*[cumprod(repmat(qk,length(p_dt)-1,1),1,'reverse');ones(1,length(qk))];
            q_dt = linspace(ip(1,1),ip(2,1),100);

            p_ut = polyfit(Qut*3600,H_ut,2);
            H_ut = @(qk) p_ut*[cumprod(repmat(qk,length(p_ut)-1,1),1,'reverse');ones(1,length(qk))];
            q_ut = linspace(ip(4,1),ip(3,1),100);
            
            obj.envelope.fig = @(aux) plot(q_35hz,H_35hz(q_35hz),':r',q_65hz,H_65hz(q_65hz),':r',q_ut,H_ut(q_ut),':r',q_dt,H_dt(q_dt),':r','LineWidth',2);
            
            fBounds = struct('H_35hz',H_35hz,'H_65hz',H_65hz,'H_dt',H_dt,'H_ut',H_ut);
            obj.envelope.Hlim = @(qk) obj.BoundHead(qk*3600,ip,fBounds);
        end
        
        function obj = OperationalEnvelop(obj)
            f = linspace(35,65,30); % Limite das frequencias
            z = linspace(15,90,30); % Limite da abertura de vavula [Down, Up]
            xss = obj.stationary.xss;
            dess = obj.stationary.uss(3:4);
            Q = zeros(length(f),length(z));
            H = Q;
            for j = 1:length(f)
                for k = 1:length(z)
                    u = [f(j);z(k);dess];
                    [~,x,~,~] = obj.SimModel(xss,u,[0,100]);
                    ss = obj.StateSolve(x(:,end),u);
                    Q(j,k) = ss.xss(3);
                    H(j,k) = ss.yss(2);
                end
            end
            %Polyfit Q vs H
            % Curva de 35 Hz
            Q35 = (Q(1,:)*3600)'; H35 = H(1,:)';
            P35 = polyfit(Q35,H35,2);
            % Curva de 65 Hz
            Q65 = (Q(end,:)*3600)'; H65 = H(end,:)';
            P65 = polyfit(Q65,H65,2);
            % Curva de Downthrust
            Qdown = (Q(:,1)*3600)'; Hdown = H(:,1)';
            Pdown = polyfit(Qdown,Hdown,2);
            % Curva de Upthrust
            Qup = (Q(:,end)*3600)'; Hup = H(:,end)';
            Pup = polyfit(Qup,Hup,2);
            
            obj.envelope.dados = struct('Q35',Q35,'Q65',Q65,'Qdown',Qdown,'Qup',Qup,'H35',H35,'H65',H65,'Hdown',Hdown,'Hup',Hup);            
            
            % Polyfit f vs H
            F35 = 35*ones(1,30); F65 = 65*ones(1,30);
             Pfdown = polyfit(f,Hdown,2); Pfup = polyfit(f,Hup,2);
            invPfdown = polyfit(Hdown,f,2); invPfup = polyfit(Hup,f,2);%polinomio inverso que calcula o fmin e fmax com base no valor de H
            
            obj.envelope.polinomio = struct('P35',P35,'P65',P65,'Pdown',Pdown,'Pup',Pup,'Pfdown',Pfdown,'Pfup',Pfup,'invPfdown',invPfdown,'invPfup',invPfup);
            
            % Envelopes
            obj.envelope.fig.qH = @(aux) plot(Q35,polyval(P35,Q35),'-.r',Q65,polyval(P65,Q65),'-.r',...
                                           Qdown,polyval(Pdown,Qdown),'-.r',Qup,polyval(Pup,Qup),'-.r');
            obj.envelope.fig.fH = @(aux) plot(F35,H35,'--r',f,Hdown,'--r',F65,H65,'--r',f,Hup,'--r');

        end
        function bounds = BoundHead(obj,x,flag)
            % Para flag = 'qH', x deve ser vazão em m3/s
            % Para flag = 'fH', x deve ser frequencia em Hz
            pol = obj.envelope.polinomio;
            switch flag
                case 'qH'
                    qk = x*3600;%x=q
                    dados = obj.envelope.dados;
                    if qk <= min(dados.Q35)
                        Hlim = [polyval(pol.P35,min(dados.Q35)),polyval(pol.P35,min(dados.Q35))];
                    elseif qk <= max(dados.Qdown)
                        Hlim = [polyval(pol.P35,qk),polyval(pol.Pdown,qk)];
                    elseif qk <= max(dados.Q35)
                        Hlim = [polyval(pol.P35,qk),polyval(pol.P65,qk)];
                    elseif qk <= max(dados.Qup)
                        Hlim = [polyval(pol.Pup,qk),polyval(pol.P65,qk)];
                    else
                        Hlim = [polyval(pol.P65,max(dados.Qup)),polyval(pol.P65,max(dados.Qup))];
                    end
                case 'fH'
                    Hlim = [polyval(pol.Pfdown,x),polyval(pol.Pfup,x)];%x=f
                case 'Hf' %calcula fmin e fmax para o valor de H atual
                    Hlim = [polyval(pol.invPfdown,x),polyval(pol.invPfup,x)];%x=H
                    
            end
            
            bounds.min = min(Hlim);
            bounds.max = max(Hlim);
        end
        
        function [Kf,ymk,xmk,Pk] = EKF(obj,xpk,ypk,upk_1,W,V,Pk,ts,pv)
            % Linearizacao a cada k Atualização da matriz A
            Ak = obj.LinearModel.A(xpk,upk_1);

            % EKF - Predicao
            [ymk,xmk] = obj.SimModel(xpk,upk_1,[0,ts]);
            % Calculo da matriz de covariancia Mk
            Phi = eye(length(xmk)) + Ak*ts + (Ak^2)*(ts^2)/2 + (Ak^3)*(ts^3)/factorial(3);
            Pk = Phi*Pk*Phi' + W;

            % Linearizacao a cada k Atualização da matriz C
            Ck = obj.LinearModel.C(xmk,upk_1);

            % EKF - correcao dos estados estimados
            % calculdo ganho
            Kf = Pk*Ck(pv,:)'/(Ck(pv,:)*Pk*Ck(pv,:)' + V) ;
            % atualizacao da matriz de variancia dos estados estimados
            Pk = (eye(length(xmk)) - Kf*Ck(pv,:))*Pk;
            % correcao do estado
            xmk = xmk + Kf*(ypk(pv) - ymk(pv));
            obj = obj.EvaluateModel(xmk,upk_1);
            ymk = obj.outputs.y;
        end
    end
end