classdef IEEE118_MOOPF_Scalability < PROBLEM
% <multi> <real> <constrained>
    properties(SetAccess=private)
        Config
        PFCalls=0
        PFFailures=0
        InitialCallCount=0
        InitializationSeed=[]
        InitialPopulation1=[]
        InitialPopulation2=[]
    end
    methods
        function Setting(obj)
            obj.Config=ieee118_fsa_config();obj.D=obj.Config.D;obj.M=2;
            obj.lower=obj.Config.variable.lower;obj.upper=obj.Config.variable.upper;
            obj.encoding=ones(1,obj.D);
        end
        function SetInitializationSeed(obj,seed)
            obj.InitializationSeed=seed;obj.InitialCallCount=0;
        end
        function P=Initialization(obj,N)
            if nargin<2,N=obj.N;end
            obj.InitialCallCount=obj.InitialCallCount+1;
            if isempty(obj.InitializationSeed)
                U=rand(N,obj.D);
            else
                seed=mod(obj.InitializationSeed+104729*(obj.InitialCallCount-1),2^32-1);
                stream=RandStream('mt19937ar','Seed',seed);U=rand(stream,N,obj.D);
            end
            X=obj.CalDec(obj.lower+U.*(obj.upper-obj.lower));
            if obj.InitialCallCount==1,obj.InitialPopulation1=X;end
            if obj.InitialCallCount==2,obj.InitialPopulation2=X;end
            P=obj.Evaluation(X);
        end
        function X=CalDec(obj,X)
            X=ieee118_fsa_decode(X,obj.Config);
        end
        function P=Evaluation(obj,varargin)
            X=varargin{1};requested=size(X,1);remaining=max(0,obj.maxFE-obj.FE);
            accepted=min(requested,remaining);
            if accepted==0
                error('PlatEMO:Termination','The common exact FE budget has been exhausted.');
            end
            varargin{1}=X(1:accepted,:);
            if accepted<requested
                for k=2:numel(varargin)
                    val=varargin{k};
                    if size(val,1)==requested,val=val(1:accepted,:);end
                    varargin{k}=val;
                end
            end
            [F,G,detail,X]=ieee118_fsa_evaluate(varargin{1},obj.Config);
            P=SOLUTION(X,F,G,varargin{2:end});
            obj.FE=obj.FE+accepted;obj.PFCalls=obj.PFCalls+accepted;
            obj.PFFailures=obj.PFFailures+sum(~detail.success);
            assert(obj.FE<=obj.maxFE && obj.PFCalls==obj.FE);
        end
        function F=CalObj(obj,X)
            [F,~]=ieee118_fsa_evaluate(X,obj.Config);
        end
        function G=CalCon(obj,X)
            [~,G]=ieee118_fsa_evaluate(X,obj.Config);
        end
        function R=GetOptimum(obj,~),R=ones(1,obj.M);end
        function R=GetPF(~),R=[];end
    end
end
