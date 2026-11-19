classdef IEEE57_MOOPF_LEGACY < PROBLEM
%IEEE57_MOOPF IEEE57 V2四目标最优潮流问题�?%
% CaseID�?成本+损耗；2成本+排放�?成本+损�?排放�?%         4成本+损�?排放+LVD�?% UseFeasibleSeed：正式比较固定为0�?仅用于显式诊断测试�?
    properties(SetAccess=private)
        Config
        CaseID
        ObjectiveMask
        UseFeasibleSeed
        SharedInitialDec
        SharedInitializationSeed
        SharedInitializationCallCount = 0
        SharedInitializationRequestSizes = []
        SharedPrimaryPoolUseCount = 0
    end

    methods
        function Setting(obj)
            obj.Config = ieee57_cdmtmo_config();
            [obj.CaseID,obj.UseFeasibleSeed] = obj.ParameterSet(4,0);
            if ~isscalar(obj.CaseID) || ~isnumeric(obj.CaseID) || ~isfinite(obj.CaseID) || obj.CaseID~=round(obj.CaseID) || obj.CaseID<1 || obj.CaseID>4
                error('IEEE57:InvalidCaseID','CaseID must be 1, 2, 3, or 4.');
            end
            if ~isscalar(obj.UseFeasibleSeed) || ~isnumeric(obj.UseFeasibleSeed) || ~isfinite(obj.UseFeasibleSeed) || ~ismember(obj.UseFeasibleSeed,[0,1])
                error('IEEE57:InvalidUseFeasibleSeed','UseFeasibleSeed must be 0 or 1.');
            end
            obj.CaseID=round(obj.CaseID);
            obj.UseFeasibleSeed=double(obj.UseFeasibleSeed);
            obj.ObjectiveMask=obj.Config.case.objective_mask(obj.CaseID,:);
            obj.M=sum(obj.ObjectiveMask);
            obj.D=obj.Config.D;
            obj.lower=obj.Config.variable.lower';
            obj.upper=obj.Config.variable.upper';
            obj.encoding=ones(1,obj.D);
            obj.SharedInitialDec=[];
            obj.SharedInitializationSeed=[];
            obj.SharedInitializationCallCount=0;
            obj.SharedInitializationRequestSizes=[];
            obj.SharedPrimaryPoolUseCount=0;
        end

        function SetSharedInitialDec(obj,X)
            if isempty(obj.Config), obj.Config=ieee57_cdmtmo_config(); end
            obj.validateSharedInitialMatrix(X);
            Xdecoded=ieee57_moopf_decode(X,obj.Config);
            if ~isequal(size(Xdecoded),[obj.N,obj.D]), error('IEEE57:SharedDecodeSizeChanged','Shared initial matrix decoded to an unexpected size.'); end
            obj.SharedInitialDec=Xdecoded;
            obj.SharedInitializationSeed=[];
            obj.SharedInitializationCallCount=0;
            obj.SharedInitializationRequestSizes=[];
            obj.SharedPrimaryPoolUseCount=0;
        end

        function SetSharedInitialStream(obj,Xprimary,randomSeed)
            if isempty(obj.Config), obj.Config=ieee57_cdmtmo_config(); end
            obj.validateSharedInitialMatrix(Xprimary);
            if ~isscalar(randomSeed) || ~isnumeric(randomSeed) || ~isfinite(randomSeed) || randomSeed~=round(randomSeed) || randomSeed<0 || randomSeed>2^32-1
                error('IEEE57:InvalidSharedSeed','Shared initialization seed must be an integer in [0,2^32-1].');
            end
            Xdecoded=ieee57_moopf_decode(Xprimary,obj.Config);
            if ~isequal(size(Xdecoded),[obj.N,obj.D]), error('IEEE57:SharedMainDecodeSizeChanged','Shared main population decoded to an unexpected size.'); end
            obj.SharedInitialDec=Xdecoded;
            obj.SharedInitializationSeed=double(randomSeed);
            obj.SharedInitializationCallCount=0;
            obj.SharedInitializationRequestSizes=[];
            obj.SharedPrimaryPoolUseCount=0;
        end

        function ClearSharedInitialDec(obj)
            obj.SharedInitialDec=[];
            obj.SharedInitializationSeed=[];
            obj.SharedInitializationCallCount=0;
            obj.SharedInitializationRequestSizes=[];
            obj.SharedPrimaryPoolUseCount=0;
        end

        function Population=Initialization(obj,N)
            if nargin<2 || isempty(N), N=obj.N; end
            if ~isscalar(N) || ~isnumeric(N) || ~isfinite(N) || N~=round(N) || N<1, error('IEEE57:InvalidN','N must be a positive integer.'); end
            N=round(N); lowerRow=reshape(obj.lower,1,obj.D); upperRow=reshape(obj.upper,1,obj.D);
            if ~isempty(obj.SharedInitialDec)
                obj.SharedInitializationCallCount=obj.SharedInitializationCallCount+1;
                callID=obj.SharedInitializationCallCount;
                obj.SharedInitializationRequestSizes(end+1,1)=N;
                if isempty(obj.SharedInitializationSeed)
                    if ~isequal(size(obj.SharedInitialDec),[N,obj.D]), error('IEEE57:SharedInitialSizeMismatch','Shared initial matrix size does not match requested N and D.'); end
                    PopDec=obj.SharedInitialDec;
                    obj.SharedPrimaryPoolUseCount=obj.SharedPrimaryPoolUseCount+1;
                elseif callID==1 && N==size(obj.SharedInitialDec,1)
                    PopDec=obj.SharedInitialDec;
                    obj.SharedPrimaryPoolUseCount=obj.SharedPrimaryPoolUseCount+1;
                else
                    auxiliarySeed=obj.calculateAuxiliarySeed(callID,N);
                    stream=RandStream('mt19937ar','Seed',auxiliarySeed);
                    PopDec=repmat(lowerRow,N,1)+rand(stream,N,obj.D).*repmat(upperRow-lowerRow,N,1);
                end
            else
                PopDec=repmat(lowerRow,N,1)+rand(N,obj.D).*repmat(upperRow-lowerRow,N,1);
                if obj.UseFeasibleSeed==1
                    if ~obj.Config.reference.feasible_seed_available
                        error('IEEE57:NoFeasibleSeed','No verified IEEE57 feasible diagnostic seed is available.');
                    end
                    PopDec(1,:)=obj.Config.reference.feasible_seed_dec;
                end
            end
            Population=obj.Evaluation(PopDec);
        end

        function PopDec=CalDec(obj,PopDec)
            if isempty(obj.Config), obj.Config=ieee57_cdmtmo_config(); end
            PopDec=ieee57_moopf_decode(PopDec,obj.Config);
        end

        function Population=Evaluation(obj,varargin)
            if isempty(varargin)
                error('IEEE57:EvaluationInputMissing','Evaluation requires at least one decision matrix.');
            end

            requestedCount = size(varargin{1},1);
            if requestedCount <= 0
                error('IEEE57:EmptyEvaluationRequest','Evaluation request contains no candidate solutions.');
            end

            if isempty(obj.maxFE) || ~isfinite(obj.maxFE)
                remainingEvaluations = requestedCount;
            else
                remainingEvaluations = max(0,obj.maxFE-obj.FE);
            end
            acceptedCount = min(requestedCount,remainingEvaluations);

            if acceptedCount < requestedCount
                if acceptedCount <= 0
                    error('IEEE57:ExactFEBudgetExhausted','Current FE has reached MaxFE, but the algorithm requested additional evaluations.');
                end
                keepIndex = 1:acceptedCount;
                varargin{1} = varargin{1}(keepIndex,:);
                for argumentIndex = 2:numel(varargin)
                    value = varargin{argumentIndex};
                    try
                        if istable(value)
                            if height(value) == requestedCount
                                value = value(keepIndex,:);
                            end
                        elseif ~isscalar(value)
                            if size(value,1) == requestedCount
                                subscripts = repmat({':'},1,ndims(value));
                                subscripts{1} = keepIndex;
                                value = value(subscripts{:});
                            elseif isvector(value) && numel(value) == requestedCount
                                value = value(keepIndex);
                            end
                        end
                    catch argumentError
                        error('IEEE57:ExactFEArgumentTruncationFailed', ...
                            'Failed to truncate Evaluation argument %d: %s', ...
                            argumentIndex,argumentError.message);
                    end
                    varargin{argumentIndex} = value;
                end
            end

            PopDec = obj.CalDec(varargin{1});
            [allObj,PopCon,~,PopDecEvaluated] = ieee57_moopf_evaluate_legacy(PopDec,obj.Config);
            PopObj = allObj(:,obj.ObjectiveMask);
            Population = SOLUTION(PopDecEvaluated,PopObj,PopCon,varargin{2:end});

            evaluatedCount = length(Population);
            if evaluatedCount ~= acceptedCount
                error('IEEE57:ExactFEPopulationCountMismatch', ...
                    'SOLUTION count %d differs from accepted evaluation count %d.', ...
                    evaluatedCount,acceptedCount);
            end
            if isfinite(obj.maxFE) && obj.FE+evaluatedCount > obj.maxFE
                error('IEEE57:ExactFEInternalCheckFailed','FE update would exceed MaxFE.');
            end
            obj.FE = obj.FE + evaluatedCount;
        end

        function PopObj=CalObj(obj,PopDec)
            PopDec=obj.CalDec(PopDec);
            [allObj,~]=ieee57_moopf_evaluate_legacy(PopDec,obj.Config);
            PopObj=allObj(:,obj.ObjectiveMask);
        end

        function PopCon=CalCon(obj,PopDec)
            PopDec=obj.CalDec(PopDec);
            [~,PopCon]=ieee57_moopf_evaluate_legacy(PopDec,obj.Config);
        end
    end

    methods(Access=private)
        function validateSharedInitialMatrix(obj,X)
            if ~isnumeric(X) || ~isreal(X) || isempty(X) || any(~isfinite(X(:)))
                error('IEEE57:InvalidSharedInitialDec','Shared initial decision matrix must be a nonempty finite real matrix.');
            end
            if ~isequal(size(X),[obj.N,obj.D])
                error('IEEE57:SharedInitialSizeMismatch','Shared initial decision matrix must be %d-by-%d, but is %d-by-%d.',obj.N,obj.D,size(X,1),size(X,2));
            end
        end

        function auxiliarySeed=calculateAuxiliarySeed(obj,callID,N)
            if isempty(obj.SharedInitializationSeed), error('IEEE57:MissingSharedInitializationSeed','Shared initialization seed is empty.'); end
            modulusValue=uint64(2^32-1);
            seed64=mod(uint64(obj.SharedInitializationSeed)+uint64(104729)*uint64(callID)+uint64(1009)*uint64(N),modulusValue);
            if seed64==0, seed64=uint64(1); end
            auxiliarySeed=double(seed64);
        end
    end
end
