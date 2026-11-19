
classdef IEEE30_MOOPF_LEGACY < PROBLEM
%IEEE30_MOOPF IEEE 30节点多目标最优潮流问�?
%
% 基础模型�?
%   - MATPOWER case_ieee30�?
%   - case30工程线路容量�?
%   - 24个控制变量；
%   - 4个可选目标；
%   - 113个状态约束�?
%
% 参数�?
%   CaseID
%       1：FuelCost + ActivePowerLoss
%       2：FuelCost + Emission
%       3：FuelCost + ActivePowerLoss + Emission
%       4：FuelCost + ActivePowerLoss + Emission + LVD
%
%   UseFeasibleSeed
%       0：完全随机初始化
%       1：第一个个体固定为已验证可行的Case 4参考解
%
% 默认参数�?
%   CaseID = 4
%   UseFeasibleSeed = 1
%
% 正式算法比较建议�?
%   CaseID = 4
%   UseFeasibleSeed = 0
%
% 共享初始化接口：
%
%   SetSharedInitialDec(X)
%       传统单种群共享初始化�?
%       每次Initialization调用都使用同一个X�?
%       主要用于NSGA-II、SPEA2等单次初始化算法的审计�?
%
%   SetSharedInitialStream(X,seed)
%       多种群约束算法共享初始化流�?
%       第一次标准规模初始化使用冻结主种群X�?
%       后续初始化调用使用由seed和调用序号确定的辅助种群�?
%
% 约束规则�?
%   PopCon <= 0 表示满足约束�?
%   PopCon > 0 表示约束违反�?

    properties(SetAccess = private)

        %% 基础配置
        Config
        CaseID
        ObjectiveMask
        UseFeasibleSeed

        %% 共享主初始种�?
        SharedInitialDec

        %% 共享初始化流审计信息
        SharedInitializationSeed
        SharedInitializationCallCount = 0
        SharedInitializationRequestSizes = []
        SharedPrimaryPoolUseCount = 0
    end

    methods

        %% ============================================================
        % 问题参数设置
        % =============================================================
        function Setting(obj)

            obj.Config = ieee30_cdmtmo_config();

            % 默认�?
            % CaseID = 4
            % UseFeasibleSeed = 1
            [obj.CaseID,obj.UseFeasibleSeed] = ...
                obj.ParameterSet(4,1);

            %% CaseID检�?

            if ~isscalar(obj.CaseID) || ...
                    ~isnumeric(obj.CaseID) || ...
                    ~isfinite(obj.CaseID) || ...
                    obj.CaseID ~= round(obj.CaseID) || ...
                    obj.CaseID < 1 || ...
                    obj.CaseID > 4

                error('IEEE30:InvalidCaseID','CaseID must be 1, 2, 3, or 4.');
            end

            obj.CaseID = round(obj.CaseID);

            %% 可行种子开关检�?

            if ~isscalar(obj.UseFeasibleSeed) || ...
                    ~isnumeric(obj.UseFeasibleSeed) || ...
                    ~isfinite(obj.UseFeasibleSeed) || ...
                    ~ismember(obj.UseFeasibleSeed,[0,1])

                error('IEEE30:InvalidUseFeasibleSeed','UseFeasibleSeed must be 0 or 1.');
            end

            obj.UseFeasibleSeed = ...
                double(obj.UseFeasibleSeed);

            %% 当前目标组合

            obj.ObjectiveMask = ...
                obj.Config.case.objective_mask( ...
                    obj.CaseID,:);

            obj.M = sum(obj.ObjectiveMask);
            obj.D = obj.Config.D;

            %% PlatEMO问题属�?

            % PlatEMO要求上下界为行向�?
            obj.lower = ...
                obj.Config.variable.lower';

            obj.upper = ...
                obj.Config.variable.upper';

            % 算法层统一使用实数编码�?
            % TAP与QC由CalDec进行离散化�?
            obj.encoding = ones(1,obj.D);

            %% 初始化共享状�?

            obj.SharedInitialDec = [];
            obj.SharedInitializationSeed = [];

            obj.SharedInitializationCallCount = 0;
            obj.SharedInitializationRequestSizes = [];
            obj.SharedPrimaryPoolUseCount = 0;
        end

        %% ============================================================
        % 设置单一共享初始决策矩阵
        % =============================================================
        function SetSharedInitialDec(obj,X)
        %SETSHAREDINITIALDEC 设置传统单种群共享初始种�?
        %
        % 输入X必须为obj.N × obj.D�?
        %
        % 设置本模式后�?
        %   - Initialization不调用rand�?
        %   - 每次Initialization调用均返回同一个共享种群；
        %   - 适用于单种群算法和共享初始化公平性审计�?
        %
        % 对于可能多次调用Initialization的多种群算法�?
        % 应改用SetSharedInitialStream�?

            if isempty(obj.Config)
                obj.Config = ieee30_cdmtmo_config();
            end

            obj.validateSharedInitialMatrix(X);

            Xdecoded = ...
                ieee30_moopf_decode(X,obj.Config);

            if ~isequal(size(Xdecoded),[obj.N,obj.D])
                error('IEEE30:SharedDecodeSizeChanged','Shared initial matrix decoded to an unexpected size.');
            end

            obj.SharedInitialDec = Xdecoded;

            % Seed为空表示采用传统固定共享种群模式
            obj.SharedInitializationSeed = [];

            obj.SharedInitializationCallCount = 0;
            obj.SharedInitializationRequestSizes = [];
            obj.SharedPrimaryPoolUseCount = 0;
        end

        %% ============================================================
        % 设置多种群共享初始化�?
        % =============================================================
        function SetSharedInitialStream(obj,Xprimary,randomSeed)
        %SETSHAREDINITIALSTREAM 设置跨算法共享初始化�?
        %
        % 第一次标准规模初始化�?
        %   使用冻结的共享主初始种群Xprimary�?
        %
        % 第二次及后续初始化：
        %   使用由randomSeed、初始化调用序号和请求规模共同决定的
        %   确定性辅助种群�?
        %
        % 这样可以�?
        %   1. 保证不同算法的第一次主种群完全一致；
        %   2. 避免多种群算法的多个内部种群完全重复�?
        %   3. 保证辅助种群可重复；
        %   4. 不改变算法自身的全局随机数状态�?

            if isempty(obj.Config)
                obj.Config = ieee30_cdmtmo_config();
            end

            obj.validateSharedInitialMatrix(Xprimary);

            if ~isscalar(randomSeed) || ...
                    ~isnumeric(randomSeed) || ...
                    ~isfinite(randomSeed) || ...
                    randomSeed ~= round(randomSeed) || ...
                    randomSeed < 0 || ...
                    randomSeed > 2^32-1

                error('IEEE30:InvalidSharedSeed','Shared initialization seed must be an integer in [0,2^32-1].');
            end

            Xdecoded = ...
                ieee30_moopf_decode( ...
                    Xprimary,obj.Config);

            if ~isequal(size(Xdecoded),[obj.N,obj.D])
                error('IEEE30:SharedMainDecodeSizeChanged','Shared main initial population decoded to an unexpected size.');
            end

            obj.SharedInitialDec = Xdecoded;

            obj.SharedInitializationSeed = ...
                double(randomSeed);

            obj.SharedInitializationCallCount = 0;
            obj.SharedInitializationRequestSizes = [];
            obj.SharedPrimaryPoolUseCount = 0;
        end

        %% ============================================================
        % 清除共享初始化设�?
        % =============================================================
        function ClearSharedInitialDec(obj)
        %CLEARSHAREDINITIALDEC 恢复问题类自身的普通随机初始化

            obj.SharedInitialDec = [];
            obj.SharedInitializationSeed = [];

            obj.SharedInitializationCallCount = 0;
            obj.SharedInitializationRequestSizes = [];
            obj.SharedPrimaryPoolUseCount = 0;
        end

        %% ============================================================
        % 统一初始�?
        % =============================================================
        function Population = Initialization(obj,N)
        %INITIALIZATION 生成IEEE 30节点MOOPF初始种群
        %
        % 模式一：未设置共享种群
        %   在上下界内均匀随机初始化；
        %   UseFeasibleSeed=1时，将第一个个体替换为Case 4�?
        %
        % 模式二：SetSharedInitialDec
        %   每次调用均使用相同共享种群�?
        %
        % 模式三：SetSharedInitialStream
        %   第一次标准规模调用使用冻结主种群�?
        %   后续调用使用确定性辅助种群�?

            %% 1. 初始化规�?

            if nargin < 2 || isempty(N)
                N = obj.N;
            end

            if isempty(N) || ...
                    ~isscalar(N) || ...
                    ~isnumeric(N) || ...
                    ~isfinite(N) || ...
                    N ~= round(N) || ...
                    N < 1

                error('IEEE30:InvalidN','Initial population size N must be a positive integer.');
            end

            N = round(N);

            lowerRow = reshape(obj.lower,1,obj.D);
            upperRow = reshape(obj.upper,1,obj.D);

            %% 2. 已设置共享种�?

            if ~isempty(obj.SharedInitialDec)

                obj.SharedInitializationCallCount = ...
                    obj.SharedInitializationCallCount + 1;

                callID = ...
                    obj.SharedInitializationCallCount;

                obj.SharedInitializationRequestSizes(end+1,1) = N;

                %% 2.1 传统固定共享种群模式

                if isempty(obj.SharedInitializationSeed)

                    if ~isequal( ...
                            size(obj.SharedInitialDec), ...
                            [N,obj.D])

                        error('IEEE30:SharedInitialSizeMismatch', ...
                            'Shared initial matrix is %d-by-%d, but requested initialization is %d-by-%d.', ...
                            size(obj.SharedInitialDec,1), ...
                            size(obj.SharedInitialDec,2), ...
                            N,obj.D);
                    end

                    PopDec = obj.SharedInitialDec;

                    obj.SharedPrimaryPoolUseCount = ...
                        obj.SharedPrimaryPoolUseCount + 1;

                else
                    %% 2.2 多种群共享初始化流模�?

                    % 第一次标准规模初始化使用冻结主种�?
                    if callID == 1 && ...
                            N == size(obj.SharedInitialDec,1)

                        PopDec = obj.SharedInitialDec;

                        obj.SharedPrimaryPoolUseCount = ...
                            obj.SharedPrimaryPoolUseCount + 1;

                    else
                        % 后续初始化调用使用确定性辅助种群�?
                        %
                        % 采用局部RandStream，不改变MATLAB全局rng状态，
                        % 因而不会扰动算法后续的交叉、变异等随机操作�?

                        auxiliarySeed = ...
                            obj.calculateAuxiliarySeed( ...
                                callID,N);

                        auxiliaryStream = RandStream( ...
                            'mt19937ar', ...
                            'Seed',auxiliarySeed);

                        PopDec = ...
                            repmat(lowerRow,N,1) + ...
                            rand(auxiliaryStream,N,obj.D) .* ...
                            repmat(upperRow-lowerRow,N,1);
                    end
                end

            else
                %% 3. 普通随机初始化

                PopDec = ...
                    repmat(lowerRow,N,1) + ...
                    rand(N,obj.D) .* ...
                    repmat(upperRow-lowerRow,N,1);

                %% 4. 可选可行种�?

                if obj.UseFeasibleSeed == 1

                    feasibleSeed = ...
                        obj.Config.reference.table3_x(4,:);

                    if ~isequal( ...
                            size(feasibleSeed), ...
                            [1,obj.D])

                        error('IEEE30:FeasibleSeedSizeMismatch', ...
                            'Reference feasible seed size is not 1-by-%d.', ...
                            obj.D);
                    end

                    PopDec(1,:) = feasibleSeed;
                end
            end

            %% 5. 统一解码与评�?

            Population = ...
                obj.Evaluation(PopDec);
        end

        %% ============================================================
        % 决策变量统一解码
        % =============================================================
        function PopDec = CalDec(obj,PopDec)

            if isempty(obj.Config)
                obj.Config = ...
                    ieee30_cdmtmo_config();
            end

            PopDec = ...
                ieee30_moopf_decode( ...
                    PopDec, ...
                    obj.Config);
        end

        %% ============================================================
        % 统一评价入口
        % =============================================================
        function Population = Evaluation(obj,varargin)
        %EVALUATION 一次潮流同时计算目标和约束
        %
        % 覆盖PROBLEM.Evaluation，避免CalObj和CalCon分别调用
        % MATPOWER，防止同一候选解重复执行潮流计算�?

            if isempty(varargin)
                error('IEEE30:EvaluationInputMissing','Evaluation requires at least one decision matrix.');
            end
            requestedCount = size(varargin{1},1);
            if isempty(obj.maxFE) || ~isfinite(obj.maxFE)
                remainingEvaluations = requestedCount;
            else
                remainingEvaluations = max(0,obj.maxFE-obj.FE);
            end
            acceptedCount = min(requestedCount,remainingEvaluations);
            if acceptedCount < requestedCount
                if acceptedCount <= 0
                    error('IEEE30:ExactFEBudgetExhausted', ...
                        'Current FE has reached MaxFE, but the algorithm requested additional evaluations.');
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
                        error('IEEE30:ExactFEArgumentTruncationFailed', ...
                            'Failed to truncate Evaluation argument %d: %s', ...
                            argumentIndex,argumentError.message);
                    end
                    varargin{argumentIndex} = value;
                end
            end


            %% 1. 统一解码

            PopDec = ...
                obj.CalDec(varargin{1});

            %% 2. 一次潮流计算全部目标和约束

            [allObjectives,PopCon,~,PopDecEvaluated] = ...
                ieee30_moopf_evaluate_legacy( ...
                    PopDec, ...
                    obj.Config);

            %% 3. 选择当前Case对应的目�?

            PopObj = ...
                allObjectives(:,obj.ObjectiveMask);

            PopDec = PopDecEvaluated;

            %% 4. 构造PlatEMO解对�?

            Population = SOLUTION( ...
                PopDec, ...
                PopObj, ...
                PopCon, ...
                varargin{2:end});

            %% 5. 更新真实评价次数

            evaluatedCount = length(Population);
            if evaluatedCount ~= acceptedCount
                error('IEEE30:ExactFEPopulationCountMismatch', ...
                    'SOLUTION count %d differs from accepted evaluation count %d.', ...
                    evaluatedCount,acceptedCount);
            end
            if isfinite(obj.maxFE) && obj.FE+evaluatedCount > obj.maxFE
                error('IEEE30:ExactFEInternalCheckFailed', ...
                    'FE update would exceed MaxFE.');
            end

            obj.FE = ...
                obj.FE + evaluatedCount;
        end

        %% ============================================================
        % 目标计算兼容接口
        % =============================================================
        function PopObj = CalObj(obj,PopDec)

            PopDec = ...
                obj.CalDec(PopDec);

            [allObjectives,~] = ...
                ieee30_moopf_evaluate_legacy( ...
                    PopDec, ...
                    obj.Config);

            PopObj = ...
                allObjectives(:,obj.ObjectiveMask);
        end

        %% ============================================================
        % 约束计算兼容接口
        % =============================================================
        function PopCon = CalCon(obj,PopDec)

            PopDec = ...
                obj.CalDec(PopDec);

            [~,PopCon] = ...
                ieee30_moopf_evaluate_legacy( ...
                    PopDec, ...
                    obj.Config);
        end
    end

    methods(Access = private)

        %% ============================================================
        % 共享初始矩阵输入检�?
        % =============================================================
        function validateSharedInitialMatrix(obj,X)

            if ~isnumeric(X) || ...
                    ~isreal(X) || ...
                    isempty(X)

                error('IEEE30:InvalidSharedInitialDec','Shared initial decision matrix must be a nonempty real matrix.');
            end

            if any(~isfinite(X(:)))
                error('IEEE30:SharedInitialDecNonfinite','Shared initial decision matrix contains NaN or Inf.');
            end

            if ~isequal(size(X),[obj.N,obj.D])

                error('IEEE30:SharedInitialSizeMismatch', ...
                    'Shared initial decision matrix must be %d-by-%d, but is %d-by-%d.', ...
                    obj.N,obj.D, ...
                    size(X,1),size(X,2));
            end
        end

        %% ============================================================
        % 计算确定性辅助种群随机种�?
        % =============================================================
        function auxiliarySeed = ...
                calculateAuxiliarySeed(obj,callID,N)
        %CALCULATEAUXILIARYSEED
        %
        % 根据�?
        %   - 主随机种子；
        %   - 初始化调用序号；
        %   - 当前请求种群规模�?
        %
        % 生成可重复的辅助种群随机种子�?

            if isempty(obj.SharedInitializationSeed)
                error('IEEE30:MissingSharedInitializationSeed','Shared initialization seed is empty.');
            end

            if ~isscalar(callID) || ...
                    callID < 1 || ...
                    callID ~= round(callID)

                error('IEEE30:InvalidInitializationCallID','Initialization call ID is invalid.');
            end

            modulusValue = uint64(2^32-1);

            seed64 = mod( ...
                uint64(obj.SharedInitializationSeed) + ...
                uint64(104729)*uint64(callID) + ...
                uint64(1009)*uint64(N), ...
                modulusValue);

            if seed64 == 0
                seed64 = uint64(1);
            end

            auxiliarySeed = double(seed64);
        end
    end
end

