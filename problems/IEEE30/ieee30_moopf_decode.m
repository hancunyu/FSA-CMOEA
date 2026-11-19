
function [Xdecoded,info] = ieee30_moopf_decode(X,cfg)
%IEEE30_MOOPF_DECODE IEEE 30节点MOOPF决策变量统一解码
%
% 用法：
%   Xdecoded = ieee30_moopf_decode(X);
%   [Xdecoded,info] = ieee30_moopf_decode(X,cfg);
%
% 输入：
%   X   - N×24决策矩阵、1×24行向量或24×1列向量
%   cfg - ieee30_cdmtmo_config返回的配置结构体
%
% 输出：
%   Xdecoded - 完成边界修复和离散化后的决策变量
%   info     - 解码诊断信息
%
% 解码规则：
%   1. 所有变量先截断到上下界；
%   2. PG、VG保持连续值；
%   3. TAP按0.02步长投影；
%   4. QC四舍五入为整数；
%   5. 最后再次执行边界保护。
%
% 本函数不执行潮流计算，不修改MATPOWER原始算例。

%% 1. 配置与输入检查

if nargin < 2 || isempty(cfg)
    cfg = ieee30_cdmtmo_config();
end

if ~isstruct(cfg) || ...
        ~isfield(cfg,'D') || ...
        ~isfield(cfg,'variable') || ...
        ~isfield(cfg,'index')
    error('cfg不是有效的IEEE 30节点MOOPF配置结构体。');
end

if ~isnumeric(X) || ~isreal(X) || isempty(X)
    error('X必须是非空实数数组。');
end

X = double(X);

if any(~isfinite(X(:)))
    error('X中存在NaN或Inf，不能执行解码。');
end

%% 2. 统一输入形状

inputWasColumnVector = false;

if isvector(X)
    if numel(X) ~= cfg.D
        error('决策向量长度必须为%d，当前为%d。', ...
            cfg.D,numel(X));
    end

    inputWasColumnVector = ...
        size(X,1) == cfg.D && size(X,2) == 1;

    Xmatrix = reshape(X,1,cfg.D);
else
    if size(X,2) ~= cfg.D
        error('决策矩阵列数必须为%d，当前为%d。', ...
            cfg.D,size(X,2));
    end

    Xmatrix = X;
end

numberOfSolutions = size(Xmatrix,1);

lower = reshape(cfg.variable.lower,1,cfg.D);
upper = reshape(cfg.variable.upper,1,cfg.D);

if any(lower > upper)
    error('配置中存在变量下界大于上界的情况。');
end

%% 3. 保存原始输入

Xoriginal = Xmatrix;

%% 4. 第一次边界截断

Xbounded = min(max(Xmatrix,lower),upper);

initialClipMask = ...
    abs(Xbounded-Xoriginal) > 1e-12;

%% 5. 变压器分接比离散化

XdecodedMatrix = Xbounded;

tapIndex = cfg.index.TAP;

for localIndex = 1:numel(tapIndex)

    variableIndex = tapIndex(localIndex);
    stepValue = cfg.variable.step(variableIndex);

    if ~isfinite(stepValue) || stepValue <= 0
        error('变量%s的离散步长无效。', ...
            cfg.variable.names(variableIndex));
    end

    XdecodedMatrix(:,variableIndex) = ...
        lower(variableIndex) + ...
        round( ...
            (XdecodedMatrix(:,variableIndex) - ...
             lower(variableIndex)) ./ stepValue) .* ...
        stepValue;

    % 消除浮点计算产生的微小尾数
    XdecodedMatrix(:,variableIndex) = ...
        round(XdecodedMatrix(:,variableIndex),12);
end

tapQuantizationMask = false(size(XdecodedMatrix));

tapQuantizationMask(:,tapIndex) = ...
    abs( ...
        XdecodedMatrix(:,tapIndex) - ...
        Xbounded(:,tapIndex)) > 1e-12;

%% 6. 并联无功补偿整数化

qcIndex = cfg.index.QC;

qcBeforeInteger = XdecodedMatrix(:,qcIndex);

XdecodedMatrix(:,qcIndex) = ...
    round(XdecodedMatrix(:,qcIndex));

qcIntegerMask = false(size(XdecodedMatrix));

qcIntegerMask(:,qcIndex) = ...
    abs( ...
        XdecodedMatrix(:,qcIndex) - ...
        qcBeforeInteger) > 1e-12;

%% 7. 最终边界保护

XdecodedMatrix = ...
    min(max(XdecodedMatrix,lower),upper);

%% 8. 最终合法性检查

if any(~isfinite(XdecodedMatrix(:)))
    error('解码结果中出现NaN或Inf。');
end

if any(any(XdecodedMatrix < lower-1e-12)) || ...
        any(any(XdecodedMatrix > upper+1e-12))
    error('解码结果仍存在边界越界。');
end

% 检查TAP是否落在离散网格上
for localIndex = 1:numel(tapIndex)

    variableIndex = tapIndex(localIndex);
    stepValue = cfg.variable.step(variableIndex);

    gridPosition = ...
        (XdecodedMatrix(:,variableIndex) - ...
         lower(variableIndex)) ./ stepValue;

    if any(abs(gridPosition-round(gridPosition)) > 1e-9)
        error('变量%s没有正确落在离散网格上。', ...
            cfg.variable.names(variableIndex));
    end
end

% 检查QC是否为整数
if any(any( ...
        abs( ...
            XdecodedMatrix(:,qcIndex) - ...
            round(XdecodedMatrix(:,qcIndex))) > 1e-12))
    error('QC变量整数化失败。');
end

%% 9. 诊断信息

finalChangeMask = ...
    abs(XdecodedMatrix-Xoriginal) > 1e-12;

info.number_of_solutions = numberOfSolutions;
info.number_of_variables = cfg.D;

info.clipped_element_count = ...
    sum(initialClipMask(:));

info.tap_quantized_element_count = ...
    sum(tapQuantizationMask(:));

info.qc_integerized_element_count = ...
    sum(qcIntegerMask(:));

info.changed_element_count = ...
    sum(finalChangeMask(:));

info.changed_solution_count = ...
    sum(any(finalChangeMask,2));

info.unchanged_solution_count = ...
    numberOfSolutions-info.changed_solution_count;

info.changed_mask = finalChangeMask;

info.maximum_absolute_change = ...
    max(abs(XdecodedMatrix(:)-Xoriginal(:)));

%% 10. 恢复单列向量输入方向

if inputWasColumnVector
    Xdecoded = XdecodedMatrix';
else
    Xdecoded = XdecodedMatrix;
end

end

