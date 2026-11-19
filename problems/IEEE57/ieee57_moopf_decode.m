function [Xdecoded,info] = ieee57_moopf_decode(X,cfg)
%IEEE57_MOOPF_DECODE IEEE57 V2决策变量统一解码。
% PG/VG保持连续；TAP按0.02投影；QC按1 MVAr整数化。

if nargin<2 || isempty(cfg), cfg = ieee57_cdmtmo_config(); end
if ~isstruct(cfg) || ~isfield(cfg,'D') || ~isfield(cfg,'variable') || ~isfield(cfg,'index')
    error('cfg不是有效的IEEE57配置。');
end
if ~isnumeric(X) || ~isreal(X) || isempty(X) || any(~isfinite(X(:)))
    error('X必须是非空有限实数数组。');
end
X = double(X);
inputWasColumnVector = false;
if isvector(X)
    if numel(X)~=cfg.D, error('决策向量长度必须为%d。',cfg.D); end
    inputWasColumnVector = size(X,1)==cfg.D && size(X,2)==1;
    Xmatrix = reshape(X,1,cfg.D);
else
    if size(X,2)~=cfg.D, error('决策矩阵列数必须为%d。',cfg.D); end
    Xmatrix = X;
end

lower = reshape(cfg.variable.lower,1,cfg.D);
upper = reshape(cfg.variable.upper,1,cfg.D);
Xoriginal = Xmatrix;
Xbounded = min(max(Xmatrix,lower),upper);
initialClipMask = abs(Xbounded-Xoriginal)>1e-12;
XdecodedMatrix = Xbounded;

for d = cfg.index.TAP
    step = cfg.variable.step(d);
    XdecodedMatrix(:,d) = lower(d)+round((XdecodedMatrix(:,d)-lower(d))./step).*step;
    XdecodedMatrix(:,d) = round(XdecodedMatrix(:,d),12);
end

tapMask = false(size(XdecodedMatrix));
tapMask(:,cfg.index.TAP) = abs(XdecodedMatrix(:,cfg.index.TAP)-Xbounded(:,cfg.index.TAP))>1e-12;
qcBefore = XdecodedMatrix(:,cfg.index.QC);
XdecodedMatrix(:,cfg.index.QC) = round(XdecodedMatrix(:,cfg.index.QC));
qcMask = false(size(XdecodedMatrix));
qcMask(:,cfg.index.QC) = abs(XdecodedMatrix(:,cfg.index.QC)-qcBefore)>1e-12;
XdecodedMatrix = min(max(XdecodedMatrix,lower),upper);

if any(~isfinite(XdecodedMatrix(:))) || any(any(XdecodedMatrix<lower-1e-12)) || any(any(XdecodedMatrix>upper+1e-12))
    error('解码结果无效。');
end
for d = cfg.index.TAP
    step = cfg.variable.step(d);
    coordinate = (XdecodedMatrix(:,d)-lower(d))./step;
    if any(abs(coordinate-round(coordinate))>1e-9), error('TAP变量没有落在离散网格上。'); end
end
if any(any(abs(XdecodedMatrix(:,cfg.index.QC)-round(XdecodedMatrix(:,cfg.index.QC)))>1e-12))
    error('QC整数化失败。');
end

changed = abs(XdecodedMatrix-Xoriginal)>1e-12;
info.number_of_solutions = size(Xmatrix,1);
info.number_of_variables = cfg.D;
info.clipped_element_count = sum(initialClipMask(:));
info.tap_quantized_element_count = sum(tapMask(:));
info.qc_integerized_element_count = sum(qcMask(:));
info.changed_element_count = sum(changed(:));
info.changed_solution_count = sum(any(changed,2));
info.unchanged_solution_count = size(Xmatrix,1)-info.changed_solution_count;
info.changed_mask = changed;
info.maximum_absolute_change = max(abs(XdecodedMatrix(:)-Xoriginal(:)));

if inputWasColumnVector, Xdecoded = XdecodedMatrix'; else, Xdecoded = XdecodedMatrix; end
end
