// TDOA定位算法实现

// 设置接收站位置和参数
function setupReceivers(receivers) {
    return receivers.map(r => ({
        x: r.x,
        y: r.y,
        timeError: r.timeError || 0,  // 时间测量误差
        positionError: r.positionError || 0  // 位置误差
    }));
}

// 设置误差参数
function setErrorParameters(params) {
    return {
        atmosphericDelay: params.atmosphericDelay || 0,  // 大气延迟
        multipathEffect: params.multipathEffect || 0,    // 多径效应
        timeError: params.timeError || 0,                // 时间测量误差
        positionError: params.positionError || 0         // 位置误差
    };
}

// 计算实际时间差（包含误差）
function calculateTimeDifferences(target, receivers, errorParams) {
    const timeDiffs = [];
    const c = 299792458; // 光速（米/秒）

    for (let i = 1; i < receivers.length; i++) {
        // 计算真实距离
        const d1 = Math.sqrt(Math.pow(target.x - receivers[0].x, 2) + 
                            Math.pow(target.y - receivers[0].y, 2));
        const d2 = Math.sqrt(Math.pow(target.x - receivers[i].x, 2) + 
                            Math.pow(target.y - receivers[i].y, 2));

        // 添加误差影响
        const timeError = (Math.random() - 0.5) * 2 * errorParams.timeError;
        const posError1 = (Math.random() - 0.5) * 2 * errorParams.positionError;
        const posError2 = (Math.random() - 0.5) * 2 * errorParams.positionError;
        const atmDelay = Math.random() * errorParams.atmosphericDelay;
        const multipath = Math.random() * errorParams.multipathEffect;

        // 计算总的时间差
        const timeDiff = (d2 - d1) / c + timeError + 
                        (posError2 - posError1) / c + 
                        atmDelay + multipath;

        timeDiffs.push(timeDiff);
    }

    return timeDiffs;
}

// 计算理想时间差（无误差）
function calculateIdealTimeDifferences(target, receivers) {
    const timeDiffs = [];
    const c = 299792458; // 光速

    for (let i = 1; i < receivers.length; i++) {
        const d1 = Math.sqrt(Math.pow(target.x - receivers[0].x, 2) + 
                            Math.pow(target.y - receivers[0].y, 2));
        const d2 = Math.sqrt(Math.pow(target.x - receivers[i].x, 2) + 
                            Math.pow(target.y - receivers[i].y, 2));
        timeDiffs.push((d2 - d1) / c);
    }

    return timeDiffs;
}

// 使用迭代最小二乘法求解TDOA定位
function solveTDOAPosition(timeDiffs, receivers, initialGuess) {
    const c = 299792458; // 光速
    const maxIterations = 100;
    const tolerance = 1e-6;
    
    let x = initialGuess.x;
    let y = initialGuess.y;

    for (let iter = 0; iter < maxIterations; iter++) {
        // 构建雅可比矩阵和残差向量
        const J = [];
        const r = [];
        const ref = receivers[0];

        for (let i = 1; i < receivers.length; i++) {
            const d1 = Math.sqrt(Math.pow(x - ref.x, 2) + Math.pow(y - ref.y, 2));
            const d2 = Math.sqrt(Math.pow(x - receivers[i].x, 2) + 
                                Math.pow(y - receivers[i].y, 2));

            // 计算雅可比矩阵元素
            const dx1 = (x - ref.x) / d1;
            const dy1 = (y - ref.y) / d1;
            const dx2 = (x - receivers[i].x) / d2;
            const dy2 = (y - receivers[i].y) / d2;

            J.push([dx2 - dx1, dy2 - dy1]);
            r.push(timeDiffs[i-1] * c - (d2 - d1));
        }

        // 求解正规方程
        const JtJ = [[0, 0], [0, 0]];
        const Jtr = [0, 0];

        for (let i = 0; i < J.length; i++) {
            JtJ[0][0] += J[i][0] * J[i][0];
            JtJ[0][1] += J[i][0] * J[i][1];
            JtJ[1][0] += J[i][1] * J[i][0];
            JtJ[1][1] += J[i][1] * J[i][1];

            Jtr[0] += J[i][0] * r[i];
            Jtr[1] += J[i][1] * r[i];
        }

        // 计算位置更新
        const det = JtJ[0][0] * JtJ[1][1] - JtJ[0][1] * JtJ[1][0];
        const dx = (JtJ[1][1] * Jtr[0] - JtJ[0][1] * Jtr[1]) / det;
        const dy = (-JtJ[1][0] * Jtr[0] + JtJ[0][0] * Jtr[1]) / det;

        // 更新位置估计
        x += dx;
        y += dy;

        // 检查收敛性
        if (Math.sqrt(dx*dx + dy*dy) < tolerance) {
            break;
        }
    }

    return { x, y };
}

// 导出函数
module.exports = {
    setupReceivers,
    setErrorParameters,
    calculateTimeDifferences,
    calculateIdealTimeDifferences,
    solveTDOAPosition
};