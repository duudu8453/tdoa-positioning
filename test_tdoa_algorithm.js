// 时差定位算法测试代码
// Time Difference of Arrival (TDOA) Positioning Algorithm Test

class TDOAPositioning {
    constructor() {
        this.c = 299792458; // 光速 m/s
        this.stations = [];
        this.errors = {
            timeError: 0,
            positionError: 0,
            atmosphericDelay: 0,
            multipathEffect: 0
        };
    }

    // 设置接收站位置
    setStations(stations) {
        this.stations = stations.map(s => ({...s}));
    }

    // 设置误差参数
    setErrors(errors) {
        this.errors = {...this.errors, ...errors};
    }

    // 计算两点间距离
    distance(p1, p2) {
        return Math.sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2);
    }

    // 添加高斯噪声
    addGaussianNoise(value, sigma) {
        const u1 = Math.random();
        const u2 = Math.random();
        const z0 = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
        return value + z0 * sigma;
    }

    // 模拟信号传播时间测量
    measurePropagationTimes(targetPosition) {
        const trueTimes = this.stations.map(station => 
            this.distance(station, targetPosition) / this.c
        );

        // 添加各种误差
        const noisyTimes = trueTimes.map((time, index) => {
            let noisyTime = time;
            
            // 时间测量误差
            noisyTime = this.addGaussianNoise(noisyTime, this.errors.timeError * 1e-6);
            
            // 大气延迟误差
            noisyTime += this.errors.atmosphericDelay * 1e-6;
            
            // 多径效应
            const multipathDelay = this.errors.multipathEffect * 
                (Math.random() - 0.5) * 0.1e-6;
            noisyTime += multipathDelay;
            
            return noisyTime;
        });

        return {
            trueTimes,
            noisyTimes,
            timeDifferences: this.calculateTimeDifferences(noisyTimes)
        };
    }

    // 计算时间差
    calculateTimeDifferences(times) {
        const diffs = {};
        for (let i = 0; i < this.stations.length; i++) {
            for (let j = i + 1; j < this.stations.length; j++) {
                const key = `${this.stations[i].id}-${this.stations[j].id}`;
                diffs[key] = times[i] - times[j];
            }
        }
        return diffs;
    }

    // TDOA定位算法 - 最小二乘法
    positionLeastSquares(timeDifferences) {
        if (this.stations.length < 3) {
            throw new Error('至少需要3个接收站');
        }

        const stations = this.stations;
        const c = this.c;
        
        // 将时间差转换为距离差
        const distanceDiffs = {};
        Object.entries(timeDifferences).forEach(([key, timeDiff]) => {
            distanceDiffs[key] = timeDiff * c;
        });

        // 使用迭代最小二乘法求解TDOA方程
        return this.iterativeLeastSquares(distanceDiffs);
    }

    // 迭代最小二乘法求解TDOA
    iterativeLeastSquares(distanceDiffs) {
        const stations = this.stations;
        
        // 初始估计位置（使用几何中心）
        let x = stations.reduce((sum, s) => sum + s.x, 0) / stations.length;
        let y = stations.reduce((sum, s) => sum + s.y, 0) / stations.length;
        
        const maxIterations = 10;
        const tolerance = 1e-6;
        
        for (let iter = 0; iter < maxIterations; iter++) {
            const A = [];
            const b = [];
            
            // 使用第一个站点作为参考
            const refStation = stations[0];
            const r0 = this.distance(refStation, {x, y});
            
            for (let i = 1; i < stations.length; i++) {
                const station = stations[i];
                const key = `${refStation.id}-${station.id}`;
                const rangeDiff = distanceDiffs[key] || 0;
                
                const ri = this.distance(station, {x, y});
                
                // 线性化的TDOA方程系数
                const a_row = [
                    (x - refStation.x) / r0 - (x - station.x) / ri,
                    (y - refStation.y) / r0 - (y - station.y) / ri
                ];
                
                // 残差
                const residual = r0 - ri - rangeDiff;
                
                A.push(a_row);
                b.push(-residual);
            }
            
            // 求解修正量
            try {
                const correction = this.solveLinearSystem(A, b);
                
                // 更新位置估计
                const newX = x + correction.x;
                const newY = y + correction.y;
                
                // 检查收敛性
                const deltaX = Math.abs(newX - x);
                const deltaY = Math.abs(newY - y);
                
                x = newX;
                y = newY;
                
                if (deltaX < tolerance && deltaY < tolerance) {
                    break;
                }
                
            } catch (e) {
                // 如果迭代失败，尝试简化的线性解法
                return this.simplifiedLinearSolution(distanceDiffs);
            }
        }
        
        return {x, y};
    }
    
    // 简化的线性解法（作为备选方案）
    simplifiedLinearSolution(distanceDiffs) {
        const stations = this.stations;
        const A = [];
        const b = [];
        
        // 使用第一个站点作为参考
        const refStation = stations[0];
        
        for (let i = 1; i < stations.length; i++) {
            const station = stations[i];
            const key = `${refStation.id}-${station.id}`;
            const rangeDiff = distanceDiffs[key] || 0;
            
            // 简化的线性方程（假设目标距离参考站较近）
            const a_row = [
                2 * (station.x - refStation.x),
                2 * (station.y - refStation.y)
            ];
            
            const b_val = (station.x ** 2 + station.y ** 2) - 
                         (refStation.x ** 2 + refStation.y ** 2) - 
                         2 * rangeDiff * this.distance(refStation, {
                             x: (station.x + refStation.x) / 2,
                             y: (station.y + refStation.y) / 2
                         });
            
            A.push(a_row);
            b.push(b_val);
        }

        return this.solveLinearSystem(A, b);
    }