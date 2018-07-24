# DEsingle_subsampling pre-experiment
## 目前功能
- 生成四种基因的数据，并执行`DEsingle`，将结果保存在`./result.csv`中
- 将用于生成数据的原始参数以及`DEsingle`生成的结果保存于`./origin_data.csv`中，可以直接通过分析R对象`record`来分析`DEsingle`的准确率了
- 输出`LoU-细胞数量`折线图`graph.png`
## TODO
- 修改随机生成四种基因数据的参数，目前生成`DEa`很有问题，生成的`DEa`数据**完全**不能被`DEsingle`发现
- __根据生成零膨胀负二项分布的两组原始参数来直接判断差异表达类型，以保证数据生成的正确性__
## 用法
修改`Script.R`中的参数，然后在R终端中执行
```
source("./Script.R")
```
