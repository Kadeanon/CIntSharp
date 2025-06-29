这是一个专为量子化学所需向量化操作设计的简单双精度数学库，包括常规线性代数和简单高维代数，旨在提供较高的性能与互操作性。

其中线性代数部分使用pure-C# 模板化技术编写1v/m/f/2m函数，对3m与lapack例程将直接使用mkl.net作为计算后端以提供高级功能。涉及功能均为按需编写（以及部分移植自以前项目的代码），保证主库最低限度要求。

高维代数部分仅实现了通用遍历方法与双张量收缩，尚未实现广播和reduce相关函数，reduce暂时依靠System.Numeric.Tensors库实现，但由于SNT库对维度步长的限制，有时需要手动进行转置。张量收缩采用TBLIS使用的bsmtc算法，未来将会引入batchmm/TTGT作为补充并实现任意数量张量收缩。

由于实现时间较短，且项目跨度较大，大量方法未经测试，可能出现编码风格不统一、边缘状况异常或者实现错误，欢迎随意纠错和提出建议。

### 使用的库：

线性代数后端：[MKL.NET ](https://github.com/MKL-NET/MKL.NET)

基础张量支持：[System.Numerics.Tensors](https://github.com/dotnet/runtime)

ufunc模板化库：[NetFabric.Numerics.Tensors](https://github.com/NetFabric/NetFabric.Numerics.Tensors)

其它辅助库：

[Microsoft.Extensions.ObjectPool](https://github.com/dotnet/runtime)

[CommunityToolkit.HighPerformance](https://github.com/CommunityToolkit/dotnet)

[ExSpans](https://github.com/zyl910/ExSpans)

### 项目参考

BLIS: [flame/blis: BLAS-like Library Instantiation Software Framework](https://github.com/flame/blis)

TBLIS: [MatthewsResearchGroup/tblis: TBLIS is a library and framework for performing tensor operations, especially tensor contraction, using efficient native algorithms.](https://github.com/MatthewsResearchGroup/tblis)