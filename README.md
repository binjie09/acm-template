# acm-template

## 模板编写规范说明

- 使用markdown语法（便于赛前打印和生成目录）
- 所有的算法模板都写在这一个md文件中
- 代码名称
- 需要有接口说明（可以看下面示例）
- 代码最小原则，不添加任何不必要的变量和函数（最小不意味着最少，不要为了减短代码而去减短代码）
- 代码块要使用如下的markdown包裹
```
    ```c++
    your code  
    ```
```
- 其他细则在实践过程中添加

## 快速幂

简介:介绍算法干啥，哪些情况用的上  


### 非递归版算法

描述:简要描述算法实现框架（便于找错）。

```c++
/*
  Version: 1.0
  Author: 徐祥昊
  Date: 2018.7.27
  Function List:
    1.q_pow(__int64 a,__int64 b)  
      a:这里是a参数的说明;  
      b:这里是b参数的说明;
      ret:这里是返回值的说明;
    2.some else func
    3.some else func
*/
bool q_pow(__int64 a,__int64 b)
{
	__int64 res=1,mod=b,carry=a;
	while(b)
	{
		if(b%2)
		{
			res=res*a%mod;
		}
		a=a*a%mod;
		b>>=1;
	}
	if(res == carry)return true;
	else return false;
}
```

## some else algorithm

### some else