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

## O(log(n)) 素数判定
目前已知的最快的单个素数判定方法

```c++
/*
  Version: 1.0
  Author: 王峰
  Date: 2018.8.26
  Function List:
	isPrime(long long n)
  返回值为n是否是质数
*/
#include <iostream> 
#include <cmath>

long long pow_mod(int a, long long d, long long mod) {
	long long res = 1, k = a;
	while(d) {
		if(d & 1) res = (res * k) % mod; 
		k = k * k % mod;
		d >>= 1;
	} 
	return res;
}
bool test(long long n, int a, long long d) {
	if(n == 2) return true;
	if(n == a) return true;
	if((n & 1) == 0) return false;
	while(!(d & 1)) d = d >> 1;
	long long t = pow_mod(a, d, n);
	while((d != n - 1) && (t != 1) && (t != n - 1)) {
		t = t * t % n;
		d <<= 1; 
	}
	return (t == n - 1 || (d & 1) == 1);
}
bool isPrime(long long n) {
	if(n < 2LL) return false;
	int a[] = {2, 3, 61};
	for (int i = 0; i <= 2; ++i) if(!test(n, a[i], n - 1)) return false;
	return true; 
} 
```

## 欧拉筛法求素数表
时间复杂度 O(n)

```
/*
  Version: 1.0
  Author: 王峰
  Date: 2018.8.26
  Function List:
	get_prime()
  check[i]表示i是否为合数
  prime[i]表示第i个质数,i从0开始计数
*/
#include<cstdio>
#include<cstring>

const int MAXN = 10000001;
const int MAXL = 10000001;
int prime[MAXN];
int check[MAXL];

int get_prime() {
	int tot = 0;
	memset(check, 0, sizeof(check));
	for (int i = 2; i < MAXL; ++i) {
		if (!check[i]){
	    	prime[tot++] = i;
		}
		for (int j = 0; j < tot; ++j) {
	    	if (i * prime[j] > MAXL) {
	    		break;
			}
	    	check[i*prime[j]] = 1;
	    	if (i % prime[j] == 0) {
	      		break;
			}
		}
	}
	return tot;
}
```

## some else algorithm

### some else
