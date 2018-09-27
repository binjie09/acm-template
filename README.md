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

```c++
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

## ST表
用于快速查询RMQ问题 建表时间复杂度O(log(n))空间复杂度O(log(n)) 查询时间复杂度O(1) 不可动态更改

```c++
/*
  Version: 1.0
  Author: 王峰
  Date: 2018.8.27
  Function List:
	init_table(n) n为原数组a的长度
	query_table(L, R) 查询原数组内区间[L,R]的最值
*/
#include <cstdio>
#include <algorithm>

using std::min;
const int N = 1000010;

int table[N][100];
int lg[N];
int a[N];

void init_table(int n) {
	for(int i=0; i<n; i++) table[i][0]=a[i];
	for(int j=1; (1<<j)<=n; j++) {
		for(int i=0; i+(1<<j)-1<n; i++) {
			table[i][j] = min(table[i][j-1],table[i+(1<<(j-1))][j-1]);
		}
	} 
	for(int i=2, k=0; i<=n; i++) {
		if(i>(1<<k)) k++;
		lg[i] = k-1;
	}
}

int query_table(int L, int R) {
	int k = lg[R-L+1];
	return min(table[L][k],table[R-(1<<k)+1][k]); 
}
```

## 树状数组（二叉索引树）模板
用于求可单点更新的前缀和，时间复杂度O(log(n))
```c++
/*
  Version: 1.0
  Author: 王峰
  Date: 2018.8.27
  Function List:
	add(k, v) 将数组中的第k个元素的值+v，注意k从1开始计数
	query(k) 查询profix[k]
*/
#include <cstdio>

const int N = 1000010;
int n; //实际数组的长度 
long long tree[4*N];

void add(int k, long long v) {
	while(k<n){
		tree[k]+=v;
		k+=(k&-k);	
	}
}
long long query(int k) {
	long long ans = 0;
	while(k>0) {
		ans+=tree[k];
		k-=k&-k;
	}
	return ans;
}
```


## dijkstra+priority_queue

用于求单源最短路径，不能处理负环
```c++
/*
  Version: 1.0
  Author: 王峰
  Date: 2018.9.3
  Function List:
	init(int n) 将带n个点的图初始化
	addEdge(int from, int to, int dist) 增加一条单向边
	dijkstra(int s) 求s为源点的最短路
	print(t, s) 打印从s为源点t为终点的一条路径 
*/
#include <cstdio>
#include <vector>
#include <queue>
#include <cstring>

using std::vector;
using std::priority_queue;
const int maxn = 1e6+7;
const int INF = 0x3f3f3f3f;

struct Edge {
	int from, to, dist;
};

struct HeapNode {
	int d, u;
	bool operator < (const HeapNode& rhs) const {
		return d > rhs.d;
	}
}; 

struct Dijkstra {
	int n, m;                //点数和边数 
	vector<Edge> edges;      //边列表 
	vector<int> G[maxn];     //每个节点出发的边编号， 从0开始编号 
	bool vis[maxn];          //是否已经标记 
	int d[maxn];             //s到各个点的距离 
	//int p[maxn];           //最短路的上一条边 
	
	void init(int n) {
		this->n = n;
		for(int i=0; i<=n; i++) {
			G[i].clear();
		}
		edges.clear();
	}
	
	/*
	void print(int t, int s) {
		int pre = edges[p[t]].from;
		if(pre==s) {
			printf("%d ", s);
			return ;
		}else{
			print(pre,s);
			printf("%d ", pre);	
		} 
	}
	*/
	
	void addEdge(int from, int to, int dist) {
		edges.push_back((Edge){from, to, dist});
		m = edges.size();
		G[from].push_back(m-1);
	}
	
	void dijkstra(int s) {
		priority_queue<HeapNode> Q;
		for(int i=0; i<=n; i++) d[i] = INF;
		d[s]=0;
		memset(vis, 0, sizeof(vis));
		Q.push((HeapNode){0, s});
		while(!Q.empty()) {
			HeapNode x = Q.top(); Q.pop();
			int u = x.u;
			if(vis[u]) continue;
			vis[u] = true;
			for(int i=0; i<G[u].size(); i++) {
				Edge & e = edges[G[u][i]];
				if(d[e.to] > d[u] + e.dist) {
					d[e.to] = d[u] + e.dist;
					//p[e.to] = G[u][i];
					Q.push((HeapNode){d[e.to], e.to});
				}
			}
		}
	}
}D;
```

## 求和线段树模板（RMQ线段树小改即可）


```c++
/*
Function List: 
    Build(1,n,1); 建树 
    Add(L,C,1,n,1); 点更新 
    Update(L,R,C,1,n,1); 区间更新 
    ANS=Query(L,R,1,n,1); 区间查询 
*/

typedef long long LL;

const int MAXN=1e5+7;
LL a[MAXN],ans[MAXN<<2],lazy[MAXN<<2];

void PushUp(int rt) {
    ans[rt]=ans[rt<<1]+ans[rt<<1|1];
}

void Build(int l,int r,int rt) {
    if (l==r) {
        ans[rt]=0;
        return;
    }
    int mid=(l+r)>>1;
    Build(l,mid,rt<<1);
    Build(mid+1,r,rt<<1|1);
    PushUp(rt);
}

void PushDown(int rt,LL ln,LL rn) {
    if (lazy[rt]) {
        lazy[rt<<1]+=lazy[rt];
        lazy[rt<<1|1]+=lazy[rt];
        ans[rt<<1]+=lazy[rt]*ln;
        ans[rt<<1|1]+=lazy[rt]*rn;
        lazy[rt]=0;
    }
}

void Add(int L,int C,int l,int r,int rt)
{
    if (l==r)
    {
        ans[rt]+=C;
        return;
    }
    int mid=(l+r)>>1;
    //PushDown(rt,mid-l+1,r-mid); 若既有点更新又有区间更新，需要这句话
    if (L<=mid)
        Add(L,C,l,mid,rt<<1);
    else
        Add(L,C,mid+1,r,rt<<1|1);
    PushUp(rt);
}

void Update(int L,int R,long long C,int l,int r,int rt) {
    if (L<=l&&r<=R) {
        ans[rt]+=C*(r-l+1);
        lazy[rt]+=C;
        return;
    }
    int mid=(l+r)>>1;
    PushDown(rt,mid-l+1,r-mid);
    if (L<=mid) Update(L,R,C,l,mid,rt<<1);
    if (R>mid) Update(L,R,C,mid+1,r,rt<<1|1);
    PushUp(rt);
}

LL Query(int L,int R,int l,int r,int rt) {
    if (L<=l&&r<=R)
        return ans[rt];
    int mid=(l+r)>>1;
    PushDown(rt,mid-l+1,r-mid);
    LL ANS=0;
    if (L<=mid) ANS+=Query(L,R,l,mid,rt<<1);
    if (R>mid) ANS+=Query(L,R,mid+1,r,rt<<1|1);
    return ANS;
}
```

## C组合数板子

时间复杂度O(2^n)

```c++
/*
  Version: 1.0
  Author: 王峰
  Date: 2018.9.17
*/
C[1][0] = C[1][1] = 1;
for (int i = 2; i < N; i++){
	C[i][0] = 1;
	for (int j = 1; j < N; j++)
		C[i][j] = (C[i - 1][j] + C[i - 1][j - 1]);
}
```

## 容斥原理

时间复杂度O(n^2)

```c++
ll ans = 0;
for(int i=0; i<(1<<tot); i++) {

	int cnt = 0;
	for(int j=0; j<tot; j++) {
		if(i&(1<<j)) cnt++;
    }

    ll tmp = 0;

        /*
            得到当前状态下的贡献值
        */

    if(cnt&1) {
    	ans ++ tmp;
    } else {
    	ans -= tmp;
    }
}
```

## 主席树

### 主席树：可持久化线段树,
### 在nlog(n)空内维护多个权值线段树
### 离线数据结构不支持修改

## 权值线段树模板未更新


```c++
#include <cstdio>
#include <vector>
#include <algorithm>



/*  
  Version: 1.0
  Author: 王峰
  Date: 2018.9.25
  Function List:
  本例为查询区间第K大,有一个基于排序去重的哈希,有时间分离出来单独成板.
*/

using std::lower_bound;
using std::vector;

const int maxn = 1e5+6;
int cnt, root[maxn], a[maxn];
struct node{
	int l, r, sum;
}T[maxn*40]; 
vector<int > v;

int getid(int x) {
	return lower_bound(v.begin(), v.end(),x)-v.begin()+1;
}

void update(int l, int r, int &x, int y, int pos) {
	T[++cnt]=T[y],T[cnt].sum++,x=cnt;
	if(l==r)return ;
	int mid=(l+r)/2;
	if(mid>=pos)update(l, mid, T[x].l, T[y].l, pos);
	else update(mid+1, r, T[x].r, T[y].r, pos);
}

int query(int l, int r, int x, int y, int k) {
	if(l==r) return l;
	int mid = (l+r)/2;
	int sum = T[T[y].l].sum - T[T[x].l].sum;
	if(sum>=k)return query(l, mid, T[x].l, T[y].l, k);
	else return query(mid+1, r, T[x].r, T[y].r, k-sum);
} 

int main() {
	int n, q;
	scanf("%d%d", &n, &q);
	for(int i=1; i<=n; i++) {
		scanf("%d", &a[i]);
		v.push_back(a[i]);
	}
	sort(v.begin(), v.end()),v.erase(unique(v.begin(), v.end()), v.end());
	for(int i=1; i<=n; i++) update(1, n, root[i], root[i-1], getid(a[i]));
	for(int i=1; i<=q; i++) {
		int x, y, k; 
		scanf("%d%d%d", &x, &y, &k);
		printf("%d\n", v[query(1,n,root[x-1],root[y],k)-1]);
	}
	return 0;
}

```

## lucas定理

描述: 直接用就行。。。
```c++
/*  
  Version: 1.0
  Author: 陈彬杰
  Date: 2018.9.27
  Function List:
  Lucas(n, m, p): 返回 C(n,m) mod p （p为素数）的值。
*/
LL Lucas(LL n, LL m, int p)
{
  return m ? Lucas(n/p, m/p, p) * comb(n%p, m%p, p) % p : 1;
}
```

## some else algorithm

### some else
