OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0776032) q[0];
sx q[0];
rz(-0.912323) q[0];
sx q[0];
rz(1.3349226) q[0];
rz(-1.1420684) q[1];
sx q[1];
rz(-0.51841441) q[1];
sx q[1];
rz(-1.4919182) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27150422) q[0];
sx q[0];
rz(-0.056500204) q[0];
sx q[0];
rz(1.9046049) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9426192) q[2];
sx q[2];
rz(-0.97731263) q[2];
sx q[2];
rz(1.3243937) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.093203737) q[1];
sx q[1];
rz(-1.2902765) q[1];
sx q[1];
rz(-0.93757052) q[1];
rz(1.8133477) q[3];
sx q[3];
rz(-1.537916) q[3];
sx q[3];
rz(1.7612847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3492744) q[2];
sx q[2];
rz(-2.8742542) q[2];
sx q[2];
rz(0.054923687) q[2];
rz(-1.6867636) q[3];
sx q[3];
rz(-0.39595404) q[3];
sx q[3];
rz(-1.6247862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2199663) q[0];
sx q[0];
rz(-1.8333789) q[0];
sx q[0];
rz(0.24800214) q[0];
rz(1.8964881) q[1];
sx q[1];
rz(-0.73687941) q[1];
sx q[1];
rz(-1.2287593) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22511521) q[0];
sx q[0];
rz(-0.18813293) q[0];
sx q[0];
rz(0.1779025) q[0];
rz(-pi) q[1];
rz(-2.6284559) q[2];
sx q[2];
rz(-2.0433132) q[2];
sx q[2];
rz(3.0233011) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5692484) q[1];
sx q[1];
rz(-2.4048785) q[1];
sx q[1];
rz(-2.8235497) q[1];
x q[2];
rz(-1.5108438) q[3];
sx q[3];
rz(-1.151181) q[3];
sx q[3];
rz(2.7762716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2174786) q[2];
sx q[2];
rz(-2.3096297) q[2];
sx q[2];
rz(-2.5577616) q[2];
rz(-2.7122688) q[3];
sx q[3];
rz(-1.2238294) q[3];
sx q[3];
rz(-1.0113641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9871224) q[0];
sx q[0];
rz(-0.38726375) q[0];
sx q[0];
rz(-1.467147) q[0];
rz(1.4779429) q[1];
sx q[1];
rz(-1.4091622) q[1];
sx q[1];
rz(-1.0999701) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3008793) q[0];
sx q[0];
rz(-1.7689118) q[0];
sx q[0];
rz(1.4639356) q[0];
x q[1];
rz(-2.2713678) q[2];
sx q[2];
rz(-0.43284518) q[2];
sx q[2];
rz(-0.74898042) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63640672) q[1];
sx q[1];
rz(-2.6224394) q[1];
sx q[1];
rz(2.4207508) q[1];
x q[2];
rz(-2.1446225) q[3];
sx q[3];
rz(-1.1862635) q[3];
sx q[3];
rz(-1.9652733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.35525068) q[2];
sx q[2];
rz(-2.2189271) q[2];
sx q[2];
rz(1.5477017) q[2];
rz(-1.8111546) q[3];
sx q[3];
rz(-1.3733613) q[3];
sx q[3];
rz(1.0734585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7610157) q[0];
sx q[0];
rz(-1.8208193) q[0];
sx q[0];
rz(-2.9837578) q[0];
rz(2.8997968) q[1];
sx q[1];
rz(-2.7561185) q[1];
sx q[1];
rz(-2.6604624) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07040992) q[0];
sx q[0];
rz(-0.86127036) q[0];
sx q[0];
rz(-1.6874403) q[0];
x q[1];
rz(1.3161737) q[2];
sx q[2];
rz(-1.2224276) q[2];
sx q[2];
rz(0.35149945) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1222904) q[1];
sx q[1];
rz(-1.9112406) q[1];
sx q[1];
rz(-2.2346418) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9322196) q[3];
sx q[3];
rz(-1.2897769) q[3];
sx q[3];
rz(-3.0176153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.87159291) q[2];
sx q[2];
rz(-0.81120482) q[2];
sx q[2];
rz(-2.8672186) q[2];
rz(1.6756049) q[3];
sx q[3];
rz(-1.2803187) q[3];
sx q[3];
rz(-2.2082641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45193732) q[0];
sx q[0];
rz(-0.87757293) q[0];
sx q[0];
rz(1.0614606) q[0];
rz(2.6929216) q[1];
sx q[1];
rz(-2.6866388) q[1];
sx q[1];
rz(1.7015069) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3685303) q[0];
sx q[0];
rz(-0.8196747) q[0];
sx q[0];
rz(-1.7858265) q[0];
x q[1];
rz(3.0147047) q[2];
sx q[2];
rz(-1.972162) q[2];
sx q[2];
rz(0.021406476) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6430407) q[1];
sx q[1];
rz(-1.3128377) q[1];
sx q[1];
rz(-0.027679701) q[1];
x q[2];
rz(-3.0331221) q[3];
sx q[3];
rz(-2.1785695) q[3];
sx q[3];
rz(-2.8911107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7777286) q[2];
sx q[2];
rz(-1.2960459) q[2];
sx q[2];
rz(1.1931922) q[2];
rz(0.39919546) q[3];
sx q[3];
rz(-1.7292855) q[3];
sx q[3];
rz(3.1138368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6804009) q[0];
sx q[0];
rz(-1.2874648) q[0];
sx q[0];
rz(2.4611018) q[0];
rz(-0.83241278) q[1];
sx q[1];
rz(-1.2544371) q[1];
sx q[1];
rz(2.9852273) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.844961) q[0];
sx q[0];
rz(-0.28066844) q[0];
sx q[0];
rz(-0.190221) q[0];
rz(-2.2920798) q[2];
sx q[2];
rz(-0.88187432) q[2];
sx q[2];
rz(-0.22077862) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.029441) q[1];
sx q[1];
rz(-1.8990771) q[1];
sx q[1];
rz(-0.29354696) q[1];
x q[2];
rz(-0.12229192) q[3];
sx q[3];
rz(-1.3548791) q[3];
sx q[3];
rz(1.2014887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51123315) q[2];
sx q[2];
rz(-2.4801621) q[2];
sx q[2];
rz(2.2516001) q[2];
rz(0.47641274) q[3];
sx q[3];
rz(-0.92372957) q[3];
sx q[3];
rz(2.711003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7156859) q[0];
sx q[0];
rz(-1.8222734) q[0];
sx q[0];
rz(2.529378) q[0];
rz(-1.3587492) q[1];
sx q[1];
rz(-0.76018676) q[1];
sx q[1];
rz(-2.8403958) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33443794) q[0];
sx q[0];
rz(-2.8557192) q[0];
sx q[0];
rz(1.8854333) q[0];
x q[1];
rz(2.0400171) q[2];
sx q[2];
rz(-2.42832) q[2];
sx q[2];
rz(-0.81339806) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4674356) q[1];
sx q[1];
rz(-1.1493719) q[1];
sx q[1];
rz(-1.9860877) q[1];
rz(-pi) q[2];
rz(-2.4047818) q[3];
sx q[3];
rz(-2.7932146) q[3];
sx q[3];
rz(0.37955561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.23413868) q[2];
sx q[2];
rz(-2.335304) q[2];
sx q[2];
rz(-1.0478728) q[2];
rz(0.35510865) q[3];
sx q[3];
rz(-0.57098782) q[3];
sx q[3];
rz(2.4560438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.909914) q[0];
sx q[0];
rz(-2.7993918) q[0];
sx q[0];
rz(2.8177596) q[0];
rz(2.553885) q[1];
sx q[1];
rz(-0.67981845) q[1];
sx q[1];
rz(-0.30276611) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97070951) q[0];
sx q[0];
rz(-1.8237178) q[0];
sx q[0];
rz(-0.33527261) q[0];
x q[1];
rz(-1.5220171) q[2];
sx q[2];
rz(-1.3338727) q[2];
sx q[2];
rz(2.0120442) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.3564612) q[1];
sx q[1];
rz(-1.0279127) q[1];
sx q[1];
rz(-2.7144542) q[1];
rz(-pi) q[2];
x q[2];
rz(1.15074) q[3];
sx q[3];
rz(-1.3358634) q[3];
sx q[3];
rz(2.3714424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6451463) q[2];
sx q[2];
rz(-2.1108997) q[2];
sx q[2];
rz(-2.5808064) q[2];
rz(-2.136611) q[3];
sx q[3];
rz(-1.2882261) q[3];
sx q[3];
rz(1.0952449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1110558) q[0];
sx q[0];
rz(-2.472214) q[0];
sx q[0];
rz(2.3916767) q[0];
rz(2.7610682) q[1];
sx q[1];
rz(-2.1546202) q[1];
sx q[1];
rz(-2.1662625) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.014054) q[0];
sx q[0];
rz(-0.97146266) q[0];
sx q[0];
rz(-2.4571193) q[0];
x q[1];
rz(1.8504233) q[2];
sx q[2];
rz(-1.6284918) q[2];
sx q[2];
rz(1.9967494) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5228032) q[1];
sx q[1];
rz(-2.0674043) q[1];
sx q[1];
rz(1.9267531) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9191429) q[3];
sx q[3];
rz(-1.2764023) q[3];
sx q[3];
rz(2.6840212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.19899496) q[2];
sx q[2];
rz(-0.90513217) q[2];
sx q[2];
rz(1.0814166) q[2];
rz(-3.0723451) q[3];
sx q[3];
rz(-1.705575) q[3];
sx q[3];
rz(-0.92250219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0104495) q[0];
sx q[0];
rz(-0.32587019) q[0];
sx q[0];
rz(0.3057873) q[0];
rz(2.8336613) q[1];
sx q[1];
rz(-1.4762069) q[1];
sx q[1];
rz(1.9889132) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.688153) q[0];
sx q[0];
rz(-1.3092293) q[0];
sx q[0];
rz(-0.75193172) q[0];
rz(-pi) q[1];
rz(-2.9342317) q[2];
sx q[2];
rz(-0.93610686) q[2];
sx q[2];
rz(0.35516741) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.252723) q[1];
sx q[1];
rz(-2.6111341) q[1];
sx q[1];
rz(0.88900755) q[1];
x q[2];
rz(3.1265101) q[3];
sx q[3];
rz(-1.6341356) q[3];
sx q[3];
rz(-0.36288211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5311188) q[2];
sx q[2];
rz(-1.0627154) q[2];
sx q[2];
rz(1.6884241) q[2];
rz(0.48062634) q[3];
sx q[3];
rz(-2.5745945) q[3];
sx q[3];
rz(-1.3948729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4881445) q[0];
sx q[0];
rz(-1.6076037) q[0];
sx q[0];
rz(-1.6758767) q[0];
rz(-2.0987971) q[1];
sx q[1];
rz(-1.4029618) q[1];
sx q[1];
rz(-0.89697368) q[1];
rz(-0.43389141) q[2];
sx q[2];
rz(-2.5627651) q[2];
sx q[2];
rz(2.4252575) q[2];
rz(3.0328209) q[3];
sx q[3];
rz(-2.7815007) q[3];
sx q[3];
rz(-0.23063767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
