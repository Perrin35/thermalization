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
rz(0.063989446) q[0];
sx q[0];
rz(-2.2292697) q[0];
sx q[0];
rz(1.8066701) q[0];
rz(-1.1420684) q[1];
sx q[1];
rz(-0.51841441) q[1];
sx q[1];
rz(-1.4919182) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60580743) q[0];
sx q[0];
rz(-1.6241747) q[0];
sx q[0];
rz(3.1230631) q[0];
x q[1];
rz(-0.69500287) q[2];
sx q[2];
rz(-2.0796513) q[2];
sx q[2];
rz(-2.5093504) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.277093) q[1];
sx q[1];
rz(-0.96597176) q[1];
sx q[1];
rz(0.34326174) q[1];
x q[2];
rz(-1.434695) q[3];
sx q[3];
rz(-2.8968662) q[3];
sx q[3];
rz(-0.32258209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7923183) q[2];
sx q[2];
rz(-0.26733843) q[2];
sx q[2];
rz(-3.086669) q[2];
rz(-1.4548291) q[3];
sx q[3];
rz(-0.39595404) q[3];
sx q[3];
rz(1.6247862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216264) q[0];
sx q[0];
rz(-1.3082137) q[0];
sx q[0];
rz(2.8935905) q[0];
rz(-1.2451046) q[1];
sx q[1];
rz(-2.4047132) q[1];
sx q[1];
rz(1.2287593) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9707391) q[0];
sx q[0];
rz(-1.6038994) q[0];
sx q[0];
rz(2.956361) q[0];
rz(-2.3360148) q[2];
sx q[2];
rz(-0.68289872) q[2];
sx q[2];
rz(-0.77308347) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5723443) q[1];
sx q[1];
rz(-0.73671417) q[1];
sx q[1];
rz(-0.31804292) q[1];
rz(3.008083) q[3];
sx q[3];
rz(-0.4236246) q[3];
sx q[3];
rz(2.6299919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2174786) q[2];
sx q[2];
rz(-2.3096297) q[2];
sx q[2];
rz(0.58383101) q[2];
rz(-2.7122688) q[3];
sx q[3];
rz(-1.9177633) q[3];
sx q[3];
rz(-2.1302285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15447021) q[0];
sx q[0];
rz(-2.7543289) q[0];
sx q[0];
rz(-1.467147) q[0];
rz(1.6636498) q[1];
sx q[1];
rz(-1.7324305) q[1];
sx q[1];
rz(-1.0999701) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80187782) q[0];
sx q[0];
rz(-2.9168282) q[0];
sx q[0];
rz(0.4883713) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8520865) q[2];
sx q[2];
rz(-1.2443674) q[2];
sx q[2];
rz(-3.141186) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.63640672) q[1];
sx q[1];
rz(-2.6224394) q[1];
sx q[1];
rz(0.7208419) q[1];
x q[2];
rz(2.2113924) q[3];
sx q[3];
rz(-0.67852321) q[3];
sx q[3];
rz(0.92032209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35525068) q[2];
sx q[2];
rz(-2.2189271) q[2];
sx q[2];
rz(1.5477017) q[2];
rz(-1.330438) q[3];
sx q[3];
rz(-1.3733613) q[3];
sx q[3];
rz(2.0681341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.380577) q[0];
sx q[0];
rz(-1.3207734) q[0];
sx q[0];
rz(2.9837578) q[0];
rz(-0.24179587) q[1];
sx q[1];
rz(-2.7561185) q[1];
sx q[1];
rz(-2.6604624) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0711827) q[0];
sx q[0];
rz(-0.86127036) q[0];
sx q[0];
rz(1.6874403) q[0];
x q[1];
rz(1.8254189) q[2];
sx q[2];
rz(-1.2224276) q[2];
sx q[2];
rz(-0.35149945) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.85192902) q[1];
sx q[1];
rz(-0.73411513) q[1];
sx q[1];
rz(-2.0925702) q[1];
rz(-pi) q[2];
rz(2.9322196) q[3];
sx q[3];
rz(-1.8518157) q[3];
sx q[3];
rz(-3.0176153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.87159291) q[2];
sx q[2];
rz(-2.3303878) q[2];
sx q[2];
rz(2.8672186) q[2];
rz(1.6756049) q[3];
sx q[3];
rz(-1.2803187) q[3];
sx q[3];
rz(0.93332851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45193732) q[0];
sx q[0];
rz(-0.87757293) q[0];
sx q[0];
rz(-2.080132) q[0];
rz(-0.44867107) q[1];
sx q[1];
rz(-2.6866388) q[1];
sx q[1];
rz(1.7015069) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3502304) q[0];
sx q[0];
rz(-1.7273979) q[0];
sx q[0];
rz(2.3788405) q[0];
rz(-pi) q[1];
x q[1];
rz(1.975073) q[2];
sx q[2];
rz(-1.4540404) q[2];
sx q[2];
rz(1.5991885) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6430407) q[1];
sx q[1];
rz(-1.828755) q[1];
sx q[1];
rz(3.113913) q[1];
x q[2];
rz(1.4164044) q[3];
sx q[3];
rz(-0.61617154) q[3];
sx q[3];
rz(2.7026724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7777286) q[2];
sx q[2];
rz(-1.2960459) q[2];
sx q[2];
rz(-1.9484005) q[2];
rz(-0.39919546) q[3];
sx q[3];
rz(-1.7292855) q[3];
sx q[3];
rz(0.027755888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6804009) q[0];
sx q[0];
rz(-1.8541279) q[0];
sx q[0];
rz(0.68049085) q[0];
rz(2.3091799) q[1];
sx q[1];
rz(-1.2544371) q[1];
sx q[1];
rz(2.9852273) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.098861) q[0];
sx q[0];
rz(-1.2953238) q[0];
sx q[0];
rz(-1.6252489) q[0];
x q[1];
rz(-0.83145468) q[2];
sx q[2];
rz(-1.0359087) q[2];
sx q[2];
rz(1.2818467) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5030843) q[1];
sx q[1];
rz(-1.8482395) q[1];
sx q[1];
rz(-1.2289407) q[1];
rz(-pi) q[2];
x q[2];
rz(1.353305) q[3];
sx q[3];
rz(-1.6902349) q[3];
sx q[3];
rz(2.7986106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.51123315) q[2];
sx q[2];
rz(-0.66143051) q[2];
sx q[2];
rz(-0.88999256) q[2];
rz(0.47641274) q[3];
sx q[3];
rz(-2.2178631) q[3];
sx q[3];
rz(-2.711003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7156859) q[0];
sx q[0];
rz(-1.8222734) q[0];
sx q[0];
rz(-2.529378) q[0];
rz(-1.7828434) q[1];
sx q[1];
rz(-2.3814059) q[1];
sx q[1];
rz(0.30119687) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33443794) q[0];
sx q[0];
rz(-2.8557192) q[0];
sx q[0];
rz(-1.8854333) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1015755) q[2];
sx q[2];
rz(-2.42832) q[2];
sx q[2];
rz(2.3281946) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.67415706) q[1];
sx q[1];
rz(-1.9922207) q[1];
sx q[1];
rz(1.9860877) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4047818) q[3];
sx q[3];
rz(-2.7932146) q[3];
sx q[3];
rz(-2.762037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.907454) q[2];
sx q[2];
rz(-0.80628866) q[2];
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
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.909914) q[0];
sx q[0];
rz(-0.34220085) q[0];
sx q[0];
rz(-2.8177596) q[0];
rz(-0.58770761) q[1];
sx q[1];
rz(-0.67981845) q[1];
sx q[1];
rz(-0.30276611) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97070951) q[0];
sx q[0];
rz(-1.3178749) q[0];
sx q[0];
rz(-0.33527261) q[0];
rz(-pi) q[1];
rz(-1.6195756) q[2];
sx q[2];
rz(-1.3338727) q[2];
sx q[2];
rz(-2.0120442) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0628793) q[1];
sx q[1];
rz(-2.4643366) q[1];
sx q[1];
rz(2.1724764) q[1];
x q[2];
rz(2.1015424) q[3];
sx q[3];
rz(-2.663739) q[3];
sx q[3];
rz(-2.8213906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6451463) q[2];
sx q[2];
rz(-2.1108997) q[2];
sx q[2];
rz(-0.56078625) q[2];
rz(-2.136611) q[3];
sx q[3];
rz(-1.8533665) q[3];
sx q[3];
rz(2.0463478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0305369) q[0];
sx q[0];
rz(-0.66937864) q[0];
sx q[0];
rz(-0.74991599) q[0];
rz(-2.7610682) q[1];
sx q[1];
rz(-0.98697248) q[1];
sx q[1];
rz(0.97533018) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98816865) q[0];
sx q[0];
rz(-1.0215217) q[0];
sx q[0];
rz(-2.2934521) q[0];
x q[1];
rz(1.777095) q[2];
sx q[2];
rz(-0.28536428) q[2];
sx q[2];
rz(2.9138164) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.9561466) q[1];
sx q[1];
rz(-0.6022033) q[1];
sx q[1];
rz(-0.57149296) q[1];
rz(-pi) q[2];
rz(-0.22244979) q[3];
sx q[3];
rz(-1.8651903) q[3];
sx q[3];
rz(2.6840212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9425977) q[2];
sx q[2];
rz(-0.90513217) q[2];
sx q[2];
rz(1.0814166) q[2];
rz(0.069247581) q[3];
sx q[3];
rz(-1.4360177) q[3];
sx q[3];
rz(0.92250219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13114318) q[0];
sx q[0];
rz(-0.32587019) q[0];
sx q[0];
rz(0.3057873) q[0];
rz(-0.30793134) q[1];
sx q[1];
rz(-1.6653857) q[1];
sx q[1];
rz(1.1526795) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15241218) q[0];
sx q[0];
rz(-2.3539641) q[0];
sx q[0];
rz(-0.37352011) q[0];
x q[1];
rz(-0.20736097) q[2];
sx q[2];
rz(-2.2054858) q[2];
sx q[2];
rz(0.35516741) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4977485) q[1];
sx q[1];
rz(-1.1670928) q[1];
sx q[1];
rz(2.7875442) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1265101) q[3];
sx q[3];
rz(-1.6341356) q[3];
sx q[3];
rz(2.7787105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5311188) q[2];
sx q[2];
rz(-1.0627154) q[2];
sx q[2];
rz(1.6884241) q[2];
rz(-2.6609663) q[3];
sx q[3];
rz(-2.5745945) q[3];
sx q[3];
rz(-1.3948729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4881445) q[0];
sx q[0];
rz(-1.5339889) q[0];
sx q[0];
rz(1.4657159) q[0];
rz(2.0987971) q[1];
sx q[1];
rz(-1.7386309) q[1];
sx q[1];
rz(2.244619) q[1];
rz(1.8389134) q[2];
sx q[2];
rz(-2.0901879) q[2];
sx q[2];
rz(-1.2218634) q[2];
rz(0.10877175) q[3];
sx q[3];
rz(-0.36009195) q[3];
sx q[3];
rz(2.910955) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
