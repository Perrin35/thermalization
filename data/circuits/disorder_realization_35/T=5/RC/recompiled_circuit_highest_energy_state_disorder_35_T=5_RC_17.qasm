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
rz(4.0539157) q[0];
sx q[0];
rz(11.231448) q[0];
rz(-1.1420684) q[1];
sx q[1];
rz(5.7647709) q[1];
sx q[1];
rz(7.9328598) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8700884) q[0];
sx q[0];
rz(-0.056500204) q[0];
sx q[0];
rz(-1.9046049) q[0];
rz(-pi) q[1];
rz(-0.9426192) q[2];
sx q[2];
rz(-2.16428) q[2];
sx q[2];
rz(-1.3243937) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.277093) q[1];
sx q[1];
rz(-0.96597176) q[1];
sx q[1];
rz(-0.34326174) q[1];
rz(-pi) q[2];
rz(-1.7068976) q[3];
sx q[3];
rz(-0.24472642) q[3];
sx q[3];
rz(-0.32258209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3492744) q[2];
sx q[2];
rz(-2.8742542) q[2];
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
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2199663) q[0];
sx q[0];
rz(-1.3082137) q[0];
sx q[0];
rz(-0.24800214) q[0];
rz(1.8964881) q[1];
sx q[1];
rz(-2.4047132) q[1];
sx q[1];
rz(-1.9128333) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4061444) q[0];
sx q[0];
rz(-1.7559253) q[0];
sx q[0];
rz(1.6044751) q[0];
rz(-2.1013772) q[2];
sx q[2];
rz(-1.1184012) q[2];
sx q[2];
rz(1.4380406) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5723443) q[1];
sx q[1];
rz(-0.73671417) q[1];
sx q[1];
rz(-0.31804292) q[1];
x q[2];
rz(3.008083) q[3];
sx q[3];
rz(-2.7179681) q[3];
sx q[3];
rz(-2.6299919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9241141) q[2];
sx q[2];
rz(-0.83196297) q[2];
sx q[2];
rz(2.5577616) q[2];
rz(2.7122688) q[3];
sx q[3];
rz(-1.2238294) q[3];
sx q[3];
rz(1.0113641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9871224) q[0];
sx q[0];
rz(-2.7543289) q[0];
sx q[0];
rz(-1.6744457) q[0];
rz(1.6636498) q[1];
sx q[1];
rz(-1.7324305) q[1];
sx q[1];
rz(-1.0999701) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80187782) q[0];
sx q[0];
rz(-0.22476443) q[0];
sx q[0];
rz(-2.6532214) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9103545) q[2];
sx q[2];
rz(-1.8446088) q[2];
sx q[2];
rz(1.6664315) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.28274714) q[1];
sx q[1];
rz(-1.2371776) q[1];
sx q[1];
rz(-0.40550123) q[1];
x q[2];
rz(-0.99697014) q[3];
sx q[3];
rz(-1.9553292) q[3];
sx q[3];
rz(-1.9652733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.786342) q[2];
sx q[2];
rz(-2.2189271) q[2];
sx q[2];
rz(1.5938909) q[2];
rz(1.330438) q[3];
sx q[3];
rz(-1.3733613) q[3];
sx q[3];
rz(-2.0681341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.380577) q[0];
sx q[0];
rz(-1.8208193) q[0];
sx q[0];
rz(-0.15783489) q[0];
rz(-0.24179587) q[1];
sx q[1];
rz(-2.7561185) q[1];
sx q[1];
rz(-2.6604624) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0711827) q[0];
sx q[0];
rz(-0.86127036) q[0];
sx q[0];
rz(1.6874403) q[0];
rz(-pi) q[1];
rz(-0.60637668) q[2];
sx q[2];
rz(-0.42841347) q[2];
sx q[2];
rz(2.1386752) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2896636) q[1];
sx q[1];
rz(-2.4074775) q[1];
sx q[1];
rz(1.0490225) q[1];
x q[2];
rz(2.1948482) q[3];
sx q[3];
rz(-2.7928154) q[3];
sx q[3];
rz(-0.52979505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2699997) q[2];
sx q[2];
rz(-2.3303878) q[2];
sx q[2];
rz(0.2743741) q[2];
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
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45193732) q[0];
sx q[0];
rz(-2.2640197) q[0];
sx q[0];
rz(2.080132) q[0];
rz(0.44867107) q[1];
sx q[1];
rz(-0.45495382) q[1];
sx q[1];
rz(1.7015069) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3502304) q[0];
sx q[0];
rz(-1.4141948) q[0];
sx q[0];
rz(-0.76275218) q[0];
rz(1.1665197) q[2];
sx q[2];
rz(-1.6875522) q[2];
sx q[2];
rz(1.5991885) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0764112) q[1];
sx q[1];
rz(-1.59756) q[1];
sx q[1];
rz(1.3127432) q[1];
rz(2.1813356) q[3];
sx q[3];
rz(-1.4818076) q[3];
sx q[3];
rz(-1.2582092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3638641) q[2];
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
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4611918) q[0];
sx q[0];
rz(-1.8541279) q[0];
sx q[0];
rz(2.4611018) q[0];
rz(0.83241278) q[1];
sx q[1];
rz(-1.2544371) q[1];
sx q[1];
rz(0.15636538) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6844821) q[0];
sx q[0];
rz(-1.5183987) q[0];
sx q[0];
rz(-2.8657317) q[0];
rz(-0.84951289) q[2];
sx q[2];
rz(-0.88187432) q[2];
sx q[2];
rz(-2.920814) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1121516) q[1];
sx q[1];
rz(-1.8990771) q[1];
sx q[1];
rz(-0.29354696) q[1];
rz(-2.0783608) q[3];
sx q[3];
rz(-0.2476736) q[3];
sx q[3];
rz(1.4192691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51123315) q[2];
sx q[2];
rz(-2.4801621) q[2];
sx q[2];
rz(0.88999256) q[2];
rz(-0.47641274) q[3];
sx q[3];
rz(-2.2178631) q[3];
sx q[3];
rz(-0.43058968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42590672) q[0];
sx q[0];
rz(-1.3193193) q[0];
sx q[0];
rz(-2.529378) q[0];
rz(-1.3587492) q[1];
sx q[1];
rz(-0.76018676) q[1];
sx q[1];
rz(-2.8403958) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33443794) q[0];
sx q[0];
rz(-0.2858735) q[0];
sx q[0];
rz(-1.2561594) q[0];
rz(-2.0400171) q[2];
sx q[2];
rz(-0.71327268) q[2];
sx q[2];
rz(2.3281946) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4674356) q[1];
sx q[1];
rz(-1.9922207) q[1];
sx q[1];
rz(-1.9860877) q[1];
x q[2];
rz(1.3314358) q[3];
sx q[3];
rz(-1.826397) q[3];
sx q[3];
rz(1.1472697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.907454) q[2];
sx q[2];
rz(-0.80628866) q[2];
sx q[2];
rz(2.0937199) q[2];
rz(0.35510865) q[3];
sx q[3];
rz(-0.57098782) q[3];
sx q[3];
rz(2.4560438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.909914) q[0];
sx q[0];
rz(-0.34220085) q[0];
sx q[0];
rz(0.32383305) q[0];
rz(-0.58770761) q[1];
sx q[1];
rz(-0.67981845) q[1];
sx q[1];
rz(-0.30276611) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9186818) q[0];
sx q[0];
rz(-0.41707539) q[0];
sx q[0];
rz(-0.66584754) q[0];
rz(2.9043973) q[2];
sx q[2];
rz(-1.5233808) q[2];
sx q[2];
rz(0.42979017) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.98339048) q[1];
sx q[1];
rz(-1.9333955) q[1];
sx q[1];
rz(0.98539488) q[1];
x q[2];
rz(-0.25637324) q[3];
sx q[3];
rz(-1.1629761) q[3];
sx q[3];
rz(2.2373509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4964464) q[2];
sx q[2];
rz(-2.1108997) q[2];
sx q[2];
rz(-0.56078625) q[2];
rz(-1.0049817) q[3];
sx q[3];
rz(-1.8533665) q[3];
sx q[3];
rz(1.0952449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
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
rz(0.97533018) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0935287) q[0];
sx q[0];
rz(-0.87665999) q[0];
sx q[0];
rz(-2.3175236) q[0];
x q[1];
rz(1.777095) q[2];
sx q[2];
rz(-0.28536428) q[2];
sx q[2];
rz(-0.22777624) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.9561466) q[1];
sx q[1];
rz(-2.5393894) q[1];
sx q[1];
rz(0.57149296) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8721902) q[3];
sx q[3];
rz(-1.3580672) q[3];
sx q[3];
rz(2.093906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19899496) q[2];
sx q[2];
rz(-2.2364605) q[2];
sx q[2];
rz(-1.0814166) q[2];
rz(3.0723451) q[3];
sx q[3];
rz(-1.705575) q[3];
sx q[3];
rz(0.92250219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13114318) q[0];
sx q[0];
rz(-2.8157225) q[0];
sx q[0];
rz(-2.8358054) q[0];
rz(-0.30793134) q[1];
sx q[1];
rz(-1.4762069) q[1];
sx q[1];
rz(1.9889132) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.688153) q[0];
sx q[0];
rz(-1.3092293) q[0];
sx q[0];
rz(-0.75193172) q[0];
x q[1];
rz(1.2981554) q[2];
sx q[2];
rz(-0.66323384) q[2];
sx q[2];
rz(-3.1273914) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.252723) q[1];
sx q[1];
rz(-2.6111341) q[1];
sx q[1];
rz(2.2525851) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6341428) q[3];
sx q[3];
rz(-1.5858486) q[3];
sx q[3];
rz(-1.2069595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61047381) q[2];
sx q[2];
rz(-1.0627154) q[2];
sx q[2];
rz(1.4531685) q[2];
rz(-0.48062634) q[3];
sx q[3];
rz(-0.56699816) q[3];
sx q[3];
rz(1.7467197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65344812) q[0];
sx q[0];
rz(-1.5339889) q[0];
sx q[0];
rz(1.4657159) q[0];
rz(-1.0427955) q[1];
sx q[1];
rz(-1.7386309) q[1];
sx q[1];
rz(2.244619) q[1];
rz(-2.6063812) q[2];
sx q[2];
rz(-1.3387398) q[2];
sx q[2];
rz(-2.6571318) q[2];
rz(-1.5299464) q[3];
sx q[3];
rz(-1.9286641) q[3];
sx q[3];
rz(2.7947938) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
