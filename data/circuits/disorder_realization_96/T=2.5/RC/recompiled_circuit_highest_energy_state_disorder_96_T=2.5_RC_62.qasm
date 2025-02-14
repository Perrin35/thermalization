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
rz(1.9198298) q[0];
sx q[0];
rz(-1.3718995) q[0];
sx q[0];
rz(1.0839809) q[0];
rz(-0.78217512) q[1];
sx q[1];
rz(3.5909619) q[1];
sx q[1];
rz(10.211791) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5575268) q[0];
sx q[0];
rz(-1.1827994) q[0];
sx q[0];
rz(2.9783413) q[0];
rz(-pi) q[1];
rz(0.89504234) q[2];
sx q[2];
rz(-1.5359725) q[2];
sx q[2];
rz(0.36961242) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2725153) q[1];
sx q[1];
rz(-2.1109588) q[1];
sx q[1];
rz(-2.9754728) q[1];
rz(-0.34073273) q[3];
sx q[3];
rz(-0.91487802) q[3];
sx q[3];
rz(-2.4994196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3419753) q[2];
sx q[2];
rz(-1.3222597) q[2];
sx q[2];
rz(2.4413617) q[2];
rz(-1.7796984) q[3];
sx q[3];
rz(-1.9340065) q[3];
sx q[3];
rz(-3.0119058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48548651) q[0];
sx q[0];
rz(-0.34657297) q[0];
sx q[0];
rz(-1.19278) q[0];
rz(-2.9876409) q[1];
sx q[1];
rz(-0.62869453) q[1];
sx q[1];
rz(1.2492294) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19156583) q[0];
sx q[0];
rz(-2.9572672) q[0];
sx q[0];
rz(-2.7937355) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2080113) q[2];
sx q[2];
rz(-1.5270479) q[2];
sx q[2];
rz(1.9160567) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.44012516) q[1];
sx q[1];
rz(-2.8258913) q[1];
sx q[1];
rz(2.8786826) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7720501) q[3];
sx q[3];
rz(-0.93825723) q[3];
sx q[3];
rz(0.12128017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.93231702) q[2];
sx q[2];
rz(-1.2195769) q[2];
sx q[2];
rz(2.2323214) q[2];
rz(3.0026109) q[3];
sx q[3];
rz(-0.2551955) q[3];
sx q[3];
rz(-0.87066686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22055498) q[0];
sx q[0];
rz(-1.9746566) q[0];
sx q[0];
rz(-1.1880818) q[0];
rz(-2.2782169) q[1];
sx q[1];
rz(-1.8977576) q[1];
sx q[1];
rz(2.1720502) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6351417) q[0];
sx q[0];
rz(-1.7291082) q[0];
sx q[0];
rz(-2.979605) q[0];
x q[1];
rz(-2.7935384) q[2];
sx q[2];
rz(-1.9241986) q[2];
sx q[2];
rz(-0.60863335) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.66268209) q[1];
sx q[1];
rz(-1.3363774) q[1];
sx q[1];
rz(-1.471651) q[1];
rz(-0.86609716) q[3];
sx q[3];
rz(-1.1148165) q[3];
sx q[3];
rz(-0.34581071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7748176) q[2];
sx q[2];
rz(-2.884951) q[2];
sx q[2];
rz(2.4859264) q[2];
rz(-1.6171148) q[3];
sx q[3];
rz(-1.3527801) q[3];
sx q[3];
rz(-2.2881983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26576385) q[0];
sx q[0];
rz(-2.0142856) q[0];
sx q[0];
rz(-1.6831336) q[0];
rz(-1.4153076) q[1];
sx q[1];
rz(-1.4348607) q[1];
sx q[1];
rz(1.8728135) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7002633) q[0];
sx q[0];
rz(-2.1887795) q[0];
sx q[0];
rz(2.000314) q[0];
rz(-0.54019467) q[2];
sx q[2];
rz(-1.0120076) q[2];
sx q[2];
rz(-1.1891728) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3038588) q[1];
sx q[1];
rz(-1.7540252) q[1];
sx q[1];
rz(0.86838037) q[1];
x q[2];
rz(0.089667356) q[3];
sx q[3];
rz(-1.7981861) q[3];
sx q[3];
rz(2.1289153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.65712523) q[2];
sx q[2];
rz(-1.7871658) q[2];
sx q[2];
rz(1.8787059) q[2];
rz(0.2028939) q[3];
sx q[3];
rz(-1.032136) q[3];
sx q[3];
rz(1.3304905) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6577067) q[0];
sx q[0];
rz(-2.0575476) q[0];
sx q[0];
rz(-0.44802353) q[0];
rz(1.9810642) q[1];
sx q[1];
rz(-2.0612165) q[1];
sx q[1];
rz(2.0652658) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58680326) q[0];
sx q[0];
rz(-2.3856253) q[0];
sx q[0];
rz(-1.2584096) q[0];
rz(2.8084023) q[2];
sx q[2];
rz(-1.6181862) q[2];
sx q[2];
rz(-1.9006001) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.35713112) q[1];
sx q[1];
rz(-1.6069043) q[1];
sx q[1];
rz(1.5090843) q[1];
rz(-3.1212937) q[3];
sx q[3];
rz(-0.95385984) q[3];
sx q[3];
rz(-2.9391409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5356323) q[2];
sx q[2];
rz(-1.863402) q[2];
sx q[2];
rz(-0.34206259) q[2];
rz(-1.131912) q[3];
sx q[3];
rz(-1.071238) q[3];
sx q[3];
rz(0.38558495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1256506) q[0];
sx q[0];
rz(-0.60369879) q[0];
sx q[0];
rz(2.0397256) q[0];
rz(2.8562538) q[1];
sx q[1];
rz(-0.97831786) q[1];
sx q[1];
rz(-2.2314821) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21134899) q[0];
sx q[0];
rz(-1.4574629) q[0];
sx q[0];
rz(1.428753) q[0];
rz(-pi) q[1];
rz(1.455972) q[2];
sx q[2];
rz(-1.6636724) q[2];
sx q[2];
rz(0.50485134) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.460798) q[1];
sx q[1];
rz(-2.6943047) q[1];
sx q[1];
rz(-1.7030925) q[1];
rz(0.85569622) q[3];
sx q[3];
rz(-0.2350546) q[3];
sx q[3];
rz(2.1657221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3204699) q[2];
sx q[2];
rz(-1.4218825) q[2];
sx q[2];
rz(-1.0538496) q[2];
rz(-2.7779135) q[3];
sx q[3];
rz(-1.4854919) q[3];
sx q[3];
rz(0.53500879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4800471) q[0];
sx q[0];
rz(-2.1286025) q[0];
sx q[0];
rz(0.15740982) q[0];
rz(2.5120381) q[1];
sx q[1];
rz(-1.5512356) q[1];
sx q[1];
rz(-0.079364337) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4422121) q[0];
sx q[0];
rz(-1.6621502) q[0];
sx q[0];
rz(-3.0708036) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39418295) q[2];
sx q[2];
rz(-1.7346343) q[2];
sx q[2];
rz(2.8904786) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.094202894) q[1];
sx q[1];
rz(-1.6999287) q[1];
sx q[1];
rz(2.6455621) q[1];
rz(-pi) q[2];
rz(1.6116051) q[3];
sx q[3];
rz(-2.5548078) q[3];
sx q[3];
rz(-2.2631078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.51084149) q[2];
sx q[2];
rz(-1.3151104) q[2];
sx q[2];
rz(0.17244478) q[2];
rz(0.6423966) q[3];
sx q[3];
rz(-2.1746217) q[3];
sx q[3];
rz(2.7914458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.073535) q[0];
sx q[0];
rz(-1.302916) q[0];
sx q[0];
rz(0.56655836) q[0];
rz(0.22816518) q[1];
sx q[1];
rz(-1.6888432) q[1];
sx q[1];
rz(-1.8295005) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.909186) q[0];
sx q[0];
rz(-2.3433422) q[0];
sx q[0];
rz(1.2212103) q[0];
rz(0.117339) q[2];
sx q[2];
rz(-2.2944231) q[2];
sx q[2];
rz(3.0229508) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.65441865) q[1];
sx q[1];
rz(-1.0680334) q[1];
sx q[1];
rz(-2.1969686) q[1];
x q[2];
rz(-1.4890461) q[3];
sx q[3];
rz(-2.7513732) q[3];
sx q[3];
rz(-1.8070328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9606026) q[2];
sx q[2];
rz(-1.484551) q[2];
sx q[2];
rz(0.20299882) q[2];
rz(2.3024043) q[3];
sx q[3];
rz(-2.8036717) q[3];
sx q[3];
rz(2.2904229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39076552) q[0];
sx q[0];
rz(-0.44742328) q[0];
sx q[0];
rz(2.9946193) q[0];
rz(-2.7087063) q[1];
sx q[1];
rz(-1.1914445) q[1];
sx q[1];
rz(1.8870707) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4449249) q[0];
sx q[0];
rz(-2.4913209) q[0];
sx q[0];
rz(2.5461322) q[0];
rz(-pi) q[1];
rz(-1.3225088) q[2];
sx q[2];
rz(-1.6035558) q[2];
sx q[2];
rz(1.133762) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7045198) q[1];
sx q[1];
rz(-0.38846782) q[1];
sx q[1];
rz(2.0411496) q[1];
x q[2];
rz(0.12740429) q[3];
sx q[3];
rz(-0.16949305) q[3];
sx q[3];
rz(-2.3338846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.36048421) q[2];
sx q[2];
rz(-2.1238748) q[2];
sx q[2];
rz(-2.600889) q[2];
rz(-2.8402719) q[3];
sx q[3];
rz(-2.7408528) q[3];
sx q[3];
rz(1.2979243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8526469) q[0];
sx q[0];
rz(-1.9444554) q[0];
sx q[0];
rz(0.72809118) q[0];
rz(0.51746619) q[1];
sx q[1];
rz(-1.4897852) q[1];
sx q[1];
rz(-0.99658406) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.655642) q[0];
sx q[0];
rz(-0.34531718) q[0];
sx q[0];
rz(1.6486069) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9884981) q[2];
sx q[2];
rz(-0.587112) q[2];
sx q[2];
rz(2.401641) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6428553) q[1];
sx q[1];
rz(-1.3142546) q[1];
sx q[1];
rz(-3.1270621) q[1];
rz(1.4703937) q[3];
sx q[3];
rz(-1.7531503) q[3];
sx q[3];
rz(-2.792458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1835798) q[2];
sx q[2];
rz(-1.5637584) q[2];
sx q[2];
rz(0.070240423) q[2];
rz(3.1349414) q[3];
sx q[3];
rz(-2.9700322) q[3];
sx q[3];
rz(0.10449617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9329279) q[0];
sx q[0];
rz(-1.4069955) q[0];
sx q[0];
rz(1.7250742) q[0];
rz(-0.76864645) q[1];
sx q[1];
rz(-1.9192764) q[1];
sx q[1];
rz(2.8142014) q[1];
rz(-1.6835536) q[2];
sx q[2];
rz(-1.590645) q[2];
sx q[2];
rz(-2.0525862) q[2];
rz(2.8715618) q[3];
sx q[3];
rz(-1.217801) q[3];
sx q[3];
rz(1.4767811) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
