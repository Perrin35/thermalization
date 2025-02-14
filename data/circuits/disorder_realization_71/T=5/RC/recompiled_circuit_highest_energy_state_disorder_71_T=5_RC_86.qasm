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
rz(1.1310391) q[0];
sx q[0];
rz(-0.57589632) q[0];
sx q[0];
rz(1.1269215) q[0];
rz(2.1908886) q[1];
sx q[1];
rz(-2.5403412) q[1];
sx q[1];
rz(0.33369219) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5883049) q[0];
sx q[0];
rz(-1.9899564) q[0];
sx q[0];
rz(0.16417154) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26716455) q[2];
sx q[2];
rz(-1.1295302) q[2];
sx q[2];
rz(1.6187606) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0036949) q[1];
sx q[1];
rz(-0.40781355) q[1];
sx q[1];
rz(-0.56038709) q[1];
rz(-pi) q[2];
rz(0.80218704) q[3];
sx q[3];
rz(-1.7628551) q[3];
sx q[3];
rz(0.036969846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2600962) q[2];
sx q[2];
rz(-1.0755971) q[2];
sx q[2];
rz(-2.0913701) q[2];
rz(0.29551926) q[3];
sx q[3];
rz(-2.3727356) q[3];
sx q[3];
rz(-1.7723005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9722026) q[0];
sx q[0];
rz(-1.0193595) q[0];
sx q[0];
rz(0.66725677) q[0];
rz(-0.15070209) q[1];
sx q[1];
rz(-2.0867911) q[1];
sx q[1];
rz(-0.31164718) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2209476) q[0];
sx q[0];
rz(-1.049982) q[0];
sx q[0];
rz(-0.08128165) q[0];
x q[1];
rz(-1.1318593) q[2];
sx q[2];
rz(-1.8510185) q[2];
sx q[2];
rz(-0.2187905) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.91193059) q[1];
sx q[1];
rz(-1.4313102) q[1];
sx q[1];
rz(-2.9565977) q[1];
x q[2];
rz(1.6436348) q[3];
sx q[3];
rz(-1.1741116) q[3];
sx q[3];
rz(-1.9997627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47505891) q[2];
sx q[2];
rz(-0.89343137) q[2];
sx q[2];
rz(-2.3739572) q[2];
rz(0.4808937) q[3];
sx q[3];
rz(-0.77156639) q[3];
sx q[3];
rz(-1.5546999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84592205) q[0];
sx q[0];
rz(-1.5376115) q[0];
sx q[0];
rz(2.3620102) q[0];
rz(0.37173158) q[1];
sx q[1];
rz(-1.0075684) q[1];
sx q[1];
rz(-0.76098162) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62255732) q[0];
sx q[0];
rz(-1.8599959) q[0];
sx q[0];
rz(-1.9054806) q[0];
x q[1];
rz(0.8849319) q[2];
sx q[2];
rz(-1.2239309) q[2];
sx q[2];
rz(0.82277966) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5058712) q[1];
sx q[1];
rz(-0.85800115) q[1];
sx q[1];
rz(0.21943032) q[1];
rz(-0.85270057) q[3];
sx q[3];
rz(-1.5745244) q[3];
sx q[3];
rz(-1.119233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3373105) q[2];
sx q[2];
rz(-0.90039841) q[2];
sx q[2];
rz(-2.6864181) q[2];
rz(-0.69194397) q[3];
sx q[3];
rz(-2.4271836) q[3];
sx q[3];
rz(-0.80879319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96875018) q[0];
sx q[0];
rz(-1.3469561) q[0];
sx q[0];
rz(-2.8853048) q[0];
rz(-1.7068656) q[1];
sx q[1];
rz(-1.7211569) q[1];
sx q[1];
rz(-0.11875471) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1349604) q[0];
sx q[0];
rz(-1.7259571) q[0];
sx q[0];
rz(0.4520024) q[0];
rz(0.1280814) q[2];
sx q[2];
rz(-1.8616779) q[2];
sx q[2];
rz(2.0258935) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9164711) q[1];
sx q[1];
rz(-2.3206704) q[1];
sx q[1];
rz(-0.31912843) q[1];
x q[2];
rz(1.9471274) q[3];
sx q[3];
rz(-2.2897165) q[3];
sx q[3];
rz(-0.014499078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.28748301) q[2];
sx q[2];
rz(-2.867925) q[2];
sx q[2];
rz(2.0859094) q[2];
rz(1.6915197) q[3];
sx q[3];
rz(-1.4472716) q[3];
sx q[3];
rz(2.292574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5797822) q[0];
sx q[0];
rz(-2.9636443) q[0];
sx q[0];
rz(1.5280888) q[0];
rz(1.9644507) q[1];
sx q[1];
rz(-1.8714995) q[1];
sx q[1];
rz(1.5632163) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10229853) q[0];
sx q[0];
rz(-1.5015409) q[0];
sx q[0];
rz(-0.62634344) q[0];
x q[1];
rz(-1.8991532) q[2];
sx q[2];
rz(-1.17982) q[2];
sx q[2];
rz(-0.41074387) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0381752) q[1];
sx q[1];
rz(-2.4929873) q[1];
sx q[1];
rz(-0.23596779) q[1];
rz(-0.40940955) q[3];
sx q[3];
rz(-1.5732523) q[3];
sx q[3];
rz(-2.5428605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13581181) q[2];
sx q[2];
rz(-1.9848738) q[2];
sx q[2];
rz(-1.8394252) q[2];
rz(0.33995134) q[3];
sx q[3];
rz(-2.3348742) q[3];
sx q[3];
rz(-0.66238856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0078916773) q[0];
sx q[0];
rz(-1.5444724) q[0];
sx q[0];
rz(1.9418465) q[0];
rz(-1.3613191) q[1];
sx q[1];
rz(-2.168455) q[1];
sx q[1];
rz(-2.5637085) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91980714) q[0];
sx q[0];
rz(-3.1263906) q[0];
sx q[0];
rz(-2.0073246) q[0];
x q[1];
rz(0.43283845) q[2];
sx q[2];
rz(-1.7146401) q[2];
sx q[2];
rz(1.270831) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.51398811) q[1];
sx q[1];
rz(-1.1179233) q[1];
sx q[1];
rz(1.4195561) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6541566) q[3];
sx q[3];
rz(-0.29988403) q[3];
sx q[3];
rz(2.5725627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.63152385) q[2];
sx q[2];
rz(-0.57454595) q[2];
sx q[2];
rz(0.4734545) q[2];
rz(-1.2724642) q[3];
sx q[3];
rz(-1.3926287) q[3];
sx q[3];
rz(-2.9191391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94721395) q[0];
sx q[0];
rz(-2.4838303) q[0];
sx q[0];
rz(0.33313242) q[0];
rz(0.22854742) q[1];
sx q[1];
rz(-1.8014149) q[1];
sx q[1];
rz(-2.5659335) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059371297) q[0];
sx q[0];
rz(-2.3239808) q[0];
sx q[0];
rz(-0.14540093) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2387462) q[2];
sx q[2];
rz(-1.1162469) q[2];
sx q[2];
rz(-0.12388307) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0927432) q[1];
sx q[1];
rz(-1.3710088) q[1];
sx q[1];
rz(1.9034991) q[1];
rz(-pi) q[2];
rz(-0.46374627) q[3];
sx q[3];
rz(-2.6131242) q[3];
sx q[3];
rz(2.1338045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8754742) q[2];
sx q[2];
rz(-0.55142752) q[2];
sx q[2];
rz(-1.6654588) q[2];
rz(2.0406145) q[3];
sx q[3];
rz(-2.0783547) q[3];
sx q[3];
rz(-2.6235918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5440893) q[0];
sx q[0];
rz(-2.2552555) q[0];
sx q[0];
rz(2.8344717) q[0];
rz(2.8111474) q[1];
sx q[1];
rz(-1.5437061) q[1];
sx q[1];
rz(1.7764567) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86548818) q[0];
sx q[0];
rz(-1.5886602) q[0];
sx q[0];
rz(-0.081327036) q[0];
rz(-pi) q[1];
rz(-2.6583985) q[2];
sx q[2];
rz(-2.7527134) q[2];
sx q[2];
rz(-2.4318397) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4933984) q[1];
sx q[1];
rz(-1.9983564) q[1];
sx q[1];
rz(-0.28458622) q[1];
x q[2];
rz(0.14162986) q[3];
sx q[3];
rz(-0.83811578) q[3];
sx q[3];
rz(-1.70873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8660628) q[2];
sx q[2];
rz(-2.5538462) q[2];
sx q[2];
rz(2.8928939) q[2];
rz(-1.2134877) q[3];
sx q[3];
rz(-1.4444193) q[3];
sx q[3];
rz(3.0799589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1520749) q[0];
sx q[0];
rz(-0.9062506) q[0];
sx q[0];
rz(0.686598) q[0];
rz(1.1129414) q[1];
sx q[1];
rz(-0.80250347) q[1];
sx q[1];
rz(-2.8259605) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7677143) q[0];
sx q[0];
rz(-1.6168103) q[0];
sx q[0];
rz(-1.2476646) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3817673) q[2];
sx q[2];
rz(-1.7526327) q[2];
sx q[2];
rz(-2.1450617) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.823608) q[1];
sx q[1];
rz(-0.66901842) q[1];
sx q[1];
rz(1.4186335) q[1];
rz(-1.6417362) q[3];
sx q[3];
rz(-2.3070863) q[3];
sx q[3];
rz(-0.77674538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45372066) q[2];
sx q[2];
rz(-0.56679711) q[2];
sx q[2];
rz(-2.4570214) q[2];
rz(1.941393) q[3];
sx q[3];
rz(-2.0122416) q[3];
sx q[3];
rz(-2.9620192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0612563) q[0];
sx q[0];
rz(-0.48114023) q[0];
sx q[0];
rz(0.0048333724) q[0];
rz(1.6239411) q[1];
sx q[1];
rz(-2.3976517) q[1];
sx q[1];
rz(2.8878816) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64297241) q[0];
sx q[0];
rz(-1.360384) q[0];
sx q[0];
rz(-0.087281373) q[0];
x q[1];
rz(2.6275835) q[2];
sx q[2];
rz(-0.80482641) q[2];
sx q[2];
rz(-2.6857306) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4082807) q[1];
sx q[1];
rz(-1.3752573) q[1];
sx q[1];
rz(0.89082373) q[1];
rz(-pi) q[2];
rz(-1.3718448) q[3];
sx q[3];
rz(-2.1124438) q[3];
sx q[3];
rz(1.0197848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.97810513) q[2];
sx q[2];
rz(-1.2349962) q[2];
sx q[2];
rz(1.182349) q[2];
rz(-2.4649418) q[3];
sx q[3];
rz(-2.7550321) q[3];
sx q[3];
rz(-2.1277908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61426281) q[0];
sx q[0];
rz(-0.78582055) q[0];
sx q[0];
rz(0.2926122) q[0];
rz(2.09265) q[1];
sx q[1];
rz(-1.7221778) q[1];
sx q[1];
rz(-2.4377951) q[1];
rz(1.9090777) q[2];
sx q[2];
rz(-2.0192374) q[2];
sx q[2];
rz(1.1981525) q[2];
rz(1.2037892) q[3];
sx q[3];
rz(-1.4784151) q[3];
sx q[3];
rz(-2.9611931) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
