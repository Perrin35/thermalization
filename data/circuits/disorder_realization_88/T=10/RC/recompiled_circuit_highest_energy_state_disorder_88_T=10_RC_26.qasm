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
rz(2.9945381) q[0];
sx q[0];
rz(-1.7400063) q[0];
sx q[0];
rz(2.2214878) q[0];
rz(-0.080634557) q[1];
sx q[1];
rz(3.7205003) q[1];
sx q[1];
rz(10.401934) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9465222) q[0];
sx q[0];
rz(-1.3597288) q[0];
sx q[0];
rz(-0.37977438) q[0];
x q[1];
rz(-1.5571655) q[2];
sx q[2];
rz(-1.7030099) q[2];
sx q[2];
rz(1.4135828) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5861533) q[1];
sx q[1];
rz(-1.7057014) q[1];
sx q[1];
rz(-2.5389266) q[1];
x q[2];
rz(-0.36698384) q[3];
sx q[3];
rz(-2.4189848) q[3];
sx q[3];
rz(-2.5740576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5795634) q[2];
sx q[2];
rz(-1.9271489) q[2];
sx q[2];
rz(-0.62082949) q[2];
rz(0.51554716) q[3];
sx q[3];
rz(-1.4225057) q[3];
sx q[3];
rz(2.2990885) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8949378) q[0];
sx q[0];
rz(-1.5351013) q[0];
sx q[0];
rz(-2.5625693) q[0];
rz(-0.82194263) q[1];
sx q[1];
rz(-1.5162946) q[1];
sx q[1];
rz(0.45375219) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61868405) q[0];
sx q[0];
rz(-2.4755619) q[0];
sx q[0];
rz(-0.36938195) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46644174) q[2];
sx q[2];
rz(-1.6714408) q[2];
sx q[2];
rz(-0.60324861) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.0099185506) q[1];
sx q[1];
rz(-1.6955396) q[1];
sx q[1];
rz(2.591955) q[1];
x q[2];
rz(2.2367495) q[3];
sx q[3];
rz(-0.65757759) q[3];
sx q[3];
rz(-0.50298728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4072676) q[2];
sx q[2];
rz(-0.13886034) q[2];
sx q[2];
rz(2.5767051) q[2];
rz(0.35483739) q[3];
sx q[3];
rz(-0.9762888) q[3];
sx q[3];
rz(1.423665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1693717) q[0];
sx q[0];
rz(-2.9300949) q[0];
sx q[0];
rz(-0.55150223) q[0];
rz(-0.36477271) q[1];
sx q[1];
rz(-2.8021937) q[1];
sx q[1];
rz(2.636748) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4992139) q[0];
sx q[0];
rz(-1.9758245) q[0];
sx q[0];
rz(-1.7444872) q[0];
rz(1.6644434) q[2];
sx q[2];
rz(-1.9818857) q[2];
sx q[2];
rz(-1.1716441) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3403459) q[1];
sx q[1];
rz(-0.97724229) q[1];
sx q[1];
rz(0.98538633) q[1];
rz(-pi) q[2];
rz(1.9471719) q[3];
sx q[3];
rz(-2.0224114) q[3];
sx q[3];
rz(-0.87290369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7079033) q[2];
sx q[2];
rz(-1.7283551) q[2];
sx q[2];
rz(-3.0943387) q[2];
rz(0.83682483) q[3];
sx q[3];
rz(-2.4182726) q[3];
sx q[3];
rz(2.1164472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4243917) q[0];
sx q[0];
rz(-2.1121139) q[0];
sx q[0];
rz(-1.8151872) q[0];
rz(1.6560417) q[1];
sx q[1];
rz(-2.0573719) q[1];
sx q[1];
rz(-0.63492376) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0795006) q[0];
sx q[0];
rz(-1.5720815) q[0];
sx q[0];
rz(1.5888402) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1702386) q[2];
sx q[2];
rz(-0.58341372) q[2];
sx q[2];
rz(0.82380481) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8532456) q[1];
sx q[1];
rz(-1.6430055) q[1];
sx q[1];
rz(-1.9068524) q[1];
x q[2];
rz(2.5136389) q[3];
sx q[3];
rz(-1.1276922) q[3];
sx q[3];
rz(-1.4301883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9012458) q[2];
sx q[2];
rz(-2.2525658) q[2];
sx q[2];
rz(1.3725613) q[2];
rz(-1.2076123) q[3];
sx q[3];
rz(-2.4977081) q[3];
sx q[3];
rz(-0.27455583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8127301) q[0];
sx q[0];
rz(-2.6565318) q[0];
sx q[0];
rz(3.0364756) q[0];
rz(-0.33991995) q[1];
sx q[1];
rz(-0.9143908) q[1];
sx q[1];
rz(-2.0786659) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30468291) q[0];
sx q[0];
rz(-1.5848918) q[0];
sx q[0];
rz(2.430116) q[0];
rz(-pi) q[1];
rz(1.3219484) q[2];
sx q[2];
rz(-2.3169059) q[2];
sx q[2];
rz(0.4363554) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2984933) q[1];
sx q[1];
rz(-1.5122422) q[1];
sx q[1];
rz(-1.5565518) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5751198) q[3];
sx q[3];
rz(-2.6427442) q[3];
sx q[3];
rz(-2.3259142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.062332705) q[2];
sx q[2];
rz(-1.3881114) q[2];
sx q[2];
rz(-1.8049392) q[2];
rz(-0.054232728) q[3];
sx q[3];
rz(-1.1905866) q[3];
sx q[3];
rz(2.5469053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17276758) q[0];
sx q[0];
rz(-2.4496138) q[0];
sx q[0];
rz(1.927595) q[0];
rz(1.0460151) q[1];
sx q[1];
rz(-2.0966625) q[1];
sx q[1];
rz(0.85375839) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19148286) q[0];
sx q[0];
rz(-1.6304509) q[0];
sx q[0];
rz(-1.7050939) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1519203) q[2];
sx q[2];
rz(-0.90384353) q[2];
sx q[2];
rz(-3.0158531) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5750372) q[1];
sx q[1];
rz(-2.2872529) q[1];
sx q[1];
rz(2.3385919) q[1];
rz(1.4183956) q[3];
sx q[3];
rz(-1.9160685) q[3];
sx q[3];
rz(-0.70122805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0115016) q[2];
sx q[2];
rz(-1.4902196) q[2];
sx q[2];
rz(-2.3940274) q[2];
rz(2.2085564) q[3];
sx q[3];
rz(-1.7739762) q[3];
sx q[3];
rz(-2.8396377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1099243) q[0];
sx q[0];
rz(-0.95369354) q[0];
sx q[0];
rz(-0.64144301) q[0];
rz(-2.172442) q[1];
sx q[1];
rz(-2.0717924) q[1];
sx q[1];
rz(-1.1136805) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9819558) q[0];
sx q[0];
rz(-0.41712077) q[0];
sx q[0];
rz(-1.3288767) q[0];
rz(0.93514438) q[2];
sx q[2];
rz(-2.0144267) q[2];
sx q[2];
rz(2.005524) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.42150527) q[1];
sx q[1];
rz(-2.0893789) q[1];
sx q[1];
rz(2.9353013) q[1];
rz(-pi) q[2];
rz(-1.5364755) q[3];
sx q[3];
rz(-1.0507108) q[3];
sx q[3];
rz(-1.8994562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.574719) q[2];
sx q[2];
rz(-0.69872624) q[2];
sx q[2];
rz(0.75801545) q[2];
rz(-2.8356683) q[3];
sx q[3];
rz(-0.35861349) q[3];
sx q[3];
rz(-2.6942286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7253983) q[0];
sx q[0];
rz(-0.60878009) q[0];
sx q[0];
rz(2.388227) q[0];
rz(-0.92195177) q[1];
sx q[1];
rz(-1.4267068) q[1];
sx q[1];
rz(0.13455483) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050722402) q[0];
sx q[0];
rz(-0.36733741) q[0];
sx q[0];
rz(-1.3274756) q[0];
x q[1];
rz(1.8250663) q[2];
sx q[2];
rz(-0.94506028) q[2];
sx q[2];
rz(-0.78142399) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0363732) q[1];
sx q[1];
rz(-1.660562) q[1];
sx q[1];
rz(-1.5898934) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5882072) q[3];
sx q[3];
rz(-0.7532922) q[3];
sx q[3];
rz(2.2224421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77384633) q[2];
sx q[2];
rz(-1.6951122) q[2];
sx q[2];
rz(-1.9473677) q[2];
rz(-1.9959244) q[3];
sx q[3];
rz(-1.1402036) q[3];
sx q[3];
rz(1.0660508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(1.0399748) q[0];
sx q[0];
rz(-0.45926738) q[0];
sx q[0];
rz(-0.38247821) q[0];
rz(-0.76599145) q[1];
sx q[1];
rz(-1.4975558) q[1];
sx q[1];
rz(1.3409748) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1857698) q[0];
sx q[0];
rz(-0.18322028) q[0];
sx q[0];
rz(-0.51143666) q[0];
rz(-pi) q[1];
rz(-0.75720301) q[2];
sx q[2];
rz(-1.6348038) q[2];
sx q[2];
rz(0.369095) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7241269) q[1];
sx q[1];
rz(-1.7977906) q[1];
sx q[1];
rz(-0.75717749) q[1];
rz(0.082499113) q[3];
sx q[3];
rz(-1.2159705) q[3];
sx q[3];
rz(0.88731836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0880903) q[2];
sx q[2];
rz(-0.92038766) q[2];
sx q[2];
rz(-0.11605334) q[2];
rz(-1.2608438) q[3];
sx q[3];
rz(-2.0033658) q[3];
sx q[3];
rz(-1.457823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55105227) q[0];
sx q[0];
rz(-0.11369471) q[0];
sx q[0];
rz(-0.1846479) q[0];
rz(0.91839904) q[1];
sx q[1];
rz(-1.5946439) q[1];
sx q[1];
rz(-1.437423) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1383863) q[0];
sx q[0];
rz(-2.0420364) q[0];
sx q[0];
rz(0.22881656) q[0];
rz(-pi) q[1];
rz(-0.080731656) q[2];
sx q[2];
rz(-0.7237607) q[2];
sx q[2];
rz(2.2249976) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9574163) q[1];
sx q[1];
rz(-2.1706958) q[1];
sx q[1];
rz(-2.9084212) q[1];
rz(-1.2380794) q[3];
sx q[3];
rz(-0.92192382) q[3];
sx q[3];
rz(1.0024662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.664428) q[2];
sx q[2];
rz(-1.8958586) q[2];
sx q[2];
rz(-2.5028382) q[2];
rz(-0.81513682) q[3];
sx q[3];
rz(-0.4709979) q[3];
sx q[3];
rz(-0.42069978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7246134) q[0];
sx q[0];
rz(-1.2495578) q[0];
sx q[0];
rz(0.15171262) q[0];
rz(-2.3607415) q[1];
sx q[1];
rz(-1.6897222) q[1];
sx q[1];
rz(2.2471468) q[1];
rz(3.0197418) q[2];
sx q[2];
rz(-1.8486628) q[2];
sx q[2];
rz(1.8589414) q[2];
rz(-2.9853447) q[3];
sx q[3];
rz(-2.8946946) q[3];
sx q[3];
rz(1.3794086) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
