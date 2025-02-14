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
rz(-0.1470546) q[0];
sx q[0];
rz(-1.4015863) q[0];
sx q[0];
rz(0.92010486) q[0];
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
rz(1.9465222) q[0];
sx q[0];
rz(-1.3597288) q[0];
sx q[0];
rz(0.37977438) q[0];
rz(1.5844272) q[2];
sx q[2];
rz(-1.4385828) q[2];
sx q[2];
rz(-1.4135828) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0339573) q[1];
sx q[1];
rz(-0.9743685) q[1];
sx q[1];
rz(-1.4075085) q[1];
rz(-pi) q[2];
rz(0.36698384) q[3];
sx q[3];
rz(-0.72260783) q[3];
sx q[3];
rz(-2.5740576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5795634) q[2];
sx q[2];
rz(-1.9271489) q[2];
sx q[2];
rz(-0.62082949) q[2];
rz(-2.6260455) q[3];
sx q[3];
rz(-1.4225057) q[3];
sx q[3];
rz(2.2990885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8949378) q[0];
sx q[0];
rz(-1.5351013) q[0];
sx q[0];
rz(-0.57902336) q[0];
rz(0.82194263) q[1];
sx q[1];
rz(-1.5162946) q[1];
sx q[1];
rz(2.6878405) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5229086) q[0];
sx q[0];
rz(-0.6660307) q[0];
sx q[0];
rz(-2.7722107) q[0];
x q[1];
rz(-2.6751509) q[2];
sx q[2];
rz(-1.6714408) q[2];
sx q[2];
rz(0.60324861) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1316741) q[1];
sx q[1];
rz(-1.6955396) q[1];
sx q[1];
rz(2.591955) q[1];
rz(-pi) q[2];
rz(-2.1165127) q[3];
sx q[3];
rz(-1.1835872) q[3];
sx q[3];
rz(-0.51147616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4072676) q[2];
sx q[2];
rz(-0.13886034) q[2];
sx q[2];
rz(-2.5767051) q[2];
rz(-0.35483739) q[3];
sx q[3];
rz(-0.9762888) q[3];
sx q[3];
rz(1.7179276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9722209) q[0];
sx q[0];
rz(-2.9300949) q[0];
sx q[0];
rz(-2.5900904) q[0];
rz(-0.36477271) q[1];
sx q[1];
rz(-0.33939895) q[1];
sx q[1];
rz(-2.636748) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0612978) q[0];
sx q[0];
rz(-2.7027931) q[0];
sx q[0];
rz(0.3831692) q[0];
rz(-pi) q[1];
rz(2.7288923) q[2];
sx q[2];
rz(-1.6566212) q[2];
sx q[2];
rz(-2.7049261) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.4144812) q[1];
sx q[1];
rz(-1.0950146) q[1];
sx q[1];
rz(-0.68051312) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1944207) q[3];
sx q[3];
rz(-1.1191812) q[3];
sx q[3];
rz(-2.268689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43368936) q[2];
sx q[2];
rz(-1.4132376) q[2];
sx q[2];
rz(-0.047253963) q[2];
rz(2.3047678) q[3];
sx q[3];
rz(-0.72332007) q[3];
sx q[3];
rz(2.1164472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.717201) q[0];
sx q[0];
rz(-2.1121139) q[0];
sx q[0];
rz(1.8151872) q[0];
rz(1.4855509) q[1];
sx q[1];
rz(-1.0842208) q[1];
sx q[1];
rz(2.5066689) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0795006) q[0];
sx q[0];
rz(-1.5720815) q[0];
sx q[0];
rz(-1.5888402) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35648326) q[2];
sx q[2];
rz(-1.0986137) q[2];
sx q[2];
rz(-3.0038578) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.078770854) q[1];
sx q[1];
rz(-0.34343849) q[1];
sx q[1];
rz(-1.354864) q[1];
rz(-pi) q[2];
rz(-0.62795378) q[3];
sx q[3];
rz(-2.0139004) q[3];
sx q[3];
rz(-1.7114044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2403468) q[2];
sx q[2];
rz(-0.88902688) q[2];
sx q[2];
rz(-1.3725613) q[2];
rz(1.2076123) q[3];
sx q[3];
rz(-2.4977081) q[3];
sx q[3];
rz(0.27455583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32886252) q[0];
sx q[0];
rz(-0.48506081) q[0];
sx q[0];
rz(-0.10511705) q[0];
rz(-0.33991995) q[1];
sx q[1];
rz(-2.2272019) q[1];
sx q[1];
rz(-1.0629268) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.88763) q[0];
sx q[0];
rz(-0.85940532) q[0];
sx q[0];
rz(1.5894058) q[0];
rz(-pi) q[1];
rz(0.26039298) q[2];
sx q[2];
rz(-2.3627334) q[2];
sx q[2];
rz(0.078140251) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8430994) q[1];
sx q[1];
rz(-1.6293504) q[1];
sx q[1];
rz(-1.5565518) q[1];
rz(-3.1392372) q[3];
sx q[3];
rz(-1.071953) q[3];
sx q[3];
rz(-2.3209907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0792599) q[2];
sx q[2];
rz(-1.7534813) q[2];
sx q[2];
rz(-1.8049392) q[2];
rz(3.0873599) q[3];
sx q[3];
rz(-1.1905866) q[3];
sx q[3];
rz(2.5469053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17276758) q[0];
sx q[0];
rz(-0.69197881) q[0];
sx q[0];
rz(-1.927595) q[0];
rz(-1.0460151) q[1];
sx q[1];
rz(-1.0449301) q[1];
sx q[1];
rz(-2.2878343) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19148286) q[0];
sx q[0];
rz(-1.5111418) q[0];
sx q[0];
rz(1.4364987) q[0];
rz(-pi) q[1];
rz(-2.3860903) q[2];
sx q[2];
rz(-1.1248338) q[2];
sx q[2];
rz(-1.310629) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5750372) q[1];
sx q[1];
rz(-2.2872529) q[1];
sx q[1];
rz(2.3385919) q[1];
x q[2];
rz(1.7231971) q[3];
sx q[3];
rz(-1.2255242) q[3];
sx q[3];
rz(2.4403646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.13009109) q[2];
sx q[2];
rz(-1.4902196) q[2];
sx q[2];
rz(-2.3940274) q[2];
rz(2.2085564) q[3];
sx q[3];
rz(-1.3676164) q[3];
sx q[3];
rz(-0.30195495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0316684) q[0];
sx q[0];
rz(-0.95369354) q[0];
sx q[0];
rz(-0.64144301) q[0];
rz(-0.96915069) q[1];
sx q[1];
rz(-1.0698003) q[1];
sx q[1];
rz(2.0279121) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4232491) q[0];
sx q[0];
rz(-1.166543) q[0];
sx q[0];
rz(3.0358311) q[0];
rz(-pi) q[1];
rz(0.93514438) q[2];
sx q[2];
rz(-1.1271659) q[2];
sx q[2];
rz(-2.005524) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.42150527) q[1];
sx q[1];
rz(-2.0893789) q[1];
sx q[1];
rz(-2.9353013) q[1];
rz(-pi) q[2];
rz(-0.05984743) q[3];
sx q[3];
rz(-0.52111292) q[3];
sx q[3];
rz(-1.9684362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.5668737) q[2];
sx q[2];
rz(-2.4428664) q[2];
sx q[2];
rz(-2.3835772) q[2];
rz(0.30592439) q[3];
sx q[3];
rz(-2.7829792) q[3];
sx q[3];
rz(2.6942286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4161943) q[0];
sx q[0];
rz(-2.5328126) q[0];
sx q[0];
rz(0.7533657) q[0];
rz(0.92195177) q[1];
sx q[1];
rz(-1.7148858) q[1];
sx q[1];
rz(-3.0070378) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20924231) q[0];
sx q[0];
rz(-1.9268231) q[0];
sx q[0];
rz(3.0491475) q[0];
rz(0.64132787) q[2];
sx q[2];
rz(-1.7761163) q[2];
sx q[2];
rz(0.63831282) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2462897) q[1];
sx q[1];
rz(-0.091769204) q[1];
sx q[1];
rz(0.20905881) q[1];
x q[2];
rz(-1.1128759) q[3];
sx q[3];
rz(-2.1919214) q[3];
sx q[3];
rz(0.21645138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.77384633) q[2];
sx q[2];
rz(-1.4464804) q[2];
sx q[2];
rz(-1.9473677) q[2];
rz(1.9959244) q[3];
sx q[3];
rz(-2.001389) q[3];
sx q[3];
rz(1.0660508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0399748) q[0];
sx q[0];
rz(-0.45926738) q[0];
sx q[0];
rz(-0.38247821) q[0];
rz(-2.3756012) q[1];
sx q[1];
rz(-1.4975558) q[1];
sx q[1];
rz(1.8006178) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0308834) q[0];
sx q[0];
rz(-1.481505) q[0];
sx q[0];
rz(-2.9813894) q[0];
x q[1];
rz(1.6587615) q[2];
sx q[2];
rz(-0.81552699) q[2];
sx q[2];
rz(1.2620827) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.080440259) q[1];
sx q[1];
rz(-2.3576479) q[1];
sx q[1];
rz(-2.817201) q[1];
x q[2];
rz(1.3519456) q[3];
sx q[3];
rz(-0.36389458) q[3];
sx q[3];
rz(-2.4879251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0535023) q[2];
sx q[2];
rz(-0.92038766) q[2];
sx q[2];
rz(3.0255393) q[2];
rz(-1.2608438) q[3];
sx q[3];
rz(-1.1382269) q[3];
sx q[3];
rz(1.457823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55105227) q[0];
sx q[0];
rz(-3.0278979) q[0];
sx q[0];
rz(-0.1846479) q[0];
rz(0.91839904) q[1];
sx q[1];
rz(-1.5469488) q[1];
sx q[1];
rz(-1.7041697) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0032063403) q[0];
sx q[0];
rz(-1.0995563) q[0];
sx q[0];
rz(0.22881656) q[0];
x q[1];
rz(-0.080731656) q[2];
sx q[2];
rz(-0.7237607) q[2];
sx q[2];
rz(-0.91659509) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.25333693) q[1];
sx q[1];
rz(-1.3789021) q[1];
sx q[1];
rz(-0.95805053) q[1];
x q[2];
rz(-0.67622306) q[3];
sx q[3];
rz(-1.3075324) q[3];
sx q[3];
rz(2.3674008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47716466) q[2];
sx q[2];
rz(-1.8958586) q[2];
sx q[2];
rz(2.5028382) q[2];
rz(-0.81513682) q[3];
sx q[3];
rz(-2.6705948) q[3];
sx q[3];
rz(-2.7208929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7246134) q[0];
sx q[0];
rz(-1.2495578) q[0];
sx q[0];
rz(0.15171262) q[0];
rz(-0.78085113) q[1];
sx q[1];
rz(-1.4518705) q[1];
sx q[1];
rz(-0.89444583) q[1];
rz(1.1679756) q[2];
sx q[2];
rz(-0.30277534) q[2];
sx q[2];
rz(1.4390611) q[2];
rz(-2.8975836) q[3];
sx q[3];
rz(-1.6088369) q[3];
sx q[3];
rz(3.101788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
