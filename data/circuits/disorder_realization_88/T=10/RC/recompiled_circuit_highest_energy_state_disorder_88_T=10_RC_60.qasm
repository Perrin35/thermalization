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
rz(4.881599) q[0];
sx q[0];
rz(10.344883) q[0];
rz(-0.080634557) q[1];
sx q[1];
rz(-2.562685) q[1];
sx q[1];
rz(0.97715598) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9465222) q[0];
sx q[0];
rz(-1.7818639) q[0];
sx q[0];
rz(2.7618183) q[0];
rz(-pi) q[1];
rz(1.5571655) q[2];
sx q[2];
rz(-1.7030099) q[2];
sx q[2];
rz(-1.4135828) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5861533) q[1];
sx q[1];
rz(-1.7057014) q[1];
sx q[1];
rz(-2.5389266) q[1];
rz(-pi) q[2];
rz(1.8771873) q[3];
sx q[3];
rz(-0.9054817) q[3];
sx q[3];
rz(3.0476335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5795634) q[2];
sx q[2];
rz(-1.2144438) q[2];
sx q[2];
rz(-0.62082949) q[2];
rz(-0.51554716) q[3];
sx q[3];
rz(-1.4225057) q[3];
sx q[3];
rz(0.8425042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8949378) q[0];
sx q[0];
rz(-1.5351013) q[0];
sx q[0];
rz(-0.57902336) q[0];
rz(-0.82194263) q[1];
sx q[1];
rz(-1.6252981) q[1];
sx q[1];
rz(2.6878405) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61868405) q[0];
sx q[0];
rz(-2.4755619) q[0];
sx q[0];
rz(-0.36938195) q[0];
rz(-pi) q[1];
rz(0.2208925) q[2];
sx q[2];
rz(-2.6651987) q[2];
sx q[2];
rz(-1.1645137) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.380341) q[1];
sx q[1];
rz(-0.56219343) q[1];
sx q[1];
rz(2.9060049) q[1];
x q[2];
rz(0.9048432) q[3];
sx q[3];
rz(-2.4840151) q[3];
sx q[3];
rz(-0.50298728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4072676) q[2];
sx q[2];
rz(-3.0027323) q[2];
sx q[2];
rz(-0.56488758) q[2];
rz(-2.7867553) q[3];
sx q[3];
rz(-0.9762888) q[3];
sx q[3];
rz(1.423665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(1.9722209) q[0];
sx q[0];
rz(-2.9300949) q[0];
sx q[0];
rz(2.5900904) q[0];
rz(2.7768199) q[1];
sx q[1];
rz(-0.33939895) q[1];
sx q[1];
rz(0.50484467) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0009815) q[0];
sx q[0];
rz(-1.4112844) q[0];
sx q[0];
rz(0.41054748) q[0];
rz(0.41270035) q[2];
sx q[2];
rz(-1.6566212) q[2];
sx q[2];
rz(2.7049261) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6712196) q[1];
sx q[1];
rz(-0.80802901) q[1];
sx q[1];
rz(0.68617448) q[1];
rz(-0.64845131) q[3];
sx q[3];
rz(-2.5621427) q[3];
sx q[3];
rz(3.0045829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7079033) q[2];
sx q[2];
rz(-1.7283551) q[2];
sx q[2];
rz(-0.047253963) q[2];
rz(-2.3047678) q[3];
sx q[3];
rz(-2.4182726) q[3];
sx q[3];
rz(-1.0251454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.717201) q[0];
sx q[0];
rz(-2.1121139) q[0];
sx q[0];
rz(-1.8151872) q[0];
rz(-1.6560417) q[1];
sx q[1];
rz(-1.0842208) q[1];
sx q[1];
rz(-0.63492376) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43760532) q[0];
sx q[0];
rz(-3.1235031) q[0];
sx q[0];
rz(1.6419069) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97135404) q[2];
sx q[2];
rz(-2.5581789) q[2];
sx q[2];
rz(2.3177878) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0628218) q[1];
sx q[1];
rz(-0.34343849) q[1];
sx q[1];
rz(1.354864) q[1];
x q[2];
rz(-0.62795378) q[3];
sx q[3];
rz(-2.0139004) q[3];
sx q[3];
rz(1.4301883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2403468) q[2];
sx q[2];
rz(-0.88902688) q[2];
sx q[2];
rz(-1.3725613) q[2];
rz(-1.9339804) q[3];
sx q[3];
rz(-2.4977081) q[3];
sx q[3];
rz(0.27455583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8127301) q[0];
sx q[0];
rz(-0.48506081) q[0];
sx q[0];
rz(-3.0364756) q[0];
rz(-0.33991995) q[1];
sx q[1];
rz(-0.9143908) q[1];
sx q[1];
rz(-2.0786659) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2539627) q[0];
sx q[0];
rz(-2.2821873) q[0];
sx q[0];
rz(1.5521868) q[0];
rz(-pi) q[1];
rz(0.7617217) q[2];
sx q[2];
rz(-1.7526547) q[2];
sx q[2];
rz(1.8363425) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.2984933) q[1];
sx q[1];
rz(-1.5122422) q[1];
sx q[1];
rz(-1.5565518) q[1];
rz(2.0696409) q[3];
sx q[3];
rz(-1.5728647) q[3];
sx q[3];
rz(-0.75132124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.062332705) q[2];
sx q[2];
rz(-1.7534813) q[2];
sx q[2];
rz(1.3366535) q[2];
rz(3.0873599) q[3];
sx q[3];
rz(-1.9510061) q[3];
sx q[3];
rz(0.59468734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9688251) q[0];
sx q[0];
rz(-2.4496138) q[0];
sx q[0];
rz(1.927595) q[0];
rz(2.0955775) q[1];
sx q[1];
rz(-2.0966625) q[1];
sx q[1];
rz(2.2878343) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19148286) q[0];
sx q[0];
rz(-1.5111418) q[0];
sx q[0];
rz(-1.7050939) q[0];
rz(2.3860903) q[2];
sx q[2];
rz(-2.0167588) q[2];
sx q[2];
rz(1.8309636) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5750372) q[1];
sx q[1];
rz(-2.2872529) q[1];
sx q[1];
rz(-2.3385919) q[1];
rz(-pi) q[2];
rz(1.7231971) q[3];
sx q[3];
rz(-1.9160685) q[3];
sx q[3];
rz(0.70122805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0115016) q[2];
sx q[2];
rz(-1.6513731) q[2];
sx q[2];
rz(-0.74756527) q[2];
rz(0.93303624) q[3];
sx q[3];
rz(-1.7739762) q[3];
sx q[3];
rz(-0.30195495) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1099243) q[0];
sx q[0];
rz(-2.1878991) q[0];
sx q[0];
rz(-2.5001496) q[0];
rz(2.172442) q[1];
sx q[1];
rz(-1.0698003) q[1];
sx q[1];
rz(-1.1136805) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1596368) q[0];
sx q[0];
rz(-2.7244719) q[0];
sx q[0];
rz(1.812716) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2064483) q[2];
sx q[2];
rz(-1.1271659) q[2];
sx q[2];
rz(-2.005524) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.021999849) q[1];
sx q[1];
rz(-2.5869859) q[1];
sx q[1];
rz(1.9153992) q[1];
rz(1.6051172) q[3];
sx q[3];
rz(-2.0908818) q[3];
sx q[3];
rz(1.8994562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.5668737) q[2];
sx q[2];
rz(-2.4428664) q[2];
sx q[2];
rz(-0.75801545) q[2];
rz(-0.30592439) q[3];
sx q[3];
rz(-0.35861349) q[3];
sx q[3];
rz(-0.44736403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-1.7253983) q[0];
sx q[0];
rz(-2.5328126) q[0];
sx q[0];
rz(-2.388227) q[0];
rz(0.92195177) q[1];
sx q[1];
rz(-1.7148858) q[1];
sx q[1];
rz(0.13455483) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9323503) q[0];
sx q[0];
rz(-1.9268231) q[0];
sx q[0];
rz(3.0491475) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64132787) q[2];
sx q[2];
rz(-1.7761163) q[2];
sx q[2];
rz(-0.63831282) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.53271097) q[1];
sx q[1];
rz(-1.5898165) q[1];
sx q[1];
rz(-0.089781922) q[1];
x q[2];
rz(0.67340019) q[3];
sx q[3];
rz(-1.2030501) q[3];
sx q[3];
rz(-1.0750225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3677463) q[2];
sx q[2];
rz(-1.4464804) q[2];
sx q[2];
rz(1.9473677) q[2];
rz(-1.1456683) q[3];
sx q[3];
rz(-1.1402036) q[3];
sx q[3];
rz(2.0755419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1016178) q[0];
sx q[0];
rz(-0.45926738) q[0];
sx q[0];
rz(0.38247821) q[0];
rz(-0.76599145) q[1];
sx q[1];
rz(-1.4975558) q[1];
sx q[1];
rz(1.3409748) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95582286) q[0];
sx q[0];
rz(-2.9583724) q[0];
sx q[0];
rz(-0.51143666) q[0];
rz(-pi) q[1];
rz(0.093042298) q[2];
sx q[2];
rz(-0.75936717) q[2];
sx q[2];
rz(2.0074646) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7786918) q[1];
sx q[1];
rz(-0.83759379) q[1];
sx q[1];
rz(-1.8785088) q[1];
rz(0.082499113) q[3];
sx q[3];
rz(-1.2159705) q[3];
sx q[3];
rz(-2.2542743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0880903) q[2];
sx q[2];
rz(-0.92038766) q[2];
sx q[2];
rz(-0.11605334) q[2];
rz(-1.2608438) q[3];
sx q[3];
rz(-1.1382269) q[3];
sx q[3];
rz(-1.6837696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.55105227) q[0];
sx q[0];
rz(-3.0278979) q[0];
sx q[0];
rz(2.9569448) q[0];
rz(2.2231936) q[1];
sx q[1];
rz(-1.5469488) q[1];
sx q[1];
rz(-1.437423) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47078314) q[0];
sx q[0];
rz(-2.6215186) q[0];
sx q[0];
rz(1.1519679) q[0];
rz(1.499648) q[2];
sx q[2];
rz(-2.2916823) q[2];
sx q[2];
rz(-0.80903731) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8882557) q[1];
sx q[1];
rz(-1.3789021) q[1];
sx q[1];
rz(-0.95805053) q[1];
rz(-pi) q[2];
rz(-0.67622306) q[3];
sx q[3];
rz(-1.3075324) q[3];
sx q[3];
rz(2.3674008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47716466) q[2];
sx q[2];
rz(-1.2457341) q[2];
sx q[2];
rz(2.5028382) q[2];
rz(-0.81513682) q[3];
sx q[3];
rz(-2.6705948) q[3];
sx q[3];
rz(0.42069978) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7246134) q[0];
sx q[0];
rz(-1.8920349) q[0];
sx q[0];
rz(-2.98988) q[0];
rz(2.3607415) q[1];
sx q[1];
rz(-1.4518705) q[1];
sx q[1];
rz(-0.89444583) q[1];
rz(3.0197418) q[2];
sx q[2];
rz(-1.8486628) q[2];
sx q[2];
rz(1.8589414) q[2];
rz(0.15624795) q[3];
sx q[3];
rz(-2.8946946) q[3];
sx q[3];
rz(1.3794086) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
