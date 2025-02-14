OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.339191) q[0];
sx q[0];
rz(-1.1561166) q[0];
sx q[0];
rz(-2.6240786) q[0];
rz(1.327688) q[1];
sx q[1];
rz(7.1751243) q[1];
sx q[1];
rz(8.8827477) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5919246) q[0];
sx q[0];
rz(-0.79093638) q[0];
sx q[0];
rz(2.8998081) q[0];
rz(-pi) q[1];
rz(0.26371358) q[2];
sx q[2];
rz(-2.3576735) q[2];
sx q[2];
rz(-0.5413407) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.65329018) q[1];
sx q[1];
rz(-0.79217109) q[1];
sx q[1];
rz(1.2912871) q[1];
x q[2];
rz(0.63707085) q[3];
sx q[3];
rz(-2.8544997) q[3];
sx q[3];
rz(1.4305103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1338256) q[2];
sx q[2];
rz(-1.0453036) q[2];
sx q[2];
rz(-2.2110151) q[2];
rz(3.0104356) q[3];
sx q[3];
rz(-0.37808642) q[3];
sx q[3];
rz(-1.650943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59830484) q[0];
sx q[0];
rz(-2.2853993) q[0];
sx q[0];
rz(1.244586) q[0];
rz(2.2159131) q[1];
sx q[1];
rz(-1.2558179) q[1];
sx q[1];
rz(-1.6474887) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6005352) q[0];
sx q[0];
rz(-2.3575258) q[0];
sx q[0];
rz(-1.1009786) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7297614) q[2];
sx q[2];
rz(-0.71626012) q[2];
sx q[2];
rz(-2.0873983) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0706961) q[1];
sx q[1];
rz(-2.5578209) q[1];
sx q[1];
rz(0.6604922) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5794287) q[3];
sx q[3];
rz(-1.2666456) q[3];
sx q[3];
rz(0.22316531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.99574789) q[2];
sx q[2];
rz(-2.6971942) q[2];
sx q[2];
rz(-2.2399529) q[2];
rz(1.7835279) q[3];
sx q[3];
rz(-1.8243022) q[3];
sx q[3];
rz(-2.5942514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3667592) q[0];
sx q[0];
rz(-1.1792553) q[0];
sx q[0];
rz(1.9955848) q[0];
rz(0.92578069) q[1];
sx q[1];
rz(-2.130276) q[1];
sx q[1];
rz(-2.944223) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5615047) q[0];
sx q[0];
rz(-2.5338123) q[0];
sx q[0];
rz(1.8124026) q[0];
x q[1];
rz(-1.2659129) q[2];
sx q[2];
rz(-1.8499057) q[2];
sx q[2];
rz(0.32372083) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.47599736) q[1];
sx q[1];
rz(-1.531812) q[1];
sx q[1];
rz(-1.798756) q[1];
rz(-pi) q[2];
rz(-0.32549627) q[3];
sx q[3];
rz(-2.5739658) q[3];
sx q[3];
rz(-2.5095255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4027412) q[2];
sx q[2];
rz(-0.98029843) q[2];
sx q[2];
rz(0.19035467) q[2];
rz(3.0800152) q[3];
sx q[3];
rz(-0.96934167) q[3];
sx q[3];
rz(2.0891345) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9848118) q[0];
sx q[0];
rz(-1.7233912) q[0];
sx q[0];
rz(-0.62830997) q[0];
rz(1.620694) q[1];
sx q[1];
rz(-1.9338806) q[1];
sx q[1];
rz(0.77888387) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6464061) q[0];
sx q[0];
rz(-0.71346006) q[0];
sx q[0];
rz(0.97795217) q[0];
x q[1];
rz(0.36036307) q[2];
sx q[2];
rz(-1.833263) q[2];
sx q[2];
rz(1.2059463) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6538222) q[1];
sx q[1];
rz(-1.5573335) q[1];
sx q[1];
rz(3.0938992) q[1];
rz(-pi) q[2];
x q[2];
rz(0.811399) q[3];
sx q[3];
rz(-2.2795942) q[3];
sx q[3];
rz(0.78042316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6598307) q[2];
sx q[2];
rz(-2.063282) q[2];
sx q[2];
rz(3.1415494) q[2];
rz(2.0569885) q[3];
sx q[3];
rz(-1.9856039) q[3];
sx q[3];
rz(0.46561766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67149177) q[0];
sx q[0];
rz(-1.6974314) q[0];
sx q[0];
rz(-1.8430365) q[0];
rz(-0.57688722) q[1];
sx q[1];
rz(-1.2737609) q[1];
sx q[1];
rz(-2.6947122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6046238) q[0];
sx q[0];
rz(-0.67567053) q[0];
sx q[0];
rz(-2.8506822) q[0];
rz(1.8700897) q[2];
sx q[2];
rz(-1.4993877) q[2];
sx q[2];
rz(-1.4324607) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.943534) q[1];
sx q[1];
rz(-1.3683649) q[1];
sx q[1];
rz(-1.979503) q[1];
x q[2];
rz(-1.6303667) q[3];
sx q[3];
rz(-1.7559663) q[3];
sx q[3];
rz(2.6310754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.4345066) q[2];
sx q[2];
rz(-2.862317) q[2];
sx q[2];
rz(-2.9217829) q[2];
rz(1.6541727) q[3];
sx q[3];
rz(-1.4498962) q[3];
sx q[3];
rz(-0.69948227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18911067) q[0];
sx q[0];
rz(-0.43287745) q[0];
sx q[0];
rz(-2.2440946) q[0];
rz(-1.7706361) q[1];
sx q[1];
rz(-2.1015344) q[1];
sx q[1];
rz(0.8955566) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61293759) q[0];
sx q[0];
rz(-2.2870758) q[0];
sx q[0];
rz(1.7578633) q[0];
rz(-pi) q[1];
rz(0.61195749) q[2];
sx q[2];
rz(-1.0494119) q[2];
sx q[2];
rz(2.251791) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1119536) q[1];
sx q[1];
rz(-1.1983419) q[1];
sx q[1];
rz(1.2122494) q[1];
rz(-pi) q[2];
rz(2.343781) q[3];
sx q[3];
rz(-0.63029248) q[3];
sx q[3];
rz(0.52784656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7728277) q[2];
sx q[2];
rz(-2.0389098) q[2];
sx q[2];
rz(0.15711288) q[2];
rz(-2.469192) q[3];
sx q[3];
rz(-1.3998569) q[3];
sx q[3];
rz(-0.11252832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86033487) q[0];
sx q[0];
rz(-0.66547886) q[0];
sx q[0];
rz(-1.459664) q[0];
rz(2.7809987) q[1];
sx q[1];
rz(-1.4692042) q[1];
sx q[1];
rz(-1.8416539) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4111709) q[0];
sx q[0];
rz(-1.9864559) q[0];
sx q[0];
rz(1.4372197) q[0];
rz(-pi) q[1];
rz(2.0327914) q[2];
sx q[2];
rz(-1.0014152) q[2];
sx q[2];
rz(1.580736) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9323401) q[1];
sx q[1];
rz(-0.31267088) q[1];
sx q[1];
rz(2.8220842) q[1];
rz(-0.54706562) q[3];
sx q[3];
rz(-0.69098847) q[3];
sx q[3];
rz(0.83741659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2099057) q[2];
sx q[2];
rz(-2.1114712) q[2];
sx q[2];
rz(0.21200655) q[2];
rz(1.0509521) q[3];
sx q[3];
rz(-1.7651599) q[3];
sx q[3];
rz(-0.088137805) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5365005) q[0];
sx q[0];
rz(-2.3516646) q[0];
sx q[0];
rz(2.9686046) q[0];
rz(1.1146924) q[1];
sx q[1];
rz(-2.098691) q[1];
sx q[1];
rz(-3.0317422) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6914586) q[0];
sx q[0];
rz(-0.5930674) q[0];
sx q[0];
rz(-2.791138) q[0];
rz(2.3187056) q[2];
sx q[2];
rz(-1.7639065) q[2];
sx q[2];
rz(-0.88179526) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6348656) q[1];
sx q[1];
rz(-2.8322729) q[1];
sx q[1];
rz(0.97815467) q[1];
rz(-pi) q[2];
rz(-0.55446378) q[3];
sx q[3];
rz(-1.672847) q[3];
sx q[3];
rz(-2.3481675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7383808) q[2];
sx q[2];
rz(-1.9387559) q[2];
sx q[2];
rz(0.66486764) q[2];
rz(-0.46857771) q[3];
sx q[3];
rz(-1.6178308) q[3];
sx q[3];
rz(2.328228) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68798962) q[0];
sx q[0];
rz(-0.60992321) q[0];
sx q[0];
rz(-2.4000121) q[0];
rz(1.954151) q[1];
sx q[1];
rz(-1.5267742) q[1];
sx q[1];
rz(0.56914079) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81108196) q[0];
sx q[0];
rz(-1.2449322) q[0];
sx q[0];
rz(1.261607) q[0];
rz(0.14573614) q[2];
sx q[2];
rz(-0.45782858) q[2];
sx q[2];
rz(-2.482058) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5829825) q[1];
sx q[1];
rz(-0.64839586) q[1];
sx q[1];
rz(-2.3984539) q[1];
rz(-pi) q[2];
rz(0.82707246) q[3];
sx q[3];
rz(-0.88445348) q[3];
sx q[3];
rz(2.4857869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.49447122) q[2];
sx q[2];
rz(-2.1529866) q[2];
sx q[2];
rz(2.3753601) q[2];
rz(1.1726441) q[3];
sx q[3];
rz(-1.6021043) q[3];
sx q[3];
rz(-2.9197781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2072993) q[0];
sx q[0];
rz(-2.793029) q[0];
sx q[0];
rz(-2.9497414) q[0];
rz(2.3244997) q[1];
sx q[1];
rz(-2.8849738) q[1];
sx q[1];
rz(-2.4339035) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1199214) q[0];
sx q[0];
rz(-1.7448398) q[0];
sx q[0];
rz(-2.8928323) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0470263) q[2];
sx q[2];
rz(-2.4942538) q[2];
sx q[2];
rz(0.29905427) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7308049) q[1];
sx q[1];
rz(-1.2025226) q[1];
sx q[1];
rz(0.94502298) q[1];
rz(-pi) q[2];
rz(-2.3419699) q[3];
sx q[3];
rz(-0.76086894) q[3];
sx q[3];
rz(1.7590673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2231458) q[2];
sx q[2];
rz(-0.94474363) q[2];
sx q[2];
rz(0.72500149) q[2];
rz(-0.21507344) q[3];
sx q[3];
rz(-0.30078617) q[3];
sx q[3];
rz(-2.422629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216777) q[0];
sx q[0];
rz(-2.324993) q[0];
sx q[0];
rz(-2.8391229) q[0];
rz(1.1014145) q[1];
sx q[1];
rz(-1.9093724) q[1];
sx q[1];
rz(-2.2473635) q[1];
rz(-1.7673372) q[2];
sx q[2];
rz(-2.5977618) q[2];
sx q[2];
rz(0.14731461) q[2];
rz(1.5433031) q[3];
sx q[3];
rz(-1.1034272) q[3];
sx q[3];
rz(-2.6381794) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
