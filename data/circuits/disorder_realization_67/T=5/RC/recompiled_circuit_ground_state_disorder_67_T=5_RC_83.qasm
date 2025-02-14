OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.23842731) q[0];
sx q[0];
rz(-2.981346) q[0];
sx q[0];
rz(1.2877553) q[0];
rz(2.2924478) q[1];
sx q[1];
rz(-1.949911) q[1];
sx q[1];
rz(-0.75872672) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903126) q[0];
sx q[0];
rz(-2.7669124) q[0];
sx q[0];
rz(-0.6775585) q[0];
rz(-pi) q[1];
rz(-2.8992462) q[2];
sx q[2];
rz(-1.0966612) q[2];
sx q[2];
rz(2.989492) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1114486) q[1];
sx q[1];
rz(-1.9322825) q[1];
sx q[1];
rz(0.8059146) q[1];
x q[2];
rz(-0.87293245) q[3];
sx q[3];
rz(-2.2651787) q[3];
sx q[3];
rz(2.5257655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8410926) q[2];
sx q[2];
rz(-0.79885834) q[2];
sx q[2];
rz(-0.50660261) q[2];
rz(-1.6182342) q[3];
sx q[3];
rz(-0.62464276) q[3];
sx q[3];
rz(0.055559572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27187207) q[0];
sx q[0];
rz(-0.46590081) q[0];
sx q[0];
rz(3.077935) q[0];
rz(-1.1977389) q[1];
sx q[1];
rz(-1.5143737) q[1];
sx q[1];
rz(1.0187842) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2301529) q[0];
sx q[0];
rz(-2.8960315) q[0];
sx q[0];
rz(-2.5018238) q[0];
x q[1];
rz(-0.59661023) q[2];
sx q[2];
rz(-2.7368174) q[2];
sx q[2];
rz(-0.34133729) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9460822) q[1];
sx q[1];
rz(-1.2942593) q[1];
sx q[1];
rz(-0.72317381) q[1];
rz(-1.7131931) q[3];
sx q[3];
rz(-2.3835858) q[3];
sx q[3];
rz(2.6482794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7055052) q[2];
sx q[2];
rz(-0.59385308) q[2];
sx q[2];
rz(-0.84211055) q[2];
rz(1.0574794) q[3];
sx q[3];
rz(-2.2130241) q[3];
sx q[3];
rz(-1.1697945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8975526) q[0];
sx q[0];
rz(-2.121448) q[0];
sx q[0];
rz(-1.1919588) q[0];
rz(0.92487088) q[1];
sx q[1];
rz(-2.4332739) q[1];
sx q[1];
rz(-1.8398197) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10629207) q[0];
sx q[0];
rz(-2.3791693) q[0];
sx q[0];
rz(-0.85233217) q[0];
x q[1];
rz(0.19916735) q[2];
sx q[2];
rz(-1.5323497) q[2];
sx q[2];
rz(3.0035915) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.30693248) q[1];
sx q[1];
rz(-1.7494769) q[1];
sx q[1];
rz(-0.4731005) q[1];
rz(-pi) q[2];
rz(-1.916094) q[3];
sx q[3];
rz(-1.32816) q[3];
sx q[3];
rz(0.43365989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7630345) q[2];
sx q[2];
rz(-0.083746567) q[2];
sx q[2];
rz(3.0540826) q[2];
rz(-1.8267501) q[3];
sx q[3];
rz(-1.0564691) q[3];
sx q[3];
rz(-2.1480613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73115504) q[0];
sx q[0];
rz(-1.9696099) q[0];
sx q[0];
rz(-2.3325969) q[0];
rz(-1.0357098) q[1];
sx q[1];
rz(-2.6782942) q[1];
sx q[1];
rz(2.1083924) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.746377) q[0];
sx q[0];
rz(-0.6653293) q[0];
sx q[0];
rz(1.6341342) q[0];
rz(-pi) q[1];
rz(-1.0339917) q[2];
sx q[2];
rz(-1.189917) q[2];
sx q[2];
rz(-1.2950667) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.08398) q[1];
sx q[1];
rz(-1.3911493) q[1];
sx q[1];
rz(2.6704726) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71113385) q[3];
sx q[3];
rz(-2.4681849) q[3];
sx q[3];
rz(0.26949689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0741299) q[2];
sx q[2];
rz(-2.9758657) q[2];
sx q[2];
rz(1.5376252) q[2];
rz(1.4175203) q[3];
sx q[3];
rz(-1.1105024) q[3];
sx q[3];
rz(1.0874776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9017225) q[0];
sx q[0];
rz(-2.9481695) q[0];
sx q[0];
rz(-1.1892009) q[0];
rz(-0.72422782) q[1];
sx q[1];
rz(-1.291357) q[1];
sx q[1];
rz(1.6667746) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9928707) q[0];
sx q[0];
rz(-1.9970511) q[0];
sx q[0];
rz(0.54963407) q[0];
rz(-pi) q[1];
rz(1.4101661) q[2];
sx q[2];
rz(-2.361627) q[2];
sx q[2];
rz(0.51973625) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.96383038) q[1];
sx q[1];
rz(-0.46198598) q[1];
sx q[1];
rz(-0.83779184) q[1];
rz(-2.0678364) q[3];
sx q[3];
rz(-1.1599564) q[3];
sx q[3];
rz(2.2548294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7597947) q[2];
sx q[2];
rz(-2.2319824) q[2];
sx q[2];
rz(-2.6749715) q[2];
rz(1.4383379) q[3];
sx q[3];
rz(-1.2983026) q[3];
sx q[3];
rz(-0.94021016) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4284215) q[0];
sx q[0];
rz(-0.82601014) q[0];
sx q[0];
rz(3.1312422) q[0];
rz(-1.5230644) q[1];
sx q[1];
rz(-1.3925939) q[1];
sx q[1];
rz(1.7656322) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3288157) q[0];
sx q[0];
rz(-1.8121908) q[0];
sx q[0];
rz(1.2230948) q[0];
rz(-pi) q[1];
rz(0.1907804) q[2];
sx q[2];
rz(-1.6990802) q[2];
sx q[2];
rz(3.0012263) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.048307) q[1];
sx q[1];
rz(-1.3897589) q[1];
sx q[1];
rz(-2.6345318) q[1];
x q[2];
rz(0.74919219) q[3];
sx q[3];
rz(-1.0674879) q[3];
sx q[3];
rz(-1.0096514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7213664) q[2];
sx q[2];
rz(-2.0927636) q[2];
sx q[2];
rz(0.72727195) q[2];
rz(-0.51491245) q[3];
sx q[3];
rz(-1.3313096) q[3];
sx q[3];
rz(-1.7236727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10776831) q[0];
sx q[0];
rz(-1.6078147) q[0];
sx q[0];
rz(-2.7346101) q[0];
rz(2.175323) q[1];
sx q[1];
rz(-0.85710183) q[1];
sx q[1];
rz(0.13872096) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16628597) q[0];
sx q[0];
rz(-1.606267) q[0];
sx q[0];
rz(1.4544288) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7916324) q[2];
sx q[2];
rz(-1.2996309) q[2];
sx q[2];
rz(1.8200995) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2695864) q[1];
sx q[1];
rz(-1.9161738) q[1];
sx q[1];
rz(-0.38265444) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2267031) q[3];
sx q[3];
rz(-1.4315245) q[3];
sx q[3];
rz(0.95668281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8653284) q[2];
sx q[2];
rz(-1.3778957) q[2];
sx q[2];
rz(2.830937) q[2];
rz(1.3079414) q[3];
sx q[3];
rz(-1.5111978) q[3];
sx q[3];
rz(-0.37180296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0528316) q[0];
sx q[0];
rz(-2.796687) q[0];
sx q[0];
rz(-0.088223591) q[0];
rz(-3.131033) q[1];
sx q[1];
rz(-1.6190395) q[1];
sx q[1];
rz(-2.6710076) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7475766) q[0];
sx q[0];
rz(-2.2758099) q[0];
sx q[0];
rz(2.2590619) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.768465) q[2];
sx q[2];
rz(-0.77357793) q[2];
sx q[2];
rz(0.3244704) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1697929) q[1];
sx q[1];
rz(-1.3108675) q[1];
sx q[1];
rz(3.0539576) q[1];
rz(0.13981522) q[3];
sx q[3];
rz(-1.2212409) q[3];
sx q[3];
rz(-0.67963723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.18870246) q[2];
sx q[2];
rz(-2.1850093) q[2];
sx q[2];
rz(2.0616222) q[2];
rz(-2.2522669) q[3];
sx q[3];
rz(-1.5998799) q[3];
sx q[3];
rz(-1.5573474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66460669) q[0];
sx q[0];
rz(-2.7598858) q[0];
sx q[0];
rz(-2.3241924) q[0];
rz(-0.74642247) q[1];
sx q[1];
rz(-2.4669929) q[1];
sx q[1];
rz(-2.2019763) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25519263) q[0];
sx q[0];
rz(-1.5879244) q[0];
sx q[0];
rz(-3.1411132) q[0];
rz(-pi) q[1];
x q[1];
rz(1.290326) q[2];
sx q[2];
rz(-1.2277516) q[2];
sx q[2];
rz(-1.5671935) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.47245644) q[1];
sx q[1];
rz(-1.8207153) q[1];
sx q[1];
rz(0.12052287) q[1];
rz(-pi) q[2];
rz(2.0321531) q[3];
sx q[3];
rz(-0.45518866) q[3];
sx q[3];
rz(-0.92092848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.32999906) q[2];
sx q[2];
rz(-0.96572319) q[2];
sx q[2];
rz(-2.7962371) q[2];
rz(1.5251478) q[3];
sx q[3];
rz(-1.7914146) q[3];
sx q[3];
rz(2.7143872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3578607) q[0];
sx q[0];
rz(-1.1564199) q[0];
sx q[0];
rz(-0.082948908) q[0];
rz(0.16231617) q[1];
sx q[1];
rz(-1.6971308) q[1];
sx q[1];
rz(-0.69581318) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0950553) q[0];
sx q[0];
rz(-2.8588606) q[0];
sx q[0];
rz(-2.9912492) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5478304) q[2];
sx q[2];
rz(-3.0122979) q[2];
sx q[2];
rz(-0.33516075) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4033234) q[1];
sx q[1];
rz(-1.5517071) q[1];
sx q[1];
rz(0.44428579) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70503321) q[3];
sx q[3];
rz(-1.789973) q[3];
sx q[3];
rz(2.8104643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7210228) q[2];
sx q[2];
rz(-0.15779725) q[2];
sx q[2];
rz(1.921462) q[2];
rz(-2.5681791) q[3];
sx q[3];
rz(-1.3308595) q[3];
sx q[3];
rz(-2.8964608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86380105) q[0];
sx q[0];
rz(-1.6391123) q[0];
sx q[0];
rz(3.0085051) q[0];
rz(-0.0028903891) q[1];
sx q[1];
rz(-2.7255701) q[1];
sx q[1];
rz(2.4400673) q[1];
rz(-2.9691545) q[2];
sx q[2];
rz(-1.8159096) q[2];
sx q[2];
rz(2.590948) q[2];
rz(-1.6297798) q[3];
sx q[3];
rz(-0.77131693) q[3];
sx q[3];
rz(-0.49280096) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
