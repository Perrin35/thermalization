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
rz(-0.92010486) q[0];
rz(-0.080634557) q[1];
sx q[1];
rz(3.7205003) q[1];
sx q[1];
rz(10.401934) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8590737) q[0];
sx q[0];
rz(-0.4319829) q[0];
sx q[0];
rz(2.6175314) q[0];
rz(1.5844272) q[2];
sx q[2];
rz(-1.7030099) q[2];
sx q[2];
rz(1.4135828) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0339573) q[1];
sx q[1];
rz(-2.1672241) q[1];
sx q[1];
rz(1.7340842) q[1];
x q[2];
rz(-0.36698384) q[3];
sx q[3];
rz(-0.72260783) q[3];
sx q[3];
rz(2.5740576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5795634) q[2];
sx q[2];
rz(-1.9271489) q[2];
sx q[2];
rz(0.62082949) q[2];
rz(-0.51554716) q[3];
sx q[3];
rz(-1.4225057) q[3];
sx q[3];
rz(0.8425042) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8949378) q[0];
sx q[0];
rz(-1.6064914) q[0];
sx q[0];
rz(-0.57902336) q[0];
rz(2.31965) q[1];
sx q[1];
rz(-1.5162946) q[1];
sx q[1];
rz(-2.6878405) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0653719) q[0];
sx q[0];
rz(-2.1848618) q[0];
sx q[0];
rz(1.2943511) q[0];
rz(-0.2208925) q[2];
sx q[2];
rz(-0.476394) q[2];
sx q[2];
rz(1.9770789) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4848043) q[1];
sx q[1];
rz(-2.1156807) q[1];
sx q[1];
rz(-1.4247895) q[1];
rz(-pi) q[2];
rz(0.445153) q[3];
sx q[3];
rz(-2.0721407) q[3];
sx q[3];
rz(-1.2847163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4072676) q[2];
sx q[2];
rz(-3.0027323) q[2];
sx q[2];
rz(0.56488758) q[2];
rz(0.35483739) q[3];
sx q[3];
rz(-2.1653039) q[3];
sx q[3];
rz(1.7179276) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
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
rz(-2.8021937) q[1];
sx q[1];
rz(2.636748) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0802949) q[0];
sx q[0];
rz(-2.7027931) q[0];
sx q[0];
rz(2.7584235) q[0];
rz(-pi) q[1];
rz(0.21130685) q[2];
sx q[2];
rz(-0.42102764) q[2];
sx q[2];
rz(2.2007934) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.4144812) q[1];
sx q[1];
rz(-2.046578) q[1];
sx q[1];
rz(-0.68051312) q[1];
rz(-2.6608493) q[3];
sx q[3];
rz(-1.2337483) q[3];
sx q[3];
rz(-2.2729006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7079033) q[2];
sx q[2];
rz(-1.4132376) q[2];
sx q[2];
rz(3.0943387) q[2];
rz(0.83682483) q[3];
sx q[3];
rz(-0.72332007) q[3];
sx q[3];
rz(1.0251454) q[3];
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
rz(-1.717201) q[0];
sx q[0];
rz(-2.1121139) q[0];
sx q[0];
rz(1.3264054) q[0];
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
rz(-2.7039873) q[0];
sx q[0];
rz(-3.1235031) q[0];
sx q[0];
rz(1.6419069) q[0];
rz(-pi) q[1];
rz(-2.1702386) q[2];
sx q[2];
rz(-0.58341372) q[2];
sx q[2];
rz(2.3177878) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.078770854) q[1];
sx q[1];
rz(-0.34343849) q[1];
sx q[1];
rz(1.7867286) q[1];
rz(-pi) q[2];
rz(2.4621099) q[3];
sx q[3];
rz(-2.3906997) q[3];
sx q[3];
rz(0.39284947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9012458) q[2];
sx q[2];
rz(-0.88902688) q[2];
sx q[2];
rz(-1.7690313) q[2];
rz(1.2076123) q[3];
sx q[3];
rz(-2.4977081) q[3];
sx q[3];
rz(-2.8670368) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(0.10511705) q[0];
rz(2.8016727) q[1];
sx q[1];
rz(-2.2272019) q[1];
sx q[1];
rz(2.0786659) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8369097) q[0];
sx q[0];
rz(-1.5848918) q[0];
sx q[0];
rz(-0.71147664) q[0];
rz(-pi) q[1];
rz(1.3219484) q[2];
sx q[2];
rz(-2.3169059) q[2];
sx q[2];
rz(-2.7052372) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2731367) q[1];
sx q[1];
rz(-1.5565762) q[1];
sx q[1];
rz(-3.0830326) q[1];
rz(-pi) q[2];
rz(1.0719518) q[3];
sx q[3];
rz(-1.5728647) q[3];
sx q[3];
rz(0.75132124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.062332705) q[2];
sx q[2];
rz(-1.3881114) q[2];
sx q[2];
rz(1.3366535) q[2];
rz(-3.0873599) q[3];
sx q[3];
rz(-1.9510061) q[3];
sx q[3];
rz(-0.59468734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9688251) q[0];
sx q[0];
rz(-0.69197881) q[0];
sx q[0];
rz(1.2139976) q[0];
rz(-2.0955775) q[1];
sx q[1];
rz(-1.0449301) q[1];
sx q[1];
rz(2.2878343) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7542242) q[0];
sx q[0];
rz(-1.7048536) q[0];
sx q[0];
rz(-3.0813974) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5327024) q[2];
sx q[2];
rz(-0.85431803) q[2];
sx q[2];
rz(2.4520055) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5665555) q[1];
sx q[1];
rz(-2.2872529) q[1];
sx q[1];
rz(-2.3385919) q[1];
rz(0.3993897) q[3];
sx q[3];
rz(-0.37617427) q[3];
sx q[3];
rz(-2.0143505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.13009109) q[2];
sx q[2];
rz(-1.6513731) q[2];
sx q[2];
rz(0.74756527) q[2];
rz(0.93303624) q[3];
sx q[3];
rz(-1.3676164) q[3];
sx q[3];
rz(-2.8396377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1099243) q[0];
sx q[0];
rz(-0.95369354) q[0];
sx q[0];
rz(0.64144301) q[0];
rz(-0.96915069) q[1];
sx q[1];
rz(-2.0717924) q[1];
sx q[1];
rz(-2.0279121) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4232491) q[0];
sx q[0];
rz(-1.9750496) q[0];
sx q[0];
rz(0.1057616) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2458057) q[2];
sx q[2];
rz(-0.75715827) q[2];
sx q[2];
rz(0.96162187) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1195928) q[1];
sx q[1];
rz(-2.5869859) q[1];
sx q[1];
rz(1.9153992) q[1];
rz(-pi) q[2];
rz(-0.05984743) q[3];
sx q[3];
rz(-0.52111292) q[3];
sx q[3];
rz(1.1731565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.574719) q[2];
sx q[2];
rz(-0.69872624) q[2];
sx q[2];
rz(-2.3835772) q[2];
rz(-2.8356683) q[3];
sx q[3];
rz(-0.35861349) q[3];
sx q[3];
rz(-2.6942286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7253983) q[0];
sx q[0];
rz(-0.60878009) q[0];
sx q[0];
rz(-2.388227) q[0];
rz(-0.92195177) q[1];
sx q[1];
rz(-1.4267068) q[1];
sx q[1];
rz(0.13455483) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9323503) q[0];
sx q[0];
rz(-1.2147696) q[0];
sx q[0];
rz(-3.0491475) q[0];
x q[1];
rz(1.3165264) q[2];
sx q[2];
rz(-0.94506028) q[2];
sx q[2];
rz(0.78142399) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.53271097) q[1];
sx q[1];
rz(-1.5517762) q[1];
sx q[1];
rz(-3.0518107) q[1];
x q[2];
rz(0.67340019) q[3];
sx q[3];
rz(-1.2030501) q[3];
sx q[3];
rz(-1.0750225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.77384633) q[2];
sx q[2];
rz(-1.4464804) q[2];
sx q[2];
rz(-1.194225) q[2];
rz(1.9959244) q[3];
sx q[3];
rz(-1.1402036) q[3];
sx q[3];
rz(2.0755419) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0399748) q[0];
sx q[0];
rz(-2.6823253) q[0];
sx q[0];
rz(-2.7591144) q[0];
rz(2.3756012) q[1];
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
rz(2.630156) q[0];
x q[1];
rz(-3.0485504) q[2];
sx q[2];
rz(-2.3822255) q[2];
sx q[2];
rz(1.134128) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0611524) q[1];
sx q[1];
rz(-0.78394475) q[1];
sx q[1];
rz(0.32439167) q[1];
rz(-pi) q[2];
rz(1.9267335) q[3];
sx q[3];
rz(-1.4934469) q[3];
sx q[3];
rz(0.71219769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0880903) q[2];
sx q[2];
rz(-2.221205) q[2];
sx q[2];
rz(0.11605334) q[2];
rz(1.2608438) q[3];
sx q[3];
rz(-1.1382269) q[3];
sx q[3];
rz(1.6837696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55105227) q[0];
sx q[0];
rz(-3.0278979) q[0];
sx q[0];
rz(-0.1846479) q[0];
rz(-2.2231936) q[1];
sx q[1];
rz(-1.5469488) q[1];
sx q[1];
rz(-1.7041697) q[1];
rz(pi/2) q[2];
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
rz(3.060861) q[2];
sx q[2];
rz(-2.417832) q[2];
sx q[2];
rz(0.91659509) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5823707) q[1];
sx q[1];
rz(-2.5032024) q[1];
sx q[1];
rz(-1.2450144) q[1];
x q[2];
rz(-2.4653696) q[3];
sx q[3];
rz(-1.8340602) q[3];
sx q[3];
rz(-0.77419188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.664428) q[2];
sx q[2];
rz(-1.8958586) q[2];
sx q[2];
rz(-2.5028382) q[2];
rz(-0.81513682) q[3];
sx q[3];
rz(-2.6705948) q[3];
sx q[3];
rz(0.42069978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
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
rz(1.9736171) q[2];
sx q[2];
rz(-2.8388173) q[2];
sx q[2];
rz(-1.7025316) q[2];
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
