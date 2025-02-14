OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.1240368) q[0];
sx q[0];
rz(-0.15260829) q[0];
sx q[0];
rz(-0.18327644) q[0];
rz(0.89219379) q[1];
sx q[1];
rz(-1.1396989) q[1];
sx q[1];
rz(-0.5782063) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86872549) q[0];
sx q[0];
rz(-1.5807165) q[0];
sx q[0];
rz(0.091477576) q[0];
rz(-pi) q[1];
rz(0.033138795) q[2];
sx q[2];
rz(-1.4862747) q[2];
sx q[2];
rz(-0.13730857) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.86204927) q[1];
sx q[1];
rz(-2.0506564) q[1];
sx q[1];
rz(-1.0399489) q[1];
rz(-pi) q[2];
rz(-0.036567612) q[3];
sx q[3];
rz(-1.1053549) q[3];
sx q[3];
rz(2.5643666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1766498) q[2];
sx q[2];
rz(-2.2735333) q[2];
sx q[2];
rz(0.098026015) q[2];
rz(2.3227504) q[3];
sx q[3];
rz(-0.10418532) q[3];
sx q[3];
rz(0.63481832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1408511) q[0];
sx q[0];
rz(-2.8241557) q[0];
sx q[0];
rz(0.66824085) q[0];
rz(0.60403281) q[1];
sx q[1];
rz(-0.030345358) q[1];
sx q[1];
rz(-2.277453) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5174422) q[0];
sx q[0];
rz(-1.6675341) q[0];
sx q[0];
rz(2.5521432) q[0];
rz(-pi) q[1];
rz(1.0572724) q[2];
sx q[2];
rz(-2.394426) q[2];
sx q[2];
rz(-0.42078373) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8616903) q[1];
sx q[1];
rz(-2.8331625) q[1];
sx q[1];
rz(0.66628404) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20012466) q[3];
sx q[3];
rz(-1.311655) q[3];
sx q[3];
rz(0.6249431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.16022564) q[2];
sx q[2];
rz(-1.9619433) q[2];
sx q[2];
rz(0.78440624) q[2];
rz(1.9085599) q[3];
sx q[3];
rz(-0.27850702) q[3];
sx q[3];
rz(-1.2109141) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3311555) q[0];
sx q[0];
rz(-1.5657319) q[0];
sx q[0];
rz(1.021215) q[0];
rz(-1.637623) q[1];
sx q[1];
rz(-2.8228788) q[1];
sx q[1];
rz(-0.57600299) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2191781) q[0];
sx q[0];
rz(-1.6873345) q[0];
sx q[0];
rz(2.7228505) q[0];
x q[1];
rz(2.8852413) q[2];
sx q[2];
rz(-1.1927989) q[2];
sx q[2];
rz(0.5452273) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0124023) q[1];
sx q[1];
rz(-2.2134212) q[1];
sx q[1];
rz(1.1430278) q[1];
rz(-pi) q[2];
rz(-1.3542451) q[3];
sx q[3];
rz(-1.3050021) q[3];
sx q[3];
rz(-0.98255052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7716498) q[2];
sx q[2];
rz(-2.6910431) q[2];
sx q[2];
rz(3.0453299) q[2];
rz(-2.7104968) q[3];
sx q[3];
rz(-2.1754706) q[3];
sx q[3];
rz(-2.9935484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8226606) q[0];
sx q[0];
rz(-3.0156101) q[0];
sx q[0];
rz(0.2952964) q[0];
rz(-2.2604306) q[1];
sx q[1];
rz(-2.9144139) q[1];
sx q[1];
rz(0.49240246) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8618362) q[0];
sx q[0];
rz(-1.7375713) q[0];
sx q[0];
rz(1.357973) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6642163) q[2];
sx q[2];
rz(-1.619941) q[2];
sx q[2];
rz(1.3256595) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3141827) q[1];
sx q[1];
rz(-0.35574177) q[1];
sx q[1];
rz(-1.1886503) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6564859) q[3];
sx q[3];
rz(-0.53561775) q[3];
sx q[3];
rz(1.519062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1262576) q[2];
sx q[2];
rz(-0.33053645) q[2];
sx q[2];
rz(0.34273657) q[2];
rz(-2.1526509) q[3];
sx q[3];
rz(-1.9804695) q[3];
sx q[3];
rz(0.69800085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88687503) q[0];
sx q[0];
rz(-0.29061341) q[0];
sx q[0];
rz(1.8732204) q[0];
rz(2.7791924) q[1];
sx q[1];
rz(-0.86530322) q[1];
sx q[1];
rz(1.5579582) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1311296) q[0];
sx q[0];
rz(-2.4005425) q[0];
sx q[0];
rz(-0.69990943) q[0];
x q[1];
rz(-2.1554016) q[2];
sx q[2];
rz(-2.8132943) q[2];
sx q[2];
rz(0.82505783) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4851361) q[1];
sx q[1];
rz(-2.167806) q[1];
sx q[1];
rz(-1.5525425) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2107047) q[3];
sx q[3];
rz(-1.0795169) q[3];
sx q[3];
rz(2.1852428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.6543988) q[2];
sx q[2];
rz(-2.0687658) q[2];
sx q[2];
rz(-2.2844592) q[2];
rz(2.1060139) q[3];
sx q[3];
rz(-2.3643957) q[3];
sx q[3];
rz(0.63151675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.421748) q[0];
sx q[0];
rz(-0.48374614) q[0];
sx q[0];
rz(-2.8068722) q[0];
rz(-0.98182976) q[1];
sx q[1];
rz(-1.5515168) q[1];
sx q[1];
rz(0.88012153) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8557601) q[0];
sx q[0];
rz(-1.8632392) q[0];
sx q[0];
rz(-2.9797735) q[0];
rz(0.17575616) q[2];
sx q[2];
rz(-1.8316016) q[2];
sx q[2];
rz(3.0523022) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0047061027) q[1];
sx q[1];
rz(-1.1019754) q[1];
sx q[1];
rz(-0.15699082) q[1];
x q[2];
rz(-2.2273034) q[3];
sx q[3];
rz(-2.5733272) q[3];
sx q[3];
rz(1.2340581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52241391) q[2];
sx q[2];
rz(-2.652707) q[2];
sx q[2];
rz(1.6467113) q[2];
rz(1.8374247) q[3];
sx q[3];
rz(-2.2924278) q[3];
sx q[3];
rz(-0.18613923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3592767) q[0];
sx q[0];
rz(-0.19141153) q[0];
sx q[0];
rz(-0.070847832) q[0];
rz(-0.62295667) q[1];
sx q[1];
rz(-0.40095913) q[1];
sx q[1];
rz(-2.7080022) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0458023) q[0];
sx q[0];
rz(-1.8754957) q[0];
sx q[0];
rz(1.6030583) q[0];
x q[1];
rz(2.8500227) q[2];
sx q[2];
rz(-1.6256412) q[2];
sx q[2];
rz(2.8758953) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.9476871) q[1];
sx q[1];
rz(-2.8031859) q[1];
sx q[1];
rz(0.29309764) q[1];
x q[2];
rz(-0.44361349) q[3];
sx q[3];
rz(-1.4673691) q[3];
sx q[3];
rz(0.87956968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5101461) q[2];
sx q[2];
rz(-2.6084709) q[2];
sx q[2];
rz(0.637429) q[2];
rz(2.0006477) q[3];
sx q[3];
rz(-0.44064042) q[3];
sx q[3];
rz(2.2250037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7618074) q[0];
sx q[0];
rz(-2.8768235) q[0];
sx q[0];
rz(2.8763212) q[0];
rz(0.028566407) q[1];
sx q[1];
rz(-0.18762372) q[1];
sx q[1];
rz(2.2144337) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87510159) q[0];
sx q[0];
rz(-1.3463839) q[0];
sx q[0];
rz(1.5427179) q[0];
rz(-1.4248542) q[2];
sx q[2];
rz(-1.9371913) q[2];
sx q[2];
rz(-1.4595569) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.68009752) q[1];
sx q[1];
rz(-2.3465183) q[1];
sx q[1];
rz(-2.9976588) q[1];
x q[2];
rz(1.6838491) q[3];
sx q[3];
rz(-1.1498442) q[3];
sx q[3];
rz(-2.497626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.022543) q[2];
sx q[2];
rz(-1.0485342) q[2];
sx q[2];
rz(2.0691464) q[2];
rz(0.33506814) q[3];
sx q[3];
rz(-0.11490331) q[3];
sx q[3];
rz(-2.9160299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91113126) q[0];
sx q[0];
rz(-1.9879531) q[0];
sx q[0];
rz(2.352584) q[0];
rz(-2.543653) q[1];
sx q[1];
rz(-2.0409248) q[1];
sx q[1];
rz(-2.423563) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9791012) q[0];
sx q[0];
rz(-1.4153061) q[0];
sx q[0];
rz(-3.0999712) q[0];
rz(1.687019) q[2];
sx q[2];
rz(-1.5614484) q[2];
sx q[2];
rz(-0.41510328) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8666229) q[1];
sx q[1];
rz(-0.94511813) q[1];
sx q[1];
rz(-0.83557202) q[1];
x q[2];
rz(0.47649033) q[3];
sx q[3];
rz(-0.19987488) q[3];
sx q[3];
rz(1.0651922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.614552) q[2];
sx q[2];
rz(-0.25200945) q[2];
sx q[2];
rz(-1.4177812) q[2];
rz(-2.0333911) q[3];
sx q[3];
rz(-2.7751444) q[3];
sx q[3];
rz(-0.38780701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80477667) q[0];
sx q[0];
rz(-2.5122061) q[0];
sx q[0];
rz(-0.25548536) q[0];
rz(-2.3707223) q[1];
sx q[1];
rz(-1.5893385) q[1];
sx q[1];
rz(-1.593387) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99799081) q[0];
sx q[0];
rz(-1.457347) q[0];
sx q[0];
rz(-1.1313637) q[0];
x q[1];
rz(-0.22996567) q[2];
sx q[2];
rz(-0.88092025) q[2];
sx q[2];
rz(1.752319) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.080049971) q[1];
sx q[1];
rz(-1.953974) q[1];
sx q[1];
rz(-1.8661168) q[1];
rz(-pi) q[2];
rz(0.53867619) q[3];
sx q[3];
rz(-2.4823501) q[3];
sx q[3];
rz(2.0672807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8678681) q[2];
sx q[2];
rz(-0.74213433) q[2];
sx q[2];
rz(-1.254427) q[2];
rz(-0.54002386) q[3];
sx q[3];
rz(-0.050914474) q[3];
sx q[3];
rz(-0.77903265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.9191606) q[0];
sx q[0];
rz(-1.5800911) q[0];
sx q[0];
rz(2.0240361) q[0];
rz(0.22668214) q[1];
sx q[1];
rz(-3.0381028) q[1];
sx q[1];
rz(1.4729952) q[1];
rz(0.37353362) q[2];
sx q[2];
rz(-1.0230156) q[2];
sx q[2];
rz(-0.50470232) q[2];
rz(-0.35861438) q[3];
sx q[3];
rz(-2.3862541) q[3];
sx q[3];
rz(-2.42498) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
