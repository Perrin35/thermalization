OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9031653) q[0];
sx q[0];
rz(2.981346) q[0];
sx q[0];
rz(7.5709406) q[0];
rz(2.2924478) q[1];
sx q[1];
rz(-1.949911) q[1];
sx q[1];
rz(-0.75872672) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0512801) q[0];
sx q[0];
rz(-2.7669124) q[0];
sx q[0];
rz(-0.6775585) q[0];
rz(-pi) q[1];
rz(2.0082194) q[2];
sx q[2];
rz(-0.52825275) q[2];
sx q[2];
rz(-2.4932441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78658653) q[1];
sx q[1];
rz(-2.2753582) q[1];
sx q[1];
rz(2.658872) q[1];
rz(2.3146199) q[3];
sx q[3];
rz(-2.0872413) q[3];
sx q[3];
rz(1.4475105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3005001) q[2];
sx q[2];
rz(-0.79885834) q[2];
sx q[2];
rz(-0.50660261) q[2];
rz(1.5233585) q[3];
sx q[3];
rz(-2.5169499) q[3];
sx q[3];
rz(-0.055559572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27187207) q[0];
sx q[0];
rz(-2.6756918) q[0];
sx q[0];
rz(3.077935) q[0];
rz(-1.9438538) q[1];
sx q[1];
rz(-1.627219) q[1];
sx q[1];
rz(1.0187842) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2301529) q[0];
sx q[0];
rz(-0.24556118) q[0];
sx q[0];
rz(2.5018238) q[0];
rz(-pi) q[1];
rz(0.59661023) q[2];
sx q[2];
rz(-0.40477529) q[2];
sx q[2];
rz(-0.34133729) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9460822) q[1];
sx q[1];
rz(-1.2942593) q[1];
sx q[1];
rz(-2.4184188) q[1];
rz(-0.13354729) q[3];
sx q[3];
rz(-2.3192647) q[3];
sx q[3];
rz(0.68822569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4360875) q[2];
sx q[2];
rz(-0.59385308) q[2];
sx q[2];
rz(0.84211055) q[2];
rz(2.0841133) q[3];
sx q[3];
rz(-0.92856854) q[3];
sx q[3];
rz(1.9717982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24404003) q[0];
sx q[0];
rz(-2.121448) q[0];
sx q[0];
rz(-1.1919588) q[0];
rz(2.2167218) q[1];
sx q[1];
rz(-2.4332739) q[1];
sx q[1];
rz(-1.301773) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.155543) q[0];
sx q[0];
rz(-1.0240004) q[0];
sx q[0];
rz(2.5803671) q[0];
rz(0.19916735) q[2];
sx q[2];
rz(-1.609243) q[2];
sx q[2];
rz(-3.0035915) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.30693248) q[1];
sx q[1];
rz(-1.3921157) q[1];
sx q[1];
rz(-0.4731005) q[1];
x q[2];
rz(0.25721217) q[3];
sx q[3];
rz(-1.9055771) q[3];
sx q[3];
rz(1.0509261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7630345) q[2];
sx q[2];
rz(-0.083746567) q[2];
sx q[2];
rz(0.087510022) q[2];
rz(-1.3148426) q[3];
sx q[3];
rz(-1.0564691) q[3];
sx q[3];
rz(2.1480613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73115504) q[0];
sx q[0];
rz(-1.1719828) q[0];
sx q[0];
rz(-0.80899578) q[0];
rz(-1.0357098) q[1];
sx q[1];
rz(-2.6782942) q[1];
sx q[1];
rz(-1.0332003) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3147723) q[0];
sx q[0];
rz(-0.90703947) q[0];
sx q[0];
rz(-3.0919667) q[0];
rz(-pi) q[1];
rz(2.7055351) q[2];
sx q[2];
rz(-1.0761217) q[2];
sx q[2];
rz(0.058017284) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82419187) q[1];
sx q[1];
rz(-0.5017952) q[1];
sx q[1];
rz(-0.3806033) q[1];
rz(2.0508893) q[3];
sx q[3];
rz(-2.0629145) q[3];
sx q[3];
rz(-1.1034213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0741299) q[2];
sx q[2];
rz(-2.9758657) q[2];
sx q[2];
rz(-1.6039675) q[2];
rz(-1.4175203) q[3];
sx q[3];
rz(-2.0310903) q[3];
sx q[3];
rz(-2.0541151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23987016) q[0];
sx q[0];
rz(-2.9481695) q[0];
sx q[0];
rz(1.9523917) q[0];
rz(2.4173648) q[1];
sx q[1];
rz(-1.291357) q[1];
sx q[1];
rz(-1.474818) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9928707) q[0];
sx q[0];
rz(-1.1445415) q[0];
sx q[0];
rz(-0.54963407) q[0];
rz(-pi) q[1];
rz(0.15691136) q[2];
sx q[2];
rz(-2.3381066) q[2];
sx q[2];
rz(0.29565865) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1777623) q[1];
sx q[1];
rz(-2.6796067) q[1];
sx q[1];
rz(0.83779184) q[1];
rz(0.46012043) q[3];
sx q[3];
rz(-1.1183721) q[3];
sx q[3];
rz(-2.2442073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7597947) q[2];
sx q[2];
rz(-0.90961027) q[2];
sx q[2];
rz(-0.46662113) q[2];
rz(-1.7032547) q[3];
sx q[3];
rz(-1.2983026) q[3];
sx q[3];
rz(2.2013825) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4284215) q[0];
sx q[0];
rz(-0.82601014) q[0];
sx q[0];
rz(0.010350479) q[0];
rz(-1.5230644) q[1];
sx q[1];
rz(-1.7489988) q[1];
sx q[1];
rz(1.3759605) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9860366) q[0];
sx q[0];
rz(-1.9080093) q[0];
sx q[0];
rz(-2.8854779) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59734663) q[2];
sx q[2];
rz(-0.22946363) q[2];
sx q[2];
rz(-2.2962183) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.79094564) q[1];
sx q[1];
rz(-0.53574359) q[1];
sx q[1];
rz(-2.78113) q[1];
x q[2];
rz(0.74919219) q[3];
sx q[3];
rz(-1.0674879) q[3];
sx q[3];
rz(2.1319413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7213664) q[2];
sx q[2];
rz(-1.0488291) q[2];
sx q[2];
rz(2.4143207) q[2];
rz(2.6266802) q[3];
sx q[3];
rz(-1.8102831) q[3];
sx q[3];
rz(1.7236727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10776831) q[0];
sx q[0];
rz(-1.6078147) q[0];
sx q[0];
rz(0.4069826) q[0];
rz(2.175323) q[1];
sx q[1];
rz(-0.85710183) q[1];
sx q[1];
rz(-3.0028717) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9753067) q[0];
sx q[0];
rz(-1.606267) q[0];
sx q[0];
rz(-1.6871638) q[0];
rz(-pi) q[1];
rz(0.66731668) q[2];
sx q[2];
rz(-2.7936068) q[2];
sx q[2];
rz(1.1225357) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1410576) q[1];
sx q[1];
rz(-2.6319017) q[1];
sx q[1];
rz(0.7668709) q[1];
rz(2.9937708) q[3];
sx q[3];
rz(-1.9114219) q[3];
sx q[3];
rz(-0.56440777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2762642) q[2];
sx q[2];
rz(-1.7636969) q[2];
sx q[2];
rz(2.830937) q[2];
rz(1.3079414) q[3];
sx q[3];
rz(-1.6303948) q[3];
sx q[3];
rz(-2.7697897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0528316) q[0];
sx q[0];
rz(-2.796687) q[0];
sx q[0];
rz(3.0533691) q[0];
rz(0.010559646) q[1];
sx q[1];
rz(-1.6190395) q[1];
sx q[1];
rz(-2.6710076) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4751401) q[0];
sx q[0];
rz(-1.0658403) q[0];
sx q[0];
rz(0.83374896) q[0];
rz(-pi) q[1];
rz(-1.3731277) q[2];
sx q[2];
rz(-2.3680147) q[2];
sx q[2];
rz(-2.8171223) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1697929) q[1];
sx q[1];
rz(-1.3108675) q[1];
sx q[1];
rz(-3.0539576) q[1];
rz(-pi) q[2];
rz(-0.13981522) q[3];
sx q[3];
rz(-1.9203517) q[3];
sx q[3];
rz(2.4619554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.18870246) q[2];
sx q[2];
rz(-2.1850093) q[2];
sx q[2];
rz(-1.0799705) q[2];
rz(0.88932577) q[3];
sx q[3];
rz(-1.5417128) q[3];
sx q[3];
rz(1.5573474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
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
rz(0.93961632) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25519263) q[0];
sx q[0];
rz(-1.5879244) q[0];
sx q[0];
rz(-3.1411132) q[0];
rz(-pi) q[1];
rz(2.7857271) q[2];
sx q[2];
rz(-1.8345222) q[2];
sx q[2];
rz(0.10019324) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6691362) q[1];
sx q[1];
rz(-1.3208773) q[1];
sx q[1];
rz(0.12052287) q[1];
rz(-pi) q[2];
rz(1.1577179) q[3];
sx q[3];
rz(-1.7677757) q[3];
sx q[3];
rz(-2.911681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8115936) q[2];
sx q[2];
rz(-2.1758695) q[2];
sx q[2];
rz(2.7962371) q[2];
rz(-1.6164448) q[3];
sx q[3];
rz(-1.7914146) q[3];
sx q[3];
rz(2.7143872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3578607) q[0];
sx q[0];
rz(-1.9851728) q[0];
sx q[0];
rz(0.082948908) q[0];
rz(0.16231617) q[1];
sx q[1];
rz(-1.6971308) q[1];
sx q[1];
rz(-0.69581318) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2029977) q[0];
sx q[0];
rz(-1.2913398) q[0];
sx q[0];
rz(-1.6142815) q[0];
x q[1];
rz(-1.7000574) q[2];
sx q[2];
rz(-1.5737572) q[2];
sx q[2];
rz(-1.9287314) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2999794) q[1];
sx q[1];
rz(-2.0149954) q[1];
sx q[1];
rz(-1.5496553) q[1];
rz(0.70503321) q[3];
sx q[3];
rz(-1.789973) q[3];
sx q[3];
rz(-0.33112835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7210228) q[2];
sx q[2];
rz(-2.9837954) q[2];
sx q[2];
rz(-1.2201307) q[2];
rz(-2.5681791) q[3];
sx q[3];
rz(-1.8107332) q[3];
sx q[3];
rz(-0.24513182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2777916) q[0];
sx q[0];
rz(-1.5024804) q[0];
sx q[0];
rz(-0.13308751) q[0];
rz(-0.0028903891) q[1];
sx q[1];
rz(-2.7255701) q[1];
sx q[1];
rz(2.4400673) q[1];
rz(-1.8194503) q[2];
sx q[2];
rz(-1.7380309) q[2];
sx q[2];
rz(-2.0792014) q[2];
rz(-1.5118128) q[3];
sx q[3];
rz(-2.3702757) q[3];
sx q[3];
rz(2.6487917) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
