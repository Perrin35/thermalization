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
rz(-1.4200014) q[0];
sx q[0];
rz(-0.95945946) q[0];
sx q[0];
rz(2.1009011) q[0];
rz(-0.33848441) q[1];
sx q[1];
rz(-1.4893293) q[1];
sx q[1];
rz(-1.4298061) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4066577) q[0];
sx q[0];
rz(-1.2996607) q[0];
sx q[0];
rz(1.6938637) q[0];
rz(1.0770887) q[2];
sx q[2];
rz(-0.66254751) q[2];
sx q[2];
rz(0.27499378) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0104645) q[1];
sx q[1];
rz(-0.87112037) q[1];
sx q[1];
rz(-1.9424428) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40994196) q[3];
sx q[3];
rz(-2.5363208) q[3];
sx q[3];
rz(-2.2469478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0040943) q[2];
sx q[2];
rz(-1.7137004) q[2];
sx q[2];
rz(-3.0924228) q[2];
rz(-0.54685012) q[3];
sx q[3];
rz(-0.85421908) q[3];
sx q[3];
rz(1.2949519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26286724) q[0];
sx q[0];
rz(-0.85064369) q[0];
sx q[0];
rz(-0.039462939) q[0];
rz(2.4354758) q[1];
sx q[1];
rz(-1.9673653) q[1];
sx q[1];
rz(1.8810693) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17968125) q[0];
sx q[0];
rz(-0.42191178) q[0];
sx q[0];
rz(1.4279813) q[0];
rz(-pi) q[1];
rz(-2.1095545) q[2];
sx q[2];
rz(-2.1646087) q[2];
sx q[2];
rz(-2.3597882) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.83585868) q[1];
sx q[1];
rz(-1.6129588) q[1];
sx q[1];
rz(-0.034138676) q[1];
x q[2];
rz(0.32716635) q[3];
sx q[3];
rz(-2.9370902) q[3];
sx q[3];
rz(0.12302264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4827106) q[2];
sx q[2];
rz(-1.005268) q[2];
sx q[2];
rz(0.91747326) q[2];
rz(2.7190123) q[3];
sx q[3];
rz(-1.249142) q[3];
sx q[3];
rz(-1.0379855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68040401) q[0];
sx q[0];
rz(-0.55484158) q[0];
sx q[0];
rz(-2.19221) q[0];
rz(-1.2481015) q[1];
sx q[1];
rz(-0.58369842) q[1];
sx q[1];
rz(-1.6277574) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7010403) q[0];
sx q[0];
rz(-1.8157542) q[0];
sx q[0];
rz(-2.7666758) q[0];
x q[1];
rz(-1.1749985) q[2];
sx q[2];
rz(-0.51766073) q[2];
sx q[2];
rz(0.84830059) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6520034) q[1];
sx q[1];
rz(-1.4966652) q[1];
sx q[1];
rz(2.8531122) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0120609) q[3];
sx q[3];
rz(-1.2848275) q[3];
sx q[3];
rz(0.27559552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7013487) q[2];
sx q[2];
rz(-1.1260208) q[2];
sx q[2];
rz(1.6181642) q[2];
rz(-2.8273888) q[3];
sx q[3];
rz(-2.115963) q[3];
sx q[3];
rz(-1.5020802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0358129) q[0];
sx q[0];
rz(-2.950225) q[0];
sx q[0];
rz(0.15709239) q[0];
rz(0.13520232) q[1];
sx q[1];
rz(-2.356485) q[1];
sx q[1];
rz(0.33321998) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5557809) q[0];
sx q[0];
rz(-0.69792047) q[0];
sx q[0];
rz(-0.74549739) q[0];
rz(-pi) q[1];
rz(2.2991247) q[2];
sx q[2];
rz(-1.7211452) q[2];
sx q[2];
rz(-1.6347234) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4823158) q[1];
sx q[1];
rz(-1.2219861) q[1];
sx q[1];
rz(-2.8608843) q[1];
rz(-0.85835056) q[3];
sx q[3];
rz(-1.3305404) q[3];
sx q[3];
rz(2.1879856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.80804431) q[2];
sx q[2];
rz(-2.4925888) q[2];
sx q[2];
rz(-1.4487079) q[2];
rz(0.73271218) q[3];
sx q[3];
rz(-1.9219739) q[3];
sx q[3];
rz(-1.5318416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6919959) q[0];
sx q[0];
rz(-1.4866328) q[0];
sx q[0];
rz(0.23155552) q[0];
rz(-2.4202276) q[1];
sx q[1];
rz(-2.2504579) q[1];
sx q[1];
rz(-2.2931633) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4554286) q[0];
sx q[0];
rz(-1.2368742) q[0];
sx q[0];
rz(2.0823821) q[0];
rz(-0.5428585) q[2];
sx q[2];
rz(-1.9221537) q[2];
sx q[2];
rz(1.5962102) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2347772) q[1];
sx q[1];
rz(-2.7473852) q[1];
sx q[1];
rz(-0.18313198) q[1];
rz(2.7293398) q[3];
sx q[3];
rz(-2.8839211) q[3];
sx q[3];
rz(-2.4564956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1344177) q[2];
sx q[2];
rz(-2.538372) q[2];
sx q[2];
rz(-1.940654) q[2];
rz(0.68931556) q[3];
sx q[3];
rz(-1.9715693) q[3];
sx q[3];
rz(-0.56572604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1597964) q[0];
sx q[0];
rz(-1.6098276) q[0];
sx q[0];
rz(1.7933886) q[0];
rz(-2.4234405) q[1];
sx q[1];
rz(-2.5703057) q[1];
sx q[1];
rz(1.8498373) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6245218) q[0];
sx q[0];
rz(-2.7291131) q[0];
sx q[0];
rz(2.5530091) q[0];
rz(-pi) q[1];
rz(0.97219418) q[2];
sx q[2];
rz(-1.2568398) q[2];
sx q[2];
rz(1.0034059) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3499432) q[1];
sx q[1];
rz(-1.4022489) q[1];
sx q[1];
rz(-1.4838292) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.060543493) q[3];
sx q[3];
rz(-2.745564) q[3];
sx q[3];
rz(-0.76914364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1527839) q[2];
sx q[2];
rz(-0.96502105) q[2];
sx q[2];
rz(-2.7394845) q[2];
rz(-2.1357338) q[3];
sx q[3];
rz(-0.22336762) q[3];
sx q[3];
rz(1.7042101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.724029) q[0];
sx q[0];
rz(-0.14823866) q[0];
sx q[0];
rz(-1.3295901) q[0];
rz(-1.3279462) q[1];
sx q[1];
rz(-1.4168394) q[1];
sx q[1];
rz(-2.6899469) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76021323) q[0];
sx q[0];
rz(-0.45256786) q[0];
sx q[0];
rz(0.20968881) q[0];
x q[1];
rz(0.30173413) q[2];
sx q[2];
rz(-0.44281755) q[2];
sx q[2];
rz(-0.91609611) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.627915) q[1];
sx q[1];
rz(-1.5202402) q[1];
sx q[1];
rz(-1.7611658) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90410952) q[3];
sx q[3];
rz(-1.0268758) q[3];
sx q[3];
rz(-2.0407807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9391276) q[2];
sx q[2];
rz(-1.3718601) q[2];
sx q[2];
rz(0.28787127) q[2];
rz(-1.4816083) q[3];
sx q[3];
rz(-2.2899254) q[3];
sx q[3];
rz(1.6847346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57956368) q[0];
sx q[0];
rz(-2.1810739) q[0];
sx q[0];
rz(0.55291837) q[0];
rz(1.8909854) q[1];
sx q[1];
rz(-2.7313373) q[1];
sx q[1];
rz(-1.3804573) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6621679) q[0];
sx q[0];
rz(-1.933455) q[0];
sx q[0];
rz(0.5525241) q[0];
rz(0.12166656) q[2];
sx q[2];
rz(-0.98699283) q[2];
sx q[2];
rz(2.9505299) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1637818) q[1];
sx q[1];
rz(-1.8076784) q[1];
sx q[1];
rz(1.6131748) q[1];
x q[2];
rz(2.9743324) q[3];
sx q[3];
rz(-0.80771192) q[3];
sx q[3];
rz(-0.72116646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.55014253) q[2];
sx q[2];
rz(-1.3478792) q[2];
sx q[2];
rz(2.4264917) q[2];
rz(3.0969369) q[3];
sx q[3];
rz(-0.5831334) q[3];
sx q[3];
rz(-2.513212) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1735246) q[0];
sx q[0];
rz(-0.76924291) q[0];
sx q[0];
rz(2.5007057) q[0];
rz(-1.6454654) q[1];
sx q[1];
rz(-2.75664) q[1];
sx q[1];
rz(-2.4900751) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7018873) q[0];
sx q[0];
rz(-0.7535075) q[0];
sx q[0];
rz(1.9508595) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6518238) q[2];
sx q[2];
rz(-1.3636053) q[2];
sx q[2];
rz(-1.734601) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2711598) q[1];
sx q[1];
rz(-1.203013) q[1];
sx q[1];
rz(0.8720542) q[1];
rz(-pi) q[2];
rz(-2.9504029) q[3];
sx q[3];
rz(-1.9131197) q[3];
sx q[3];
rz(-2.1514307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4064101) q[2];
sx q[2];
rz(-2.0969022) q[2];
sx q[2];
rz(-2.6835105) q[2];
rz(-0.011890751) q[3];
sx q[3];
rz(-1.0606822) q[3];
sx q[3];
rz(0.021765821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7387725) q[0];
sx q[0];
rz(-2.9467376) q[0];
sx q[0];
rz(-0.028976945) q[0];
rz(2.2150691) q[1];
sx q[1];
rz(-1.69311) q[1];
sx q[1];
rz(-0.40625939) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0615912) q[0];
sx q[0];
rz(-1.5621981) q[0];
sx q[0];
rz(-2.8614869) q[0];
x q[1];
rz(-2.788061) q[2];
sx q[2];
rz(-0.95530704) q[2];
sx q[2];
rz(-0.36223499) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.528842) q[1];
sx q[1];
rz(-1.6019197) q[1];
sx q[1];
rz(-1.2791388) q[1];
x q[2];
rz(0.80244949) q[3];
sx q[3];
rz(-0.93367773) q[3];
sx q[3];
rz(-0.040590057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3709092) q[2];
sx q[2];
rz(-1.4981046) q[2];
sx q[2];
rz(-0.59477273) q[2];
rz(0.70991436) q[3];
sx q[3];
rz(-0.9570595) q[3];
sx q[3];
rz(-0.94533581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57472537) q[0];
sx q[0];
rz(-2.18676) q[0];
sx q[0];
rz(-2.2796897) q[0];
rz(0.30543874) q[1];
sx q[1];
rz(-1.8157235) q[1];
sx q[1];
rz(-2.9095412) q[1];
rz(2.7965056) q[2];
sx q[2];
rz(-1.9660334) q[2];
sx q[2];
rz(2.6031969) q[2];
rz(-1.7954682) q[3];
sx q[3];
rz(-1.9763038) q[3];
sx q[3];
rz(1.4122813) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
