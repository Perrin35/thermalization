OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.5716612) q[0];
sx q[0];
rz(-0.067458955) q[0];
sx q[0];
rz(10.098739) q[0];
rz(-0.31710467) q[1];
sx q[1];
rz(-1.6333406) q[1];
sx q[1];
rz(-0.79467264) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7793517) q[0];
sx q[0];
rz(-1.2138214) q[0];
sx q[0];
rz(-0.39874052) q[0];
rz(-2.652466) q[2];
sx q[2];
rz(-1.7009652) q[2];
sx q[2];
rz(-2.7444073) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0104048) q[1];
sx q[1];
rz(-0.63001761) q[1];
sx q[1];
rz(2.5351934) q[1];
rz(-pi) q[2];
rz(3.0029293) q[3];
sx q[3];
rz(-2.0141467) q[3];
sx q[3];
rz(-1.5594304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6661466) q[2];
sx q[2];
rz(-2.011516) q[2];
sx q[2];
rz(1.8165992) q[2];
rz(-2.825286) q[3];
sx q[3];
rz(-0.91747147) q[3];
sx q[3];
rz(2.2345208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1871724) q[0];
sx q[0];
rz(-0.81048727) q[0];
sx q[0];
rz(2.1818838) q[0];
rz(1.6540487) q[1];
sx q[1];
rz(-0.68351775) q[1];
sx q[1];
rz(-2.5141292) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7213631) q[0];
sx q[0];
rz(-0.010275928) q[0];
sx q[0];
rz(2.4179439) q[0];
x q[1];
rz(-2.9299514) q[2];
sx q[2];
rz(-1.5492808) q[2];
sx q[2];
rz(0.053552901) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5761766) q[1];
sx q[1];
rz(-1.4285354) q[1];
sx q[1];
rz(-1.8131282) q[1];
rz(1.1153568) q[3];
sx q[3];
rz(-3.1133828) q[3];
sx q[3];
rz(1.3970323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.54005694) q[2];
sx q[2];
rz(-1.4011708) q[2];
sx q[2];
rz(2.1616948) q[2];
rz(1.4510138) q[3];
sx q[3];
rz(-0.7106978) q[3];
sx q[3];
rz(1.6987945) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5805805) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(1.3797492) q[0];
rz(2.9648932) q[1];
sx q[1];
rz(-1.2172164) q[1];
sx q[1];
rz(-1.3476936) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4189258) q[0];
sx q[0];
rz(-1.8113266) q[0];
sx q[0];
rz(-2.6618883) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9413207) q[2];
sx q[2];
rz(-1.592604) q[2];
sx q[2];
rz(-2.3238316) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2587535) q[1];
sx q[1];
rz(-1.6707318) q[1];
sx q[1];
rz(3.019481) q[1];
rz(-1.3860116) q[3];
sx q[3];
rz(-1.4818958) q[3];
sx q[3];
rz(2.256306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1220864) q[2];
sx q[2];
rz(-2.2821189) q[2];
sx q[2];
rz(1.1536095) q[2];
rz(2.7518318) q[3];
sx q[3];
rz(-0.4699769) q[3];
sx q[3];
rz(0.78891689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47238123) q[0];
sx q[0];
rz(-0.62494576) q[0];
sx q[0];
rz(-2.5770082) q[0];
rz(-3.0584884) q[1];
sx q[1];
rz(-1.9080947) q[1];
sx q[1];
rz(-2.9262503) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1911083) q[0];
sx q[0];
rz(-1.2380301) q[0];
sx q[0];
rz(2.3118408) q[0];
rz(-0.16673659) q[2];
sx q[2];
rz(-2.1519289) q[2];
sx q[2];
rz(2.2819448) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6268057) q[1];
sx q[1];
rz(-1.6759911) q[1];
sx q[1];
rz(-1.766596) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1057042) q[3];
sx q[3];
rz(-2.0715189) q[3];
sx q[3];
rz(-0.53946686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5061491) q[2];
sx q[2];
rz(-2.052867) q[2];
sx q[2];
rz(-2.5396458) q[2];
rz(2.4510032) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(2.5155892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2895806) q[0];
sx q[0];
rz(-0.99442011) q[0];
sx q[0];
rz(-0.96486282) q[0];
rz(-3.124974) q[1];
sx q[1];
rz(-0.42533541) q[1];
sx q[1];
rz(-0.66666493) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0913038) q[0];
sx q[0];
rz(-1.4698403) q[0];
sx q[0];
rz(1.9712649) q[0];
rz(-0.35258106) q[2];
sx q[2];
rz(-2.0902299) q[2];
sx q[2];
rz(2.7378766) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1866245) q[1];
sx q[1];
rz(-0.72997) q[1];
sx q[1];
rz(-2.7093922) q[1];
rz(-0.93249647) q[3];
sx q[3];
rz(-2.7780048) q[3];
sx q[3];
rz(-2.5225085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2255286) q[2];
sx q[2];
rz(-1.3801882) q[2];
sx q[2];
rz(1.8886245) q[2];
rz(-1.8270252) q[3];
sx q[3];
rz(-1.9606749) q[3];
sx q[3];
rz(1.130828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0261633) q[0];
sx q[0];
rz(-2.207343) q[0];
sx q[0];
rz(-1.3622267) q[0];
rz(1.7419787) q[1];
sx q[1];
rz(-1.5618961) q[1];
sx q[1];
rz(1.5302352) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8996457) q[0];
sx q[0];
rz(-2.7149704) q[0];
sx q[0];
rz(-2.9645779) q[0];
rz(-pi) q[1];
rz(2.7026664) q[2];
sx q[2];
rz(-1.8669087) q[2];
sx q[2];
rz(0.93026464) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1215026) q[1];
sx q[1];
rz(-1.1768747) q[1];
sx q[1];
rz(2.9083088) q[1];
rz(0.050614428) q[3];
sx q[3];
rz(-1.8236056) q[3];
sx q[3];
rz(-1.8265754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.35649148) q[2];
sx q[2];
rz(-2.1337528) q[2];
sx q[2];
rz(-0.029416857) q[2];
rz(-1.4908837) q[3];
sx q[3];
rz(-2.1865032) q[3];
sx q[3];
rz(1.9194549) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4239663) q[0];
sx q[0];
rz(-0.76537776) q[0];
sx q[0];
rz(0.18280612) q[0];
rz(2.706066) q[1];
sx q[1];
rz(-0.85414779) q[1];
sx q[1];
rz(-0.51087728) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8683838) q[0];
sx q[0];
rz(-1.3084992) q[0];
sx q[0];
rz(-1.457085) q[0];
rz(-pi) q[1];
rz(0.81803825) q[2];
sx q[2];
rz(-2.15089) q[2];
sx q[2];
rz(-1.8028508) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.70799815) q[1];
sx q[1];
rz(-2.3494548) q[1];
sx q[1];
rz(0.57452332) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0914518) q[3];
sx q[3];
rz(-0.58590349) q[3];
sx q[3];
rz(1.0019765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2304948) q[2];
sx q[2];
rz(-0.7926597) q[2];
sx q[2];
rz(2.5210209) q[2];
rz(2.7190322) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(1.2703936) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7997416) q[0];
sx q[0];
rz(-0.36062476) q[0];
sx q[0];
rz(2.0525232) q[0];
rz(2.0607121) q[1];
sx q[1];
rz(-1.8573449) q[1];
sx q[1];
rz(-0.63527766) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23112049) q[0];
sx q[0];
rz(-1.4164682) q[0];
sx q[0];
rz(-1.2088103) q[0];
x q[1];
rz(-2.9674171) q[2];
sx q[2];
rz(-2.0617699) q[2];
sx q[2];
rz(0.27506405) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0465225) q[1];
sx q[1];
rz(-0.97192837) q[1];
sx q[1];
rz(1.2477161) q[1];
x q[2];
rz(-2.1463296) q[3];
sx q[3];
rz(-1.4008153) q[3];
sx q[3];
rz(0.02804027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.13715956) q[2];
sx q[2];
rz(-1.4617504) q[2];
sx q[2];
rz(3.1414202) q[2];
rz(0.6862644) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(2.2837158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8453318) q[0];
sx q[0];
rz(-2.2289472) q[0];
sx q[0];
rz(2.4849179) q[0];
rz(-2.7776921) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(2.231853) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1014935) q[0];
sx q[0];
rz(-1.3001406) q[0];
sx q[0];
rz(0.33466848) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7889195) q[2];
sx q[2];
rz(-1.2303196) q[2];
sx q[2];
rz(-3.043963) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7418356) q[1];
sx q[1];
rz(-1.4474807) q[1];
sx q[1];
rz(2.9305305) q[1];
rz(-pi) q[2];
rz(0.25760381) q[3];
sx q[3];
rz(-1.2649049) q[3];
sx q[3];
rz(0.56759873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.83071128) q[2];
sx q[2];
rz(-0.41822663) q[2];
sx q[2];
rz(-0.48661423) q[2];
rz(-0.1085554) q[3];
sx q[3];
rz(-0.71043772) q[3];
sx q[3];
rz(-1.9177115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20031032) q[0];
sx q[0];
rz(-2.0599984) q[0];
sx q[0];
rz(-2.6249028) q[0];
rz(-1.5962881) q[1];
sx q[1];
rz(-2.089112) q[1];
sx q[1];
rz(-0.89458481) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45458083) q[0];
sx q[0];
rz(-2.5068388) q[0];
sx q[0];
rz(-1.4995585) q[0];
rz(-pi) q[1];
rz(-1.1012494) q[2];
sx q[2];
rz(-1.38648) q[2];
sx q[2];
rz(-1.4354401) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4641061) q[1];
sx q[1];
rz(-1.9312526) q[1];
sx q[1];
rz(0.067671138) q[1];
x q[2];
rz(1.3668725) q[3];
sx q[3];
rz(-0.66019928) q[3];
sx q[3];
rz(3.077293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9937925) q[2];
sx q[2];
rz(-0.96269572) q[2];
sx q[2];
rz(1.0220035) q[2];
rz(-1.3379124) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(2.5125304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.548303) q[0];
sx q[0];
rz(-1.6175445) q[0];
sx q[0];
rz(2.2486726) q[0];
rz(-1.0409566) q[1];
sx q[1];
rz(-3.0491842) q[1];
sx q[1];
rz(-1.4728117) q[1];
rz(-2.2330877) q[2];
sx q[2];
rz(-1.1237332) q[2];
sx q[2];
rz(-0.087469812) q[2];
rz(-1.7265504) q[3];
sx q[3];
rz(-0.30434904) q[3];
sx q[3];
rz(0.54868922) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
