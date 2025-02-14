OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11290057) q[0];
sx q[0];
rz(2.8960462) q[0];
sx q[0];
rz(9.4774376) q[0];
rz(-2.0242937) q[1];
sx q[1];
rz(-0.3408365) q[1];
sx q[1];
rz(-1.0885106) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65050012) q[0];
sx q[0];
rz(-2.740747) q[0];
sx q[0];
rz(-1.7328788) q[0];
rz(-2.229548) q[2];
sx q[2];
rz(-0.17870644) q[2];
sx q[2];
rz(0.9389053) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3237985) q[1];
sx q[1];
rz(-1.4060188) q[1];
sx q[1];
rz(0.64357693) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5616712) q[3];
sx q[3];
rz(-2.7219028) q[3];
sx q[3];
rz(1.4008824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6344305) q[2];
sx q[2];
rz(-1.3619225) q[2];
sx q[2];
rz(2.5334899) q[2];
rz(-2.0173343) q[3];
sx q[3];
rz(-1.2048771) q[3];
sx q[3];
rz(2.2500136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7267777) q[0];
sx q[0];
rz(-1.5730653) q[0];
sx q[0];
rz(2.538105) q[0];
rz(-2.6314349) q[1];
sx q[1];
rz(-0.77168232) q[1];
sx q[1];
rz(1.0268432) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48906836) q[0];
sx q[0];
rz(-1.4155672) q[0];
sx q[0];
rz(0.17431577) q[0];
rz(-pi) q[1];
rz(2.3281664) q[2];
sx q[2];
rz(-1.0283089) q[2];
sx q[2];
rz(2.7493966) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3852147) q[1];
sx q[1];
rz(-2.2526388) q[1];
sx q[1];
rz(2.4663976) q[1];
rz(-pi) q[2];
rz(2.2193416) q[3];
sx q[3];
rz(-0.77510683) q[3];
sx q[3];
rz(-2.4668281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5828731) q[2];
sx q[2];
rz(-1.9931953) q[2];
sx q[2];
rz(2.4951475) q[2];
rz(1.2785814) q[3];
sx q[3];
rz(-2.1918112) q[3];
sx q[3];
rz(-1.4518552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.932514) q[0];
sx q[0];
rz(-0.13020733) q[0];
sx q[0];
rz(1.5721488) q[0];
rz(2.7945844) q[1];
sx q[1];
rz(-1.4748814) q[1];
sx q[1];
rz(-2.2775547) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0542674) q[0];
sx q[0];
rz(-2.1313166) q[0];
sx q[0];
rz(1.3232687) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8225602) q[2];
sx q[2];
rz(-1.1472817) q[2];
sx q[2];
rz(2.8734795) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.62926312) q[1];
sx q[1];
rz(-2.1844668) q[1];
sx q[1];
rz(-2.3113109) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3898464) q[3];
sx q[3];
rz(-2.6084508) q[3];
sx q[3];
rz(-2.0651434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0916834) q[2];
sx q[2];
rz(-2.9361762) q[2];
sx q[2];
rz(-1.4570215) q[2];
rz(-2.125804) q[3];
sx q[3];
rz(-2.6088645) q[3];
sx q[3];
rz(2.8252025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9006573) q[0];
sx q[0];
rz(-2.7641986) q[0];
sx q[0];
rz(-1.578791) q[0];
rz(-2.0511625) q[1];
sx q[1];
rz(-0.54160392) q[1];
sx q[1];
rz(1.9900367) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3131696) q[0];
sx q[0];
rz(-2.1439664) q[0];
sx q[0];
rz(-1.2970807) q[0];
x q[1];
rz(-2.6461965) q[2];
sx q[2];
rz(-2.317111) q[2];
sx q[2];
rz(-1.3935312) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.87428625) q[1];
sx q[1];
rz(-1.6055672) q[1];
sx q[1];
rz(1.0465996) q[1];
x q[2];
rz(-1.5308446) q[3];
sx q[3];
rz(-0.60725799) q[3];
sx q[3];
rz(2.0366397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84732032) q[2];
sx q[2];
rz(-0.5449833) q[2];
sx q[2];
rz(1.6418183) q[2];
rz(-0.13449399) q[3];
sx q[3];
rz(-1.2765063) q[3];
sx q[3];
rz(-0.78181481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0624369) q[0];
sx q[0];
rz(-0.9117313) q[0];
sx q[0];
rz(-0.4656747) q[0];
rz(1.8648719) q[1];
sx q[1];
rz(-1.0036422) q[1];
sx q[1];
rz(-1.0328971) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5860234) q[0];
sx q[0];
rz(-0.15842552) q[0];
sx q[0];
rz(1.6371631) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14012248) q[2];
sx q[2];
rz(-1.6701785) q[2];
sx q[2];
rz(-1.4750208) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3001662) q[1];
sx q[1];
rz(-0.71117095) q[1];
sx q[1];
rz(2.7513642) q[1];
rz(-pi) q[2];
rz(-0.9744076) q[3];
sx q[3];
rz(-0.63035175) q[3];
sx q[3];
rz(0.74385616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5501962) q[2];
sx q[2];
rz(-1.952221) q[2];
sx q[2];
rz(2.2097394) q[2];
rz(-0.094680928) q[3];
sx q[3];
rz(-2.2414312) q[3];
sx q[3];
rz(-1.8315106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4513627) q[0];
sx q[0];
rz(-0.19915038) q[0];
sx q[0];
rz(3.0354101) q[0];
rz(-2.5871318) q[1];
sx q[1];
rz(-0.78840557) q[1];
sx q[1];
rz(-1.6000481) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7035737) q[0];
sx q[0];
rz(-2.7503563) q[0];
sx q[0];
rz(-3.0783669) q[0];
rz(1.1136069) q[2];
sx q[2];
rz(-1.8181674) q[2];
sx q[2];
rz(2.8988225) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20075837) q[1];
sx q[1];
rz(-1.4486099) q[1];
sx q[1];
rz(-2.4442375) q[1];
rz(2.9227436) q[3];
sx q[3];
rz(-0.25222029) q[3];
sx q[3];
rz(0.24009934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9083378) q[2];
sx q[2];
rz(-0.7408064) q[2];
sx q[2];
rz(1.4976658) q[2];
rz(2.6805367) q[3];
sx q[3];
rz(-0.65270972) q[3];
sx q[3];
rz(1.8259995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6719565) q[0];
sx q[0];
rz(-2.3905583) q[0];
sx q[0];
rz(-2.0627956) q[0];
rz(0.59393334) q[1];
sx q[1];
rz(-2.1682231) q[1];
sx q[1];
rz(0.22720164) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0086454) q[0];
sx q[0];
rz(-1.4221749) q[0];
sx q[0];
rz(3.0287543) q[0];
rz(-pi) q[1];
x q[1];
rz(1.237713) q[2];
sx q[2];
rz(-1.7041612) q[2];
sx q[2];
rz(-1.9072717) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6727184) q[1];
sx q[1];
rz(-2.1203845) q[1];
sx q[1];
rz(-2.7764715) q[1];
x q[2];
rz(-1.8232089) q[3];
sx q[3];
rz(-0.38312437) q[3];
sx q[3];
rz(0.33465222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6195153) q[2];
sx q[2];
rz(-0.41831133) q[2];
sx q[2];
rz(2.412001) q[2];
rz(-1.4531762) q[3];
sx q[3];
rz(-1.4540693) q[3];
sx q[3];
rz(2.5280473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59298092) q[0];
sx q[0];
rz(-2.225003) q[0];
sx q[0];
rz(-2.7899637) q[0];
rz(0.31373203) q[1];
sx q[1];
rz(-1.2064563) q[1];
sx q[1];
rz(1.7650013) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54395318) q[0];
sx q[0];
rz(-0.40297976) q[0];
sx q[0];
rz(-3.0209994) q[0];
rz(1.4046762) q[2];
sx q[2];
rz(-1.4567647) q[2];
sx q[2];
rz(1.5064552) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.87159705) q[1];
sx q[1];
rz(-2.13378) q[1];
sx q[1];
rz(1.7558805) q[1];
rz(-pi) q[2];
rz(-1.8217161) q[3];
sx q[3];
rz(-0.38057113) q[3];
sx q[3];
rz(-1.5151092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54479638) q[2];
sx q[2];
rz(-0.72303855) q[2];
sx q[2];
rz(1.5210305) q[2];
rz(0.14791402) q[3];
sx q[3];
rz(-1.3444129) q[3];
sx q[3];
rz(1.773905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3677597) q[0];
sx q[0];
rz(-1.2667043) q[0];
sx q[0];
rz(1.2441147) q[0];
rz(-1.7077839) q[1];
sx q[1];
rz(-1.5807187) q[1];
sx q[1];
rz(-0.22020766) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7871683) q[0];
sx q[0];
rz(-2.1666013) q[0];
sx q[0];
rz(-2.6902249) q[0];
x q[1];
rz(1.2559863) q[2];
sx q[2];
rz(-2.6546835) q[2];
sx q[2];
rz(-1.032342) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9088194) q[1];
sx q[1];
rz(-2.4202883) q[1];
sx q[1];
rz(2.2386293) q[1];
x q[2];
rz(3.0433398) q[3];
sx q[3];
rz(-1.6643401) q[3];
sx q[3];
rz(-2.6152809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1844909) q[2];
sx q[2];
rz(-1.1471006) q[2];
sx q[2];
rz(1.5000878) q[2];
rz(-0.084065048) q[3];
sx q[3];
rz(-1.2596687) q[3];
sx q[3];
rz(-3.0715004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67749196) q[0];
sx q[0];
rz(-0.75834948) q[0];
sx q[0];
rz(-2.7229934) q[0];
rz(-0.94340008) q[1];
sx q[1];
rz(-1.6523596) q[1];
sx q[1];
rz(-2.8387866) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3269791) q[0];
sx q[0];
rz(-1.7075065) q[0];
sx q[0];
rz(-1.4025406) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74818989) q[2];
sx q[2];
rz(-2.8837969) q[2];
sx q[2];
rz(-0.92743826) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.329656) q[1];
sx q[1];
rz(-1.6534263) q[1];
sx q[1];
rz(-1.5272115) q[1];
x q[2];
rz(1.8859768) q[3];
sx q[3];
rz(-2.7806296) q[3];
sx q[3];
rz(2.7415467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6454978) q[2];
sx q[2];
rz(-2.2690513) q[2];
sx q[2];
rz(-2.3910451) q[2];
rz(-1.4613072) q[3];
sx q[3];
rz(-0.76992005) q[3];
sx q[3];
rz(0.51699483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1207598) q[0];
sx q[0];
rz(-1.5789541) q[0];
sx q[0];
rz(-1.5776237) q[0];
rz(-1.9137406) q[1];
sx q[1];
rz(-2.6422983) q[1];
sx q[1];
rz(-1.9849389) q[1];
rz(1.7534134) q[2];
sx q[2];
rz(-1.0116521) q[2];
sx q[2];
rz(0.73965363) q[2];
rz(2.078656) q[3];
sx q[3];
rz(-0.36527562) q[3];
sx q[3];
rz(1.9648413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
