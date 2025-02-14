OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.469874) q[0];
sx q[0];
rz(-1.9865541) q[0];
sx q[0];
rz(1.7116829) q[0];
rz(-2.6861796) q[1];
sx q[1];
rz(-2.5502584) q[1];
sx q[1];
rz(1.1270181) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9874632) q[0];
sx q[0];
rz(-2.9107981) q[0];
sx q[0];
rz(1.5202281) q[0];
x q[1];
rz(2.7989109) q[2];
sx q[2];
rz(-1.8215873) q[2];
sx q[2];
rz(0.33736526) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0309567) q[1];
sx q[1];
rz(-2.1345761) q[1];
sx q[1];
rz(-0.33302506) q[1];
rz(-pi) q[2];
rz(0.19076244) q[3];
sx q[3];
rz(-1.5883633) q[3];
sx q[3];
rz(2.3004408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11898018) q[2];
sx q[2];
rz(-0.70465124) q[2];
sx q[2];
rz(-1.653778) q[2];
rz(-0.37114272) q[3];
sx q[3];
rz(-1.1266339) q[3];
sx q[3];
rz(-0.77905542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.30598518) q[0];
sx q[0];
rz(-1.7636517) q[0];
sx q[0];
rz(2.4775179) q[0];
rz(2.5750419) q[1];
sx q[1];
rz(-1.605875) q[1];
sx q[1];
rz(-0.18633349) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5390426) q[0];
sx q[0];
rz(-1.7930174) q[0];
sx q[0];
rz(-2.1150388) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8736224) q[2];
sx q[2];
rz(-1.2023965) q[2];
sx q[2];
rz(-0.024193833) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.948958) q[1];
sx q[1];
rz(-1.7342469) q[1];
sx q[1];
rz(0.88731097) q[1];
x q[2];
rz(1.0496907) q[3];
sx q[3];
rz(-1.6288886) q[3];
sx q[3];
rz(-2.5863877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8191007) q[2];
sx q[2];
rz(-1.842061) q[2];
sx q[2];
rz(1.5619649) q[2];
rz(1.9942079) q[3];
sx q[3];
rz(-2.8473144) q[3];
sx q[3];
rz(-1.7779721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.086562432) q[0];
sx q[0];
rz(-1.4921621) q[0];
sx q[0];
rz(-0.30558875) q[0];
rz(-1.2809523) q[1];
sx q[1];
rz(-2.8260904) q[1];
sx q[1];
rz(1.618128) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20210914) q[0];
sx q[0];
rz(-1.3393219) q[0];
sx q[0];
rz(-0.37136308) q[0];
rz(-pi) q[1];
rz(0.32714897) q[2];
sx q[2];
rz(-3.0939348) q[2];
sx q[2];
rz(2.7313155) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.198632) q[1];
sx q[1];
rz(-2.3395798) q[1];
sx q[1];
rz(-0.50870163) q[1];
rz(-0.23766808) q[3];
sx q[3];
rz(-2.1690627) q[3];
sx q[3];
rz(0.054758398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65172255) q[2];
sx q[2];
rz(-1.2469331) q[2];
sx q[2];
rz(2.9215802) q[2];
rz(-0.079719933) q[3];
sx q[3];
rz(-2.1646175) q[3];
sx q[3];
rz(-0.93130934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43059573) q[0];
sx q[0];
rz(-1.3871223) q[0];
sx q[0];
rz(-0.0084477607) q[0];
rz(-1.1620109) q[1];
sx q[1];
rz(-1.9879257) q[1];
sx q[1];
rz(-2.0268424) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8793068) q[0];
sx q[0];
rz(-1.7588076) q[0];
sx q[0];
rz(1.8858869) q[0];
rz(2.5774245) q[2];
sx q[2];
rz(-2.0551066) q[2];
sx q[2];
rz(-2.6251305) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.47059688) q[1];
sx q[1];
rz(-0.45491114) q[1];
sx q[1];
rz(2.2347593) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9825806) q[3];
sx q[3];
rz(-1.1819541) q[3];
sx q[3];
rz(-1.6568615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4343425) q[2];
sx q[2];
rz(-2.12314) q[2];
sx q[2];
rz(-0.84558359) q[2];
rz(-2.8905458) q[3];
sx q[3];
rz(-1.8699402) q[3];
sx q[3];
rz(0.74560753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4276328) q[0];
sx q[0];
rz(-2.7290955) q[0];
sx q[0];
rz(-1.0300256) q[0];
rz(-0.59459844) q[1];
sx q[1];
rz(-1.7183036) q[1];
sx q[1];
rz(1.7880218) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0145688) q[0];
sx q[0];
rz(-1.5480124) q[0];
sx q[0];
rz(1.5915074) q[0];
rz(-2.6313303) q[2];
sx q[2];
rz(-2.7760421) q[2];
sx q[2];
rz(-2.9175959) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8882519) q[1];
sx q[1];
rz(-0.8903114) q[1];
sx q[1];
rz(-2.2254105) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0849801) q[3];
sx q[3];
rz(-1.0407036) q[3];
sx q[3];
rz(0.3584273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2305962) q[2];
sx q[2];
rz(-1.5065008) q[2];
sx q[2];
rz(-0.050749151) q[2];
rz(2.0283608) q[3];
sx q[3];
rz(-1.166393) q[3];
sx q[3];
rz(-0.39065233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0942866) q[0];
sx q[0];
rz(-1.0599437) q[0];
sx q[0];
rz(2.768709) q[0];
rz(1.8376384) q[1];
sx q[1];
rz(-1.8770437) q[1];
sx q[1];
rz(-2.1162927) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7367497) q[0];
sx q[0];
rz(-1.4944634) q[0];
sx q[0];
rz(0.45053225) q[0];
rz(-pi) q[1];
rz(-2.0228902) q[2];
sx q[2];
rz(-2.7774924) q[2];
sx q[2];
rz(-2.5293179) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.655788) q[1];
sx q[1];
rz(-0.88690573) q[1];
sx q[1];
rz(-0.25836103) q[1];
rz(-pi) q[2];
rz(1.9944836) q[3];
sx q[3];
rz(-1.3653339) q[3];
sx q[3];
rz(-2.4567571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9103526) q[2];
sx q[2];
rz(-2.9515036) q[2];
sx q[2];
rz(-1.4300038) q[2];
rz(1.6010239) q[3];
sx q[3];
rz(-0.68683306) q[3];
sx q[3];
rz(-1.471126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98825276) q[0];
sx q[0];
rz(-1.3000458) q[0];
sx q[0];
rz(2.7100995) q[0];
rz(-1.2639812) q[1];
sx q[1];
rz(-1.5411721) q[1];
sx q[1];
rz(-0.65013179) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7509814) q[0];
sx q[0];
rz(-1.6529875) q[0];
sx q[0];
rz(2.85336) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0372978) q[2];
sx q[2];
rz(-0.89020911) q[2];
sx q[2];
rz(-2.7460737) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3172315) q[1];
sx q[1];
rz(-1.6859833) q[1];
sx q[1];
rz(-0.86800602) q[1];
x q[2];
rz(1.5863083) q[3];
sx q[3];
rz(-0.45792327) q[3];
sx q[3];
rz(1.0481038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.55107895) q[2];
sx q[2];
rz(-1.7503909) q[2];
sx q[2];
rz(-3.1371269) q[2];
rz(-0.0095857754) q[3];
sx q[3];
rz(-1.3349345) q[3];
sx q[3];
rz(-0.34346223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9850605) q[0];
sx q[0];
rz(-0.28230202) q[0];
sx q[0];
rz(-2.9811133) q[0];
rz(-1.7566682) q[1];
sx q[1];
rz(-2.3424708) q[1];
sx q[1];
rz(0.30881944) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7910116) q[0];
sx q[0];
rz(-1.6948943) q[0];
sx q[0];
rz(-1.7589633) q[0];
rz(-pi) q[1];
rz(-0.71134348) q[2];
sx q[2];
rz(-1.2738196) q[2];
sx q[2];
rz(-2.1893196) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.05986) q[1];
sx q[1];
rz(-0.55559056) q[1];
sx q[1];
rz(-0.76738552) q[1];
x q[2];
rz(-2.4229523) q[3];
sx q[3];
rz(-2.6860533) q[3];
sx q[3];
rz(-0.18265858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7979692) q[2];
sx q[2];
rz(-2.5324731) q[2];
sx q[2];
rz(1.5416175) q[2];
rz(2.4566417) q[3];
sx q[3];
rz(-1.5545605) q[3];
sx q[3];
rz(-1.6182914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5110382) q[0];
sx q[0];
rz(-1.7031952) q[0];
sx q[0];
rz(-0.55139971) q[0];
rz(-2.664387) q[1];
sx q[1];
rz(-2.2255032) q[1];
sx q[1];
rz(-0.83470693) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3480417) q[0];
sx q[0];
rz(-2.0299596) q[0];
sx q[0];
rz(2.5870917) q[0];
x q[1];
rz(3.1006591) q[2];
sx q[2];
rz(-0.70128317) q[2];
sx q[2];
rz(-2.3778039) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7013423) q[1];
sx q[1];
rz(-1.5363438) q[1];
sx q[1];
rz(-0.95878307) q[1];
rz(2.3605462) q[3];
sx q[3];
rz(-0.78069842) q[3];
sx q[3];
rz(2.4311284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7029552) q[2];
sx q[2];
rz(-2.1795887) q[2];
sx q[2];
rz(0.23923242) q[2];
rz(-0.3244102) q[3];
sx q[3];
rz(-1.7041465) q[3];
sx q[3];
rz(1.602406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5156373) q[0];
sx q[0];
rz(-3.0037168) q[0];
sx q[0];
rz(2.4249518) q[0];
rz(-0.58249885) q[1];
sx q[1];
rz(-1.7889675) q[1];
sx q[1];
rz(-2.8315721) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6577756) q[0];
sx q[0];
rz(-0.80379009) q[0];
sx q[0];
rz(2.0652384) q[0];
rz(-0.76002174) q[2];
sx q[2];
rz(-2.5993002) q[2];
sx q[2];
rz(-3.077728) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0359874) q[1];
sx q[1];
rz(-1.2140555) q[1];
sx q[1];
rz(0.8143266) q[1];
x q[2];
rz(-0.72978348) q[3];
sx q[3];
rz(-0.68396362) q[3];
sx q[3];
rz(-1.7548949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.50905716) q[2];
sx q[2];
rz(-0.57764235) q[2];
sx q[2];
rz(1.7217815) q[2];
rz(-1.696473) q[3];
sx q[3];
rz(-0.71075478) q[3];
sx q[3];
rz(-0.76326171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9025018) q[0];
sx q[0];
rz(-0.9466753) q[0];
sx q[0];
rz(-1.7539903) q[0];
rz(1.4801964) q[1];
sx q[1];
rz(-1.4085242) q[1];
sx q[1];
rz(-2.9396802) q[1];
rz(-0.26726503) q[2];
sx q[2];
rz(-1.4980734) q[2];
sx q[2];
rz(0.43792133) q[2];
rz(2.673951) q[3];
sx q[3];
rz(-0.94398879) q[3];
sx q[3];
rz(2.332609) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
