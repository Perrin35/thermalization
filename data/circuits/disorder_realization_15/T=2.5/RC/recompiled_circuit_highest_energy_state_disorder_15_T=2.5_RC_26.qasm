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
rz(0.61007845) q[0];
sx q[0];
rz(3.371513) q[0];
sx q[0];
rz(11.846677) q[0];
rz(2.5201058) q[1];
sx q[1];
rz(-1.7401594) q[1];
sx q[1];
rz(-3.016234) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7490112) q[0];
sx q[0];
rz(-1.1454937) q[0];
sx q[0];
rz(2.0122737) q[0];
rz(-1.6654451) q[2];
sx q[2];
rz(-1.2834335) q[2];
sx q[2];
rz(-0.3269302) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.6290063) q[1];
sx q[1];
rz(-0.0014481469) q[1];
sx q[1];
rz(-1.8667029) q[1];
x q[2];
rz(-1.9699279) q[3];
sx q[3];
rz(-2.3775775) q[3];
sx q[3];
rz(-1.7872052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7980935) q[2];
sx q[2];
rz(-2.7332879) q[2];
sx q[2];
rz(-2.2957392) q[2];
rz(-0.79126233) q[3];
sx q[3];
rz(-3.1281804) q[3];
sx q[3];
rz(-3.0748034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4132408) q[0];
sx q[0];
rz(-0.46717307) q[0];
sx q[0];
rz(-3.0979284) q[0];
rz(-1.5664258) q[1];
sx q[1];
rz(-1.7685726) q[1];
sx q[1];
rz(-1.498819) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9696143) q[0];
sx q[0];
rz(-0.96706264) q[0];
sx q[0];
rz(-0.0049334299) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0129406) q[2];
sx q[2];
rz(-2.1421551) q[2];
sx q[2];
rz(-0.017627942) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1825492) q[1];
sx q[1];
rz(-1.5397397) q[1];
sx q[1];
rz(-3.0654383) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2263857) q[3];
sx q[3];
rz(-1.5395172) q[3];
sx q[3];
rz(2.7117604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0955536) q[2];
sx q[2];
rz(-2.9914896) q[2];
sx q[2];
rz(-2.6138439) q[2];
rz(2.2857417) q[3];
sx q[3];
rz(-0.0015043613) q[3];
sx q[3];
rz(1.2214448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7432778) q[0];
sx q[0];
rz(-2.1687431) q[0];
sx q[0];
rz(1.995218) q[0];
rz(-1.758681) q[1];
sx q[1];
rz(-2.8490366) q[1];
sx q[1];
rz(0.10398277) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2706174) q[0];
sx q[0];
rz(-1.7954602) q[0];
sx q[0];
rz(-1.490834) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5095897) q[2];
sx q[2];
rz(-2.8687697) q[2];
sx q[2];
rz(-0.83977985) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.79499352) q[1];
sx q[1];
rz(-1.3665587) q[1];
sx q[1];
rz(-0.7612106) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29071112) q[3];
sx q[3];
rz(-1.4658648) q[3];
sx q[3];
rz(-0.19874979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9581703) q[2];
sx q[2];
rz(-3.1347745) q[2];
sx q[2];
rz(0.56162322) q[2];
rz(-3.0567452) q[3];
sx q[3];
rz(-3.1361129) q[3];
sx q[3];
rz(0.00076278846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4168004) q[0];
sx q[0];
rz(-0.264072) q[0];
sx q[0];
rz(-0.52077878) q[0];
rz(0.15781038) q[1];
sx q[1];
rz(-0.66693711) q[1];
sx q[1];
rz(-3.0657943) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4593845) q[0];
sx q[0];
rz(-0.59218107) q[0];
sx q[0];
rz(-0.20488744) q[0];
rz(-pi) q[1];
rz(1.130777) q[2];
sx q[2];
rz(-0.001507757) q[2];
sx q[2];
rz(-1.8720117) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.27386777) q[1];
sx q[1];
rz(-0.40479044) q[1];
sx q[1];
rz(1.2048436) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10785477) q[3];
sx q[3];
rz(-2.0216935) q[3];
sx q[3];
rz(-1.0552561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51125222) q[2];
sx q[2];
rz(-3.1255836) q[2];
sx q[2];
rz(-1.7208257) q[2];
rz(-0.00682791) q[3];
sx q[3];
rz(-0.029191645) q[3];
sx q[3];
rz(1.6543057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1348949) q[0];
sx q[0];
rz(-0.20773523) q[0];
sx q[0];
rz(-0.2625221) q[0];
rz(0.94995704) q[1];
sx q[1];
rz(-0.078153178) q[1];
sx q[1];
rz(2.9238759) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9757468) q[0];
sx q[0];
rz(-1.4474156) q[0];
sx q[0];
rz(1.1749772) q[0];
x q[1];
rz(-0.76518671) q[2];
sx q[2];
rz(-0.12066855) q[2];
sx q[2];
rz(2.3007002) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9432944) q[1];
sx q[1];
rz(-1.5729331) q[1];
sx q[1];
rz(1.398477) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0988555) q[3];
sx q[3];
rz(-1.510757) q[3];
sx q[3];
rz(2.754902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9265284) q[2];
sx q[2];
rz(-0.0096409163) q[2];
sx q[2];
rz(2.8936774) q[2];
rz(1.3748112) q[3];
sx q[3];
rz(-3.091076) q[3];
sx q[3];
rz(-1.7063399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0404102) q[0];
sx q[0];
rz(-1.8721767) q[0];
sx q[0];
rz(2.5293479) q[0];
rz(2.9712037) q[1];
sx q[1];
rz(-0.082516106) q[1];
sx q[1];
rz(1.6498529) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1338324) q[0];
sx q[0];
rz(-2.2585224) q[0];
sx q[0];
rz(2.4439815) q[0];
rz(-pi) q[1];
rz(-3.1270624) q[2];
sx q[2];
rz(-1.5749914) q[2];
sx q[2];
rz(2.1452877) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.41497624) q[1];
sx q[1];
rz(-1.5041989) q[1];
sx q[1];
rz(-2.8482751) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.065656) q[3];
sx q[3];
rz(-2.7759984) q[3];
sx q[3];
rz(2.5557809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6461569) q[2];
sx q[2];
rz(-3.1313681) q[2];
sx q[2];
rz(-2.9025027) q[2];
rz(2.332989) q[3];
sx q[3];
rz(-3.1294398) q[3];
sx q[3];
rz(1.3330207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7072153) q[0];
sx q[0];
rz(-0.093952976) q[0];
sx q[0];
rz(0.77675003) q[0];
rz(3.1306733) q[1];
sx q[1];
rz(-0.24591406) q[1];
sx q[1];
rz(1.4728665) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4368255) q[0];
sx q[0];
rz(-0.87427199) q[0];
sx q[0];
rz(0.27946194) q[0];
rz(-pi) q[1];
rz(1.5644349) q[2];
sx q[2];
rz(-1.5943051) q[2];
sx q[2];
rz(0.17447093) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.90172988) q[1];
sx q[1];
rz(-2.3119676) q[1];
sx q[1];
rz(-1.7077441) q[1];
rz(-pi) q[2];
rz(1.8223354) q[3];
sx q[3];
rz(-2.8467379) q[3];
sx q[3];
rz(0.40094646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9933068) q[2];
sx q[2];
rz(-3.1243656) q[2];
sx q[2];
rz(-0.41603184) q[2];
rz(-2.3038583) q[3];
sx q[3];
rz(-3.0031524) q[3];
sx q[3];
rz(1.4778888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96529043) q[0];
sx q[0];
rz(-3.1017113) q[0];
sx q[0];
rz(-2.1862929) q[0];
rz(1.2166066) q[1];
sx q[1];
rz(-2.80426) q[1];
sx q[1];
rz(0.40562707) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0829725) q[0];
sx q[0];
rz(-2.342413) q[0];
sx q[0];
rz(0.45244502) q[0];
rz(-2.9826016) q[2];
sx q[2];
rz(-1.482725) q[2];
sx q[2];
rz(-0.51364567) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7136399) q[1];
sx q[1];
rz(-1.7312839) q[1];
sx q[1];
rz(-1.6482216) q[1];
x q[2];
rz(-0.21307034) q[3];
sx q[3];
rz(-1.4858507) q[3];
sx q[3];
rz(1.5379048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2847298) q[2];
sx q[2];
rz(-3.1046107) q[2];
sx q[2];
rz(0.82376087) q[2];
rz(-3.01037) q[3];
sx q[3];
rz(-2.8516912) q[3];
sx q[3];
rz(-2.8141865) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71309483) q[0];
sx q[0];
rz(-2.9976124) q[0];
sx q[0];
rz(-0.71260989) q[0];
rz(-0.54829848) q[1];
sx q[1];
rz(-0.24990853) q[1];
sx q[1];
rz(0.20864329) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31748235) q[0];
sx q[0];
rz(-1.1718996) q[0];
sx q[0];
rz(2.0092416) q[0];
rz(-1.6054262) q[2];
sx q[2];
rz(-1.5896665) q[2];
sx q[2];
rz(0.46771995) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.123173) q[1];
sx q[1];
rz(-0.55572617) q[1];
sx q[1];
rz(-1.5617784) q[1];
rz(-pi) q[2];
rz(-2.3712158) q[3];
sx q[3];
rz(-2.340014) q[3];
sx q[3];
rz(-1.0018536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.448552) q[2];
sx q[2];
rz(-3.0604477) q[2];
sx q[2];
rz(-2.9244259) q[2];
rz(0.36755696) q[3];
sx q[3];
rz(-3.1072072) q[3];
sx q[3];
rz(2.1448081) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17598584) q[0];
sx q[0];
rz(-0.088986926) q[0];
sx q[0];
rz(-2.9678645) q[0];
rz(-1.4088176) q[1];
sx q[1];
rz(-1.4585739) q[1];
sx q[1];
rz(1.6857612) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30688863) q[0];
sx q[0];
rz(-1.406989) q[0];
sx q[0];
rz(2.6192946) q[0];
rz(-pi) q[1];
rz(-2.1591804) q[2];
sx q[2];
rz(-0.2478226) q[2];
sx q[2];
rz(0.10977015) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5816702) q[1];
sx q[1];
rz(-1.6877203) q[1];
sx q[1];
rz(-0.70150866) q[1];
rz(1.9718658) q[3];
sx q[3];
rz(-2.8725343) q[3];
sx q[3];
rz(2.6869686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.17560656) q[2];
sx q[2];
rz(-2.6484428) q[2];
sx q[2];
rz(-1.7612339) q[2];
rz(0.68584758) q[3];
sx q[3];
rz(-3.139747) q[3];
sx q[3];
rz(2.4628911) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3995517) q[0];
sx q[0];
rz(-2.4492332) q[0];
sx q[0];
rz(-1.5204182) q[0];
rz(-1.5834658) q[1];
sx q[1];
rz(-1.5065267) q[1];
sx q[1];
rz(-2.9361257) q[1];
rz(1.4420527) q[2];
sx q[2];
rz(-3.0151571) q[2];
sx q[2];
rz(-2.7966316) q[2];
rz(-1.2711502) q[3];
sx q[3];
rz(-1.4574403) q[3];
sx q[3];
rz(1.311655) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
