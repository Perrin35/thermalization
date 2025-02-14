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
rz(2.1287542) q[0];
sx q[0];
rz(-1.1829809) q[0];
sx q[0];
rz(0.24721375) q[0];
rz(1.5565058) q[1];
sx q[1];
rz(-1.0143919) q[1];
sx q[1];
rz(-0.95451626) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1314976) q[0];
sx q[0];
rz(-0.37535497) q[0];
sx q[0];
rz(-1.6502871) q[0];
rz(2.688301) q[2];
sx q[2];
rz(-1.6214091) q[2];
sx q[2];
rz(2.5091189) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.87434972) q[1];
sx q[1];
rz(-1.8863719) q[1];
sx q[1];
rz(-2.3467031) q[1];
x q[2];
rz(0.77311109) q[3];
sx q[3];
rz(-1.793981) q[3];
sx q[3];
rz(3.0497568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9556094) q[2];
sx q[2];
rz(-2.535203) q[2];
sx q[2];
rz(2.8601904) q[2];
rz(1.9143117) q[3];
sx q[3];
rz(-1.3325007) q[3];
sx q[3];
rz(1.9277771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94673741) q[0];
sx q[0];
rz(-0.40732107) q[0];
sx q[0];
rz(2.5545004) q[0];
rz(-2.6214444) q[1];
sx q[1];
rz(-1.9303493) q[1];
sx q[1];
rz(-2.928226) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2706695) q[0];
sx q[0];
rz(-1.6324703) q[0];
sx q[0];
rz(-1.7212409) q[0];
rz(-pi) q[1];
rz(-0.93307067) q[2];
sx q[2];
rz(-2.3847849) q[2];
sx q[2];
rz(0.30100664) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0679464) q[1];
sx q[1];
rz(-1.6151307) q[1];
sx q[1];
rz(-1.4774051) q[1];
rz(2.0265686) q[3];
sx q[3];
rz(-2.5682862) q[3];
sx q[3];
rz(2.440883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1063891) q[2];
sx q[2];
rz(-2.8530402) q[2];
sx q[2];
rz(0.2717379) q[2];
rz(0.64618293) q[3];
sx q[3];
rz(-1.8281507) q[3];
sx q[3];
rz(1.9338231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3216517) q[0];
sx q[0];
rz(-0.8701179) q[0];
sx q[0];
rz(-2.8735549) q[0];
rz(-0.23621121) q[1];
sx q[1];
rz(-1.7832489) q[1];
sx q[1];
rz(1.7265629) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9504323) q[0];
sx q[0];
rz(-1.6163905) q[0];
sx q[0];
rz(3.0395503) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7790952) q[2];
sx q[2];
rz(-1.7300743) q[2];
sx q[2];
rz(0.25937072) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.336795) q[1];
sx q[1];
rz(-2.6141254) q[1];
sx q[1];
rz(2.5673366) q[1];
x q[2];
rz(0.3424267) q[3];
sx q[3];
rz(-2.4391616) q[3];
sx q[3];
rz(0.46063761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1872824) q[2];
sx q[2];
rz(-0.3287181) q[2];
sx q[2];
rz(1.3817361) q[2];
rz(2.7749744) q[3];
sx q[3];
rz(-1.4724052) q[3];
sx q[3];
rz(-0.097675145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5496552) q[0];
sx q[0];
rz(-2.9357935) q[0];
sx q[0];
rz(3.1254712) q[0];
rz(-1.4664117) q[1];
sx q[1];
rz(-1.5212395) q[1];
sx q[1];
rz(2.256567) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8649053) q[0];
sx q[0];
rz(-1.7430796) q[0];
sx q[0];
rz(-0.40934632) q[0];
rz(1.9375291) q[2];
sx q[2];
rz(-1.7842494) q[2];
sx q[2];
rz(2.9652924) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72674561) q[1];
sx q[1];
rz(-1.9512842) q[1];
sx q[1];
rz(0.49729113) q[1];
rz(0.4039558) q[3];
sx q[3];
rz(-1.0896297) q[3];
sx q[3];
rz(-1.4766525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5625988) q[2];
sx q[2];
rz(-2.6851974) q[2];
sx q[2];
rz(-1.558051) q[2];
rz(-1.3715749) q[3];
sx q[3];
rz(-1.2746425) q[3];
sx q[3];
rz(2.9461327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9850995) q[0];
sx q[0];
rz(-1.766196) q[0];
sx q[0];
rz(-0.14955713) q[0];
rz(-2.0984446) q[1];
sx q[1];
rz(-2.1688192) q[1];
sx q[1];
rz(-2.3057888) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5468109) q[0];
sx q[0];
rz(-1.8203837) q[0];
sx q[0];
rz(-1.9052192) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1321394) q[2];
sx q[2];
rz(-0.88047709) q[2];
sx q[2];
rz(-3.0005665) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0440825) q[1];
sx q[1];
rz(-2.9335576) q[1];
sx q[1];
rz(-1.5670304) q[1];
x q[2];
rz(-0.18330611) q[3];
sx q[3];
rz(-1.3245031) q[3];
sx q[3];
rz(0.016427996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0189556) q[2];
sx q[2];
rz(-1.476172) q[2];
sx q[2];
rz(2.7254851) q[2];
rz(-0.082402669) q[3];
sx q[3];
rz(-0.96937886) q[3];
sx q[3];
rz(2.940322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06269726) q[0];
sx q[0];
rz(-2.1128928) q[0];
sx q[0];
rz(2.0507574) q[0];
rz(1.1385607) q[1];
sx q[1];
rz(-1.9416315) q[1];
sx q[1];
rz(-0.60637766) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17601062) q[0];
sx q[0];
rz(-2.2041498) q[0];
sx q[0];
rz(1.8853582) q[0];
rz(2.3396083) q[2];
sx q[2];
rz(-1.238566) q[2];
sx q[2];
rz(1.3114245) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6003638) q[1];
sx q[1];
rz(-2.0906587) q[1];
sx q[1];
rz(-0.71345774) q[1];
x q[2];
rz(1.6929469) q[3];
sx q[3];
rz(-2.2854317) q[3];
sx q[3];
rz(2.2981616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.31327569) q[2];
sx q[2];
rz(-0.54855359) q[2];
sx q[2];
rz(-3.0307148) q[2];
rz(-1.5634792) q[3];
sx q[3];
rz(-2.9302247) q[3];
sx q[3];
rz(-1.3255239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(2.9368847) q[0];
sx q[0];
rz(-0.81975833) q[0];
sx q[0];
rz(-2.4844266) q[0];
rz(1.3232629) q[1];
sx q[1];
rz(-0.75138775) q[1];
sx q[1];
rz(1.125186) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96699636) q[0];
sx q[0];
rz(-1.5856992) q[0];
sx q[0];
rz(3.1257939) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29203307) q[2];
sx q[2];
rz(-0.32675535) q[2];
sx q[2];
rz(0.20694831) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4954056) q[1];
sx q[1];
rz(-2.0584724) q[1];
sx q[1];
rz(1.133092) q[1];
rz(-pi) q[2];
rz(-1.2631038) q[3];
sx q[3];
rz(-2.2224226) q[3];
sx q[3];
rz(3.0857112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7323759) q[2];
sx q[2];
rz(-1.9595307) q[2];
sx q[2];
rz(-2.1803975) q[2];
rz(-3.1327278) q[3];
sx q[3];
rz(-2.2564087) q[3];
sx q[3];
rz(-1.2113021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(0.081414118) q[0];
sx q[0];
rz(-3.0280805) q[0];
sx q[0];
rz(0.42763448) q[0];
rz(1.112452) q[1];
sx q[1];
rz(-1.3149657) q[1];
sx q[1];
rz(2.1449259) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0286197) q[0];
sx q[0];
rz(-2.1899208) q[0];
sx q[0];
rz(2.1834247) q[0];
x q[1];
rz(1.8420503) q[2];
sx q[2];
rz(-0.82640663) q[2];
sx q[2];
rz(0.24559337) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5298047) q[1];
sx q[1];
rz(-2.5080097) q[1];
sx q[1];
rz(-2.0287013) q[1];
rz(-2.5708267) q[3];
sx q[3];
rz(-1.3955355) q[3];
sx q[3];
rz(-1.4890763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4619649) q[2];
sx q[2];
rz(-0.5373911) q[2];
sx q[2];
rz(2.4930387) q[2];
rz(-0.31351659) q[3];
sx q[3];
rz(-2.2532012) q[3];
sx q[3];
rz(-0.052481767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.736883) q[0];
sx q[0];
rz(-0.22060224) q[0];
sx q[0];
rz(2.1746461) q[0];
rz(2.1838358) q[1];
sx q[1];
rz(-1.0916595) q[1];
sx q[1];
rz(1.0831833) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0592902) q[0];
sx q[0];
rz(-1.3246312) q[0];
sx q[0];
rz(2.7013426) q[0];
rz(-pi) q[1];
rz(-0.73019256) q[2];
sx q[2];
rz(-1.0456523) q[2];
sx q[2];
rz(2.5333135) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7249917) q[1];
sx q[1];
rz(-2.3146475) q[1];
sx q[1];
rz(1.0241072) q[1];
rz(-0.47934909) q[3];
sx q[3];
rz(-1.5707301) q[3];
sx q[3];
rz(2.655055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12072418) q[2];
sx q[2];
rz(-2.3126297) q[2];
sx q[2];
rz(-1.6797569) q[2];
rz(0.36137897) q[3];
sx q[3];
rz(-1.8166108) q[3];
sx q[3];
rz(2.1593275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7892889) q[0];
sx q[0];
rz(-0.73065773) q[0];
sx q[0];
rz(-2.5590382) q[0];
rz(1.7763058) q[1];
sx q[1];
rz(-2.1047695) q[1];
sx q[1];
rz(-2.9127311) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6270646) q[0];
sx q[0];
rz(-2.702003) q[0];
sx q[0];
rz(2.702781) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9624356) q[2];
sx q[2];
rz(-1.8300042) q[2];
sx q[2];
rz(-1.1863866) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0693384) q[1];
sx q[1];
rz(-1.4120585) q[1];
sx q[1];
rz(-3.0760514) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7515342) q[3];
sx q[3];
rz(-1.164468) q[3];
sx q[3];
rz(-1.4618946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2710496) q[2];
sx q[2];
rz(-0.55176631) q[2];
sx q[2];
rz(-1.40847) q[2];
rz(0.72077858) q[3];
sx q[3];
rz(-1.9920789) q[3];
sx q[3];
rz(-1.6225947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085717399) q[0];
sx q[0];
rz(-1.7039104) q[0];
sx q[0];
rz(-1.3088551) q[0];
rz(-1.5695288) q[1];
sx q[1];
rz(-0.61620284) q[1];
sx q[1];
rz(-2.9175704) q[1];
rz(-2.8926579) q[2];
sx q[2];
rz(-2.1215083) q[2];
sx q[2];
rz(-0.10100867) q[2];
rz(-0.15275501) q[3];
sx q[3];
rz(-1.7293617) q[3];
sx q[3];
rz(0.21241906) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
