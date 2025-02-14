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
rz(-2.8943789) q[0];
rz(1.5565058) q[1];
sx q[1];
rz(-1.0143919) q[1];
sx q[1];
rz(-0.95451626) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2169089) q[0];
sx q[0];
rz(-1.1966853) q[0];
sx q[0];
rz(-3.1103136) q[0];
rz(-1.5145095) q[2];
sx q[2];
rz(-2.0234642) q[2];
sx q[2];
rz(0.91368352) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2672429) q[1];
sx q[1];
rz(-1.8863719) q[1];
sx q[1];
rz(0.79488956) q[1];
rz(-1.263721) q[3];
sx q[3];
rz(-0.82160866) q[3];
sx q[3];
rz(1.6916569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9556094) q[2];
sx q[2];
rz(-2.535203) q[2];
sx q[2];
rz(-2.8601904) q[2];
rz(-1.227281) q[3];
sx q[3];
rz(-1.3325007) q[3];
sx q[3];
rz(1.9277771) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1948552) q[0];
sx q[0];
rz(-0.40732107) q[0];
sx q[0];
rz(2.5545004) q[0];
rz(-2.6214444) q[1];
sx q[1];
rz(-1.9303493) q[1];
sx q[1];
rz(0.21336666) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8321229) q[0];
sx q[0];
rz(-1.42064) q[0];
sx q[0];
rz(3.0792159) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92170282) q[2];
sx q[2];
rz(-1.1496759) q[2];
sx q[2];
rz(-0.77563167) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.060886968) q[1];
sx q[1];
rz(-0.10335246) q[1];
sx q[1];
rz(-1.1267613) q[1];
rz(-pi) q[2];
rz(-2.8647086) q[3];
sx q[3];
rz(-2.0794984) q[3];
sx q[3];
rz(0.17252082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0352036) q[2];
sx q[2];
rz(-2.8530402) q[2];
sx q[2];
rz(2.8698548) q[2];
rz(-0.64618293) q[3];
sx q[3];
rz(-1.8281507) q[3];
sx q[3];
rz(-1.9338231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3216517) q[0];
sx q[0];
rz(-0.8701179) q[0];
sx q[0];
rz(-0.2680378) q[0];
rz(-2.9053814) q[1];
sx q[1];
rz(-1.3583438) q[1];
sx q[1];
rz(-1.4150298) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9504323) q[0];
sx q[0];
rz(-1.5252021) q[0];
sx q[0];
rz(3.0395503) q[0];
rz(0.36249749) q[2];
sx q[2];
rz(-1.4115184) q[2];
sx q[2];
rz(0.25937072) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9795003) q[1];
sx q[1];
rz(-1.1344754) q[1];
sx q[1];
rz(-1.2643344) q[1];
x q[2];
rz(-0.67309849) q[3];
sx q[3];
rz(-1.7894701) q[3];
sx q[3];
rz(0.84450561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.95431027) q[2];
sx q[2];
rz(-0.3287181) q[2];
sx q[2];
rz(-1.3817361) q[2];
rz(0.36661822) q[3];
sx q[3];
rz(-1.4724052) q[3];
sx q[3];
rz(0.097675145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.5919375) q[0];
sx q[0];
rz(-2.9357935) q[0];
sx q[0];
rz(0.016121443) q[0];
rz(-1.4664117) q[1];
sx q[1];
rz(-1.5212395) q[1];
sx q[1];
rz(-0.88502562) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21986905) q[0];
sx q[0];
rz(-1.1678639) q[0];
sx q[0];
rz(1.7582488) q[0];
rz(-2.1145203) q[2];
sx q[2];
rz(-2.7197065) q[2];
sx q[2];
rz(2.2510901) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.72674561) q[1];
sx q[1];
rz(-1.9512842) q[1];
sx q[1];
rz(2.6443015) q[1];
rz(-pi) q[2];
rz(-2.0871986) q[3];
sx q[3];
rz(-1.9266911) q[3];
sx q[3];
rz(-3.0404224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.57899388) q[2];
sx q[2];
rz(-0.4563953) q[2];
sx q[2];
rz(1.558051) q[2];
rz(1.3715749) q[3];
sx q[3];
rz(-1.8669502) q[3];
sx q[3];
rz(-0.19545999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9850995) q[0];
sx q[0];
rz(-1.3753966) q[0];
sx q[0];
rz(-0.14955713) q[0];
rz(2.0984446) q[1];
sx q[1];
rz(-2.1688192) q[1];
sx q[1];
rz(2.3057888) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1096032) q[0];
sx q[0];
rz(-1.8944724) q[0];
sx q[0];
rz(2.8780185) q[0];
rz(-pi) q[1];
rz(0.0094532254) q[2];
sx q[2];
rz(-0.88047709) q[2];
sx q[2];
rz(-0.1410262) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52302915) q[1];
sx q[1];
rz(-1.5715741) q[1];
sx q[1];
rz(1.3627627) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1981524) q[3];
sx q[3];
rz(-0.30590484) q[3];
sx q[3];
rz(2.4750575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.122637) q[2];
sx q[2];
rz(-1.476172) q[2];
sx q[2];
rz(2.7254851) q[2];
rz(0.082402669) q[3];
sx q[3];
rz(-2.1722138) q[3];
sx q[3];
rz(-0.20127067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0788954) q[0];
sx q[0];
rz(-1.0286999) q[0];
sx q[0];
rz(2.0507574) q[0];
rz(-2.003032) q[1];
sx q[1];
rz(-1.9416315) q[1];
sx q[1];
rz(-0.60637766) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9370458) q[0];
sx q[0];
rz(-1.318745) q[0];
sx q[0];
rz(-2.4840647) q[0];
rz(-pi) q[1];
rz(-1.1101704) q[2];
sx q[2];
rz(-2.3177882) q[2];
sx q[2];
rz(-0.065814171) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.37651) q[1];
sx q[1];
rz(-0.9667337) q[1];
sx q[1];
rz(2.2187697) q[1];
rz(-3.0020671) q[3];
sx q[3];
rz(-2.4184113) q[3];
sx q[3];
rz(-2.4833402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.828317) q[2];
sx q[2];
rz(-2.5930391) q[2];
sx q[2];
rz(0.11087785) q[2];
rz(-1.5781135) q[3];
sx q[3];
rz(-0.21136798) q[3];
sx q[3];
rz(1.8160688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.204708) q[0];
sx q[0];
rz(-0.81975833) q[0];
sx q[0];
rz(2.4844266) q[0];
rz(1.3232629) q[1];
sx q[1];
rz(-2.3902049) q[1];
sx q[1];
rz(2.0164067) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1745963) q[0];
sx q[0];
rz(-1.5856992) q[0];
sx q[0];
rz(3.1257939) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4735339) q[2];
sx q[2];
rz(-1.8832369) q[2];
sx q[2];
rz(-2.6272803) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.64618707) q[1];
sx q[1];
rz(-2.0584724) q[1];
sx q[1];
rz(-1.133092) q[1];
rz(-pi) q[2];
rz(2.7636307) q[3];
sx q[3];
rz(-0.71092793) q[3];
sx q[3];
rz(2.7148249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7323759) q[2];
sx q[2];
rz(-1.9595307) q[2];
sx q[2];
rz(2.1803975) q[2];
rz(-0.0088648908) q[3];
sx q[3];
rz(-2.2564087) q[3];
sx q[3];
rz(1.2113021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081414118) q[0];
sx q[0];
rz(-3.0280805) q[0];
sx q[0];
rz(2.7139582) q[0];
rz(2.0291406) q[1];
sx q[1];
rz(-1.3149657) q[1];
sx q[1];
rz(0.99666673) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2964945) q[0];
sx q[0];
rz(-2.0581492) q[0];
sx q[0];
rz(2.4250406) q[0];
rz(-pi) q[1];
rz(-2.3786167) q[2];
sx q[2];
rz(-1.7691649) q[2];
sx q[2];
rz(-1.5114443) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0786619) q[1];
sx q[1];
rz(-2.1306296) q[1];
sx q[1];
rz(-2.8275851) q[1];
rz(2.8248722) q[3];
sx q[3];
rz(-2.5473928) q[3];
sx q[3];
rz(-2.7946928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6796278) q[2];
sx q[2];
rz(-2.6042016) q[2];
sx q[2];
rz(0.648554) q[2];
rz(2.8280761) q[3];
sx q[3];
rz(-0.88839141) q[3];
sx q[3];
rz(0.052481767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40470966) q[0];
sx q[0];
rz(-0.22060224) q[0];
sx q[0];
rz(0.96694651) q[0];
rz(-2.1838358) q[1];
sx q[1];
rz(-2.0499332) q[1];
sx q[1];
rz(1.0831833) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0592902) q[0];
sx q[0];
rz(-1.8169615) q[0];
sx q[0];
rz(2.7013426) q[0];
x q[1];
rz(2.2317829) q[2];
sx q[2];
rz(-2.1860115) q[2];
sx q[2];
rz(1.7571952) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.41660096) q[1];
sx q[1];
rz(-0.82694516) q[1];
sx q[1];
rz(-1.0241072) q[1];
rz(1.5708709) q[3];
sx q[3];
rz(-1.0914472) q[3];
sx q[3];
rz(1.0842931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.12072418) q[2];
sx q[2];
rz(-0.82896295) q[2];
sx q[2];
rz(-1.6797569) q[2];
rz(-2.7802137) q[3];
sx q[3];
rz(-1.8166108) q[3];
sx q[3];
rz(2.1593275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35230377) q[0];
sx q[0];
rz(-2.4109349) q[0];
sx q[0];
rz(-2.5590382) q[0];
rz(1.7763058) q[1];
sx q[1];
rz(-2.1047695) q[1];
sx q[1];
rz(0.22886151) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.992975) q[0];
sx q[0];
rz(-1.1753191) q[0];
sx q[0];
rz(-1.7680042) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8340298) q[2];
sx q[2];
rz(-1.7439067) q[2];
sx q[2];
rz(0.43079475) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0693384) q[1];
sx q[1];
rz(-1.4120585) q[1];
sx q[1];
rz(-3.0760514) q[1];
rz(-pi) q[2];
rz(-2.2945461) q[3];
sx q[3];
rz(-0.55560613) q[3];
sx q[3];
rz(0.65680066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87054306) q[2];
sx q[2];
rz(-2.5898263) q[2];
sx q[2];
rz(1.40847) q[2];
rz(2.4208141) q[3];
sx q[3];
rz(-1.9920789) q[3];
sx q[3];
rz(-1.518998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0558753) q[0];
sx q[0];
rz(-1.7039104) q[0];
sx q[0];
rz(-1.3088551) q[0];
rz(1.5720639) q[1];
sx q[1];
rz(-0.61620284) q[1];
sx q[1];
rz(-2.9175704) q[1];
rz(2.8926579) q[2];
sx q[2];
rz(-1.0200844) q[2];
sx q[2];
rz(3.040584) q[2];
rz(-2.3313777) q[3];
sx q[3];
rz(-0.21972909) q[3];
sx q[3];
rz(-0.56032203) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
