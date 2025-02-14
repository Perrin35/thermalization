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
rz(-1.5850868) q[1];
sx q[1];
rz(-2.1272008) q[1];
sx q[1];
rz(-2.1870764) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2169089) q[0];
sx q[0];
rz(-1.9449073) q[0];
sx q[0];
rz(-3.1103136) q[0];
rz(1.6270832) q[2];
sx q[2];
rz(-1.1181284) q[2];
sx q[2];
rz(-0.91368352) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.87434972) q[1];
sx q[1];
rz(-1.8863719) q[1];
sx q[1];
rz(-2.3467031) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3684816) q[3];
sx q[3];
rz(-1.793981) q[3];
sx q[3];
rz(3.0497568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9556094) q[2];
sx q[2];
rz(-2.535203) q[2];
sx q[2];
rz(0.28140226) q[2];
rz(-1.227281) q[3];
sx q[3];
rz(-1.3325007) q[3];
sx q[3];
rz(-1.2138155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94673741) q[0];
sx q[0];
rz(-0.40732107) q[0];
sx q[0];
rz(2.5545004) q[0];
rz(-0.52014822) q[1];
sx q[1];
rz(-1.9303493) q[1];
sx q[1];
rz(-0.21336666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086047202) q[0];
sx q[0];
rz(-2.9790857) q[0];
sx q[0];
rz(-1.1799728) q[0];
rz(-pi) q[1];
rz(0.92170282) q[2];
sx q[2];
rz(-1.9919167) q[2];
sx q[2];
rz(-2.365961) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6345916) q[1];
sx q[1];
rz(-1.4774972) q[1];
sx q[1];
rz(-3.0970645) q[1];
rz(-1.0454093) q[3];
sx q[3];
rz(-1.3297218) q[3];
sx q[3];
rz(1.8808533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1063891) q[2];
sx q[2];
rz(-0.28855244) q[2];
sx q[2];
rz(2.8698548) q[2];
rz(2.4954097) q[3];
sx q[3];
rz(-1.313442) q[3];
sx q[3];
rz(-1.2077695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81994098) q[0];
sx q[0];
rz(-2.2714748) q[0];
sx q[0];
rz(-0.2680378) q[0];
rz(-0.23621121) q[1];
sx q[1];
rz(-1.3583438) q[1];
sx q[1];
rz(1.4150298) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9504323) q[0];
sx q[0];
rz(-1.5252021) q[0];
sx q[0];
rz(-3.0395503) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4006538) q[2];
sx q[2];
rz(-1.2130951) q[2];
sx q[2];
rz(-1.3715054) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9795003) q[1];
sx q[1];
rz(-2.0071173) q[1];
sx q[1];
rz(-1.8772582) q[1];
rz(-pi) q[2];
rz(-0.67309849) q[3];
sx q[3];
rz(-1.3521225) q[3];
sx q[3];
rz(-0.84450561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.95431027) q[2];
sx q[2];
rz(-0.3287181) q[2];
sx q[2];
rz(1.7598565) q[2];
rz(-0.36661822) q[3];
sx q[3];
rz(-1.6691875) q[3];
sx q[3];
rz(0.097675145) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5919375) q[0];
sx q[0];
rz(-0.20579919) q[0];
sx q[0];
rz(-3.1254712) q[0];
rz(1.6751809) q[1];
sx q[1];
rz(-1.5212395) q[1];
sx q[1];
rz(-0.88502562) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4711759) q[0];
sx q[0];
rz(-0.44222853) q[0];
sx q[0];
rz(2.7294374) q[0];
rz(-1.0270723) q[2];
sx q[2];
rz(-2.7197065) q[2];
sx q[2];
rz(0.89050259) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4440892) q[1];
sx q[1];
rz(-0.61629811) q[1];
sx q[1];
rz(2.4438436) q[1];
rz(-pi) q[2];
rz(-0.92547379) q[3];
sx q[3];
rz(-2.5237115) q[3];
sx q[3];
rz(-2.2223652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9850995) q[0];
sx q[0];
rz(-1.766196) q[0];
sx q[0];
rz(0.14955713) q[0];
rz(-2.0984446) q[1];
sx q[1];
rz(-0.97277343) q[1];
sx q[1];
rz(-0.8358039) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5468109) q[0];
sx q[0];
rz(-1.321209) q[0];
sx q[0];
rz(-1.2363734) q[0];
x q[1];
rz(3.1321394) q[2];
sx q[2];
rz(-2.2611156) q[2];
sx q[2];
rz(3.0005665) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0936613) q[1];
sx q[1];
rz(-1.3627628) q[1];
sx q[1];
rz(-3.1407977) q[1];
rz(-pi) q[2];
rz(1.8211143) q[3];
sx q[3];
rz(-1.7485108) q[3];
sx q[3];
rz(-1.632393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.122637) q[2];
sx q[2];
rz(-1.476172) q[2];
sx q[2];
rz(-2.7254851) q[2];
rz(-3.05919) q[3];
sx q[3];
rz(-0.96937886) q[3];
sx q[3];
rz(0.20127067) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06269726) q[0];
sx q[0];
rz(-2.1128928) q[0];
sx q[0];
rz(-2.0507574) q[0];
rz(-2.003032) q[1];
sx q[1];
rz(-1.1999612) q[1];
sx q[1];
rz(0.60637766) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17601062) q[0];
sx q[0];
rz(-0.93744288) q[0];
sx q[0];
rz(-1.8853582) q[0];
rz(1.1101704) q[2];
sx q[2];
rz(-2.3177882) q[2];
sx q[2];
rz(0.065814171) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55058897) q[1];
sx q[1];
rz(-0.85500756) q[1];
sx q[1];
rz(2.4229933) q[1];
x q[2];
rz(-0.71833937) q[3];
sx q[3];
rz(-1.6629617) q[3];
sx q[3];
rz(0.80764333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.828317) q[2];
sx q[2];
rz(-0.54855359) q[2];
sx q[2];
rz(-3.0307148) q[2];
rz(1.5781135) q[3];
sx q[3];
rz(-0.21136798) q[3];
sx q[3];
rz(-1.8160688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9368847) q[0];
sx q[0];
rz(-0.81975833) q[0];
sx q[0];
rz(-0.65716609) q[0];
rz(1.3232629) q[1];
sx q[1];
rz(-2.3902049) q[1];
sx q[1];
rz(-1.125186) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.781625) q[0];
sx q[0];
rz(-0.021718135) q[0];
sx q[0];
rz(-2.3853073) q[0];
x q[1];
rz(-1.6680587) q[2];
sx q[2];
rz(-1.8832369) q[2];
sx q[2];
rz(0.51431235) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7106428) q[1];
sx q[1];
rz(-0.64326566) q[1];
sx q[1];
rz(0.67420723) q[1];
rz(-pi) q[2];
rz(-1.2631038) q[3];
sx q[3];
rz(-2.2224226) q[3];
sx q[3];
rz(3.0857112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.40921673) q[2];
sx q[2];
rz(-1.9595307) q[2];
sx q[2];
rz(0.96119514) q[2];
rz(0.0088648908) q[3];
sx q[3];
rz(-2.2564087) q[3];
sx q[3];
rz(-1.2113021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.081414118) q[0];
sx q[0];
rz(-3.0280805) q[0];
sx q[0];
rz(2.7139582) q[0];
rz(-2.0291406) q[1];
sx q[1];
rz(-1.8266269) q[1];
sx q[1];
rz(-2.1449259) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11297299) q[0];
sx q[0];
rz(-0.95167187) q[0];
sx q[0];
rz(-2.1834247) q[0];
rz(2.8585343) q[2];
sx q[2];
rz(-2.3583204) q[2];
sx q[2];
rz(0.1439134) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.06293077) q[1];
sx q[1];
rz(-1.0109631) q[1];
sx q[1];
rz(0.31400756) q[1];
x q[2];
rz(0.31672041) q[3];
sx q[3];
rz(-2.5473928) q[3];
sx q[3];
rz(2.7946928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6796278) q[2];
sx q[2];
rz(-2.6042016) q[2];
sx q[2];
rz(-0.648554) q[2];
rz(-2.8280761) q[3];
sx q[3];
rz(-2.2532012) q[3];
sx q[3];
rz(-3.0891109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.736883) q[0];
sx q[0];
rz(-0.22060224) q[0];
sx q[0];
rz(-0.96694651) q[0];
rz(0.95775682) q[1];
sx q[1];
rz(-1.0916595) q[1];
sx q[1];
rz(-1.0831833) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0111128) q[0];
sx q[0];
rz(-0.50043538) q[0];
sx q[0];
rz(-2.6088663) q[0];
rz(-0.90980977) q[2];
sx q[2];
rz(-2.1860115) q[2];
sx q[2];
rz(1.7571952) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5965278) q[1];
sx q[1];
rz(-1.178243) q[1];
sx q[1];
rz(2.3190581) q[1];
rz(3.1414491) q[3];
sx q[3];
rz(-0.4793491) q[3];
sx q[3];
rz(1.0841313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.12072418) q[2];
sx q[2];
rz(-2.3126297) q[2];
sx q[2];
rz(1.6797569) q[2];
rz(0.36137897) q[3];
sx q[3];
rz(-1.3249818) q[3];
sx q[3];
rz(0.98226515) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35230377) q[0];
sx q[0];
rz(-0.73065773) q[0];
sx q[0];
rz(2.5590382) q[0];
rz(-1.3652868) q[1];
sx q[1];
rz(-2.1047695) q[1];
sx q[1];
rz(-2.9127311) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34535566) q[0];
sx q[0];
rz(-1.7526049) q[0];
sx q[0];
rz(0.40248351) q[0];
rz(-pi) q[1];
rz(-0.17915704) q[2];
sx q[2];
rz(-1.3115885) q[2];
sx q[2];
rz(-1.1863866) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6534253) q[1];
sx q[1];
rz(-1.6355124) q[1];
sx q[1];
rz(-1.72987) q[1];
x q[2];
rz(1.1353605) q[3];
sx q[3];
rz(-1.2139911) q[3];
sx q[3];
rz(-0.26998587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.87054306) q[2];
sx q[2];
rz(-2.5898263) q[2];
sx q[2];
rz(1.7331227) q[2];
rz(2.4208141) q[3];
sx q[3];
rz(-1.9920789) q[3];
sx q[3];
rz(-1.518998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0558753) q[0];
sx q[0];
rz(-1.7039104) q[0];
sx q[0];
rz(-1.3088551) q[0];
rz(-1.5695288) q[1];
sx q[1];
rz(-0.61620284) q[1];
sx q[1];
rz(-2.9175704) q[1];
rz(-1.0060251) q[2];
sx q[2];
rz(-1.3592764) q[2];
sx q[2];
rz(-1.5395561) q[2];
rz(1.7311981) q[3];
sx q[3];
rz(-1.7216202) q[3];
sx q[3];
rz(-1.3826821) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
