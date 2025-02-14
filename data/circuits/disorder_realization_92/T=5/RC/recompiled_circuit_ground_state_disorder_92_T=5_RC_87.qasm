OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7654704) q[0];
sx q[0];
rz(-2.0630615) q[0];
sx q[0];
rz(-1.7662788) q[0];
rz(-2.7219661) q[1];
sx q[1];
rz(-1.3803177) q[1];
sx q[1];
rz(2.3461411) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4856483) q[0];
sx q[0];
rz(-2.7599381) q[0];
sx q[0];
rz(-0.21393805) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8654557) q[2];
sx q[2];
rz(-0.41414663) q[2];
sx q[2];
rz(-0.93246952) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7401984) q[1];
sx q[1];
rz(-0.90448442) q[1];
sx q[1];
rz(1.3191651) q[1];
rz(-1.8303371) q[3];
sx q[3];
rz(-1.653977) q[3];
sx q[3];
rz(-0.2587425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2835283) q[2];
sx q[2];
rz(-2.8499446) q[2];
sx q[2];
rz(0.38823271) q[2];
rz(-2.9325716) q[3];
sx q[3];
rz(-1.4744604) q[3];
sx q[3];
rz(2.2504263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46880576) q[0];
sx q[0];
rz(-0.46872941) q[0];
sx q[0];
rz(-0.76954532) q[0];
rz(2.1107213) q[1];
sx q[1];
rz(-2.2831235) q[1];
sx q[1];
rz(1.7147725) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010410943) q[0];
sx q[0];
rz(-0.37574925) q[0];
sx q[0];
rz(0.6300169) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67536036) q[2];
sx q[2];
rz(-1.5262233) q[2];
sx q[2];
rz(1.7365484) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.64411341) q[1];
sx q[1];
rz(-2.1783324) q[1];
sx q[1];
rz(-1.2588905) q[1];
rz(2.2852422) q[3];
sx q[3];
rz(-1.6110241) q[3];
sx q[3];
rz(0.86733199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1398805) q[2];
sx q[2];
rz(-2.6931245) q[2];
sx q[2];
rz(-1.9072388) q[2];
rz(-2.6185696) q[3];
sx q[3];
rz(-1.1250027) q[3];
sx q[3];
rz(-2.4856868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6958385) q[0];
sx q[0];
rz(-0.0025175968) q[0];
sx q[0];
rz(-2.5939831) q[0];
rz(-2.2721263) q[1];
sx q[1];
rz(-2.165386) q[1];
sx q[1];
rz(-1.4617823) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6610091) q[0];
sx q[0];
rz(-2.5877773) q[0];
sx q[0];
rz(0.69032918) q[0];
rz(2.3794054) q[2];
sx q[2];
rz(-1.3831639) q[2];
sx q[2];
rz(2.4248276) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.96986497) q[1];
sx q[1];
rz(-1.1369924) q[1];
sx q[1];
rz(2.8126905) q[1];
rz(-pi) q[2];
rz(-0.84949228) q[3];
sx q[3];
rz(-0.85334001) q[3];
sx q[3];
rz(0.8161374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.72472858) q[2];
sx q[2];
rz(-0.29837307) q[2];
sx q[2];
rz(-2.5073063) q[2];
rz(2.0554845) q[3];
sx q[3];
rz(-2.6606798) q[3];
sx q[3];
rz(-2.0079131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5241765) q[0];
sx q[0];
rz(-0.15707459) q[0];
sx q[0];
rz(-1.6198535) q[0];
rz(-0.98328439) q[1];
sx q[1];
rz(-1.1055999) q[1];
sx q[1];
rz(0.081667893) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0591563) q[0];
sx q[0];
rz(-0.56917039) q[0];
sx q[0];
rz(1.3630483) q[0];
rz(-1.5310202) q[2];
sx q[2];
rz(-0.33333862) q[2];
sx q[2];
rz(-3.1266629) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.28934902) q[1];
sx q[1];
rz(-2.056958) q[1];
sx q[1];
rz(0.48883171) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6713167) q[3];
sx q[3];
rz(-2.649622) q[3];
sx q[3];
rz(0.63382404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6433158) q[2];
sx q[2];
rz(-1.4772819) q[2];
sx q[2];
rz(-0.91900438) q[2];
rz(2.8194341) q[3];
sx q[3];
rz(-1.6291658) q[3];
sx q[3];
rz(-0.34644103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0568986) q[0];
sx q[0];
rz(-2.0094805) q[0];
sx q[0];
rz(1.803501) q[0];
rz(-0.40254205) q[1];
sx q[1];
rz(-1.891338) q[1];
sx q[1];
rz(-1.58443) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65619233) q[0];
sx q[0];
rz(-1.3043405) q[0];
sx q[0];
rz(1.043826) q[0];
x q[1];
rz(2.3473559) q[2];
sx q[2];
rz(-1.4068687) q[2];
sx q[2];
rz(1.9203223) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.63478032) q[1];
sx q[1];
rz(-2.1012329) q[1];
sx q[1];
rz(-1.7407559) q[1];
x q[2];
rz(-1.1013056) q[3];
sx q[3];
rz(-1.3302505) q[3];
sx q[3];
rz(2.8194373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0255787) q[2];
sx q[2];
rz(-1.8639001) q[2];
sx q[2];
rz(-2.4063827) q[2];
rz(-1.0544581) q[3];
sx q[3];
rz(-2.1078883) q[3];
sx q[3];
rz(-1.2247156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9758107) q[0];
sx q[0];
rz(-0.42401341) q[0];
sx q[0];
rz(-1.8884416) q[0];
rz(1.5755298) q[1];
sx q[1];
rz(-2.0063446) q[1];
sx q[1];
rz(2.6083686) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3778702) q[0];
sx q[0];
rz(-1.8061751) q[0];
sx q[0];
rz(-0.81416582) q[0];
rz(0.27006939) q[2];
sx q[2];
rz(-0.52723072) q[2];
sx q[2];
rz(-2.466885) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3371083) q[1];
sx q[1];
rz(-1.7586367) q[1];
sx q[1];
rz(-0.2262745) q[1];
rz(0.76119555) q[3];
sx q[3];
rz(-0.71925876) q[3];
sx q[3];
rz(-1.0750036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.48625654) q[2];
sx q[2];
rz(-2.76261) q[2];
sx q[2];
rz(0.8302702) q[2];
rz(1.8298979) q[3];
sx q[3];
rz(-2.0528841) q[3];
sx q[3];
rz(-2.8970498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.526392) q[0];
sx q[0];
rz(-1.4854234) q[0];
sx q[0];
rz(-0.8335337) q[0];
rz(2.1444881) q[1];
sx q[1];
rz(-1.763696) q[1];
sx q[1];
rz(1.5901784) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7187945) q[0];
sx q[0];
rz(-1.4961406) q[0];
sx q[0];
rz(-0.43674998) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0710414) q[2];
sx q[2];
rz(-0.65325551) q[2];
sx q[2];
rz(2.5219301) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7128858) q[1];
sx q[1];
rz(-0.16413153) q[1];
sx q[1];
rz(-0.81961373) q[1];
rz(-pi) q[2];
rz(1.3209913) q[3];
sx q[3];
rz(-2.0916998) q[3];
sx q[3];
rz(-1.5438207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3414063) q[2];
sx q[2];
rz(-2.1087346) q[2];
sx q[2];
rz(1.3979073) q[2];
rz(-2.7684815) q[3];
sx q[3];
rz(-0.92445508) q[3];
sx q[3];
rz(-1.0188518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43847346) q[0];
sx q[0];
rz(-1.6699474) q[0];
sx q[0];
rz(-0.31914172) q[0];
rz(-1.1646264) q[1];
sx q[1];
rz(-1.9705557) q[1];
sx q[1];
rz(3.107531) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8933786) q[0];
sx q[0];
rz(-2.3015018) q[0];
sx q[0];
rz(2.7284751) q[0];
x q[1];
rz(-0.77000846) q[2];
sx q[2];
rz(-1.3427375) q[2];
sx q[2];
rz(2.8213825) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0314591) q[1];
sx q[1];
rz(-1.5601399) q[1];
sx q[1];
rz(-1.5926804) q[1];
x q[2];
rz(-0.56896992) q[3];
sx q[3];
rz(-2.5317857) q[3];
sx q[3];
rz(0.43492521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.065319149) q[2];
sx q[2];
rz(-1.7460145) q[2];
sx q[2];
rz(3.0714387) q[2];
rz(1.8326727) q[3];
sx q[3];
rz(-0.85246837) q[3];
sx q[3];
rz(-3.0064228) q[3];
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
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0604621) q[0];
sx q[0];
rz(-1.9996996) q[0];
sx q[0];
rz(-2.6666226) q[0];
rz(1.9604856) q[1];
sx q[1];
rz(-2.2468552) q[1];
sx q[1];
rz(-2.2244577) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.006893039) q[0];
sx q[0];
rz(-1.0118766) q[0];
sx q[0];
rz(-1.8447645) q[0];
rz(-pi) q[1];
rz(-1.9462637) q[2];
sx q[2];
rz(-2.4175543) q[2];
sx q[2];
rz(2.6065116) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0168382) q[1];
sx q[1];
rz(-1.4887344) q[1];
sx q[1];
rz(1.7793307) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9539509) q[3];
sx q[3];
rz(-2.1360548) q[3];
sx q[3];
rz(1.3711898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8721623) q[2];
sx q[2];
rz(-1.3698801) q[2];
sx q[2];
rz(-1.5577462) q[2];
rz(-2.9861084) q[3];
sx q[3];
rz(-1.0289861) q[3];
sx q[3];
rz(2.1492929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.48142049) q[0];
sx q[0];
rz(-1.6654797) q[0];
sx q[0];
rz(2.7760264) q[0];
rz(1.2008601) q[1];
sx q[1];
rz(-1.5301842) q[1];
sx q[1];
rz(-1.0467122) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7378254) q[0];
sx q[0];
rz(-1.0407454) q[0];
sx q[0];
rz(2.5395495) q[0];
rz(-0.44136845) q[2];
sx q[2];
rz(-1.0689931) q[2];
sx q[2];
rz(-2.8428915) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9871708) q[1];
sx q[1];
rz(-0.94564309) q[1];
sx q[1];
rz(-2.465017) q[1];
rz(0.58308122) q[3];
sx q[3];
rz(-1.1016535) q[3];
sx q[3];
rz(-2.6173979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4552348) q[2];
sx q[2];
rz(-2.0265667) q[2];
sx q[2];
rz(2.4380747) q[2];
rz(-2.5840058) q[3];
sx q[3];
rz(-2.360354) q[3];
sx q[3];
rz(-0.47344661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8483509) q[0];
sx q[0];
rz(-1.0307172) q[0];
sx q[0];
rz(0.95681554) q[0];
rz(2.1059857) q[1];
sx q[1];
rz(-2.5685812) q[1];
sx q[1];
rz(-2.0249637) q[1];
rz(-2.507092) q[2];
sx q[2];
rz(-2.366323) q[2];
sx q[2];
rz(0.97699588) q[2];
rz(-1.5286636) q[3];
sx q[3];
rz(-1.7880472) q[3];
sx q[3];
rz(0.40732297) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
