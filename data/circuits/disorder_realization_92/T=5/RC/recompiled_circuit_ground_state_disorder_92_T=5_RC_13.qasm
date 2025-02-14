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
rz(0.41962656) q[1];
sx q[1];
rz(-1.761275) q[1];
sx q[1];
rz(0.79545155) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6559443) q[0];
sx q[0];
rz(-2.7599381) q[0];
sx q[0];
rz(0.21393805) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8654557) q[2];
sx q[2];
rz(-0.41414663) q[2];
sx q[2];
rz(0.93246952) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7955543) q[1];
sx q[1];
rz(-2.4361894) q[1];
sx q[1];
rz(2.8349174) q[1];
x q[2];
rz(1.8849089) q[3];
sx q[3];
rz(-0.2722578) q[3];
sx q[3];
rz(-1.5262878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2835283) q[2];
sx q[2];
rz(-2.8499446) q[2];
sx q[2];
rz(-2.7533599) q[2];
rz(-2.9325716) q[3];
sx q[3];
rz(-1.4744604) q[3];
sx q[3];
rz(-0.89116636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6727869) q[0];
sx q[0];
rz(-0.46872941) q[0];
sx q[0];
rz(0.76954532) q[0];
rz(1.0308713) q[1];
sx q[1];
rz(-0.85846916) q[1];
sx q[1];
rz(1.7147725) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4663812) q[0];
sx q[0];
rz(-1.2697518) q[0];
sx q[0];
rz(-1.3424323) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0703697) q[2];
sx q[2];
rz(-0.67659934) q[2];
sx q[2];
rz(0.22132193) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0328889) q[1];
sx q[1];
rz(-1.3160909) q[1];
sx q[1];
rz(0.6309066) q[1];
x q[2];
rz(1.5094425) q[3];
sx q[3];
rz(-0.7153782) q[3];
sx q[3];
rz(-0.74980199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1398805) q[2];
sx q[2];
rz(-2.6931245) q[2];
sx q[2];
rz(1.9072388) q[2];
rz(0.52302304) q[3];
sx q[3];
rz(-2.0165899) q[3];
sx q[3];
rz(-0.6559059) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6958385) q[0];
sx q[0];
rz(-3.1390751) q[0];
sx q[0];
rz(2.5939831) q[0];
rz(-2.2721263) q[1];
sx q[1];
rz(-2.165386) q[1];
sx q[1];
rz(-1.4617823) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47785178) q[0];
sx q[0];
rz(-1.2292851) q[0];
sx q[0];
rz(0.44490258) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76218729) q[2];
sx q[2];
rz(-1.7584287) q[2];
sx q[2];
rz(2.4248276) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1717277) q[1];
sx q[1];
rz(-2.0046003) q[1];
sx q[1];
rz(2.8126905) q[1];
x q[2];
rz(-0.84949228) q[3];
sx q[3];
rz(-2.2882526) q[3];
sx q[3];
rz(2.3254553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.72472858) q[2];
sx q[2];
rz(-2.8432196) q[2];
sx q[2];
rz(2.5073063) q[2];
rz(-2.0554845) q[3];
sx q[3];
rz(-0.48091286) q[3];
sx q[3];
rz(-2.0079131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.61741614) q[0];
sx q[0];
rz(-2.9845181) q[0];
sx q[0];
rz(-1.6198535) q[0];
rz(-2.1583083) q[1];
sx q[1];
rz(-2.0359928) q[1];
sx q[1];
rz(0.081667893) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0591563) q[0];
sx q[0];
rz(-2.5724223) q[0];
sx q[0];
rz(-1.3630483) q[0];
rz(-pi) q[1];
rz(-1.5310202) q[2];
sx q[2];
rz(-2.808254) q[2];
sx q[2];
rz(-0.014929742) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.52503) q[1];
sx q[1];
rz(-1.9989415) q[1];
sx q[1];
rz(-1.031428) q[1];
rz(-1.6713167) q[3];
sx q[3];
rz(-0.49197061) q[3];
sx q[3];
rz(-0.63382404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6433158) q[2];
sx q[2];
rz(-1.4772819) q[2];
sx q[2];
rz(0.91900438) q[2];
rz(0.32215858) q[3];
sx q[3];
rz(-1.5124269) q[3];
sx q[3];
rz(-0.34644103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0568986) q[0];
sx q[0];
rz(-2.0094805) q[0];
sx q[0];
rz(1.3380916) q[0];
rz(-0.40254205) q[1];
sx q[1];
rz(-1.2502547) q[1];
sx q[1];
rz(1.58443) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76257574) q[0];
sx q[0];
rz(-2.077353) q[0];
sx q[0];
rz(2.8357154) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9137325) q[2];
sx q[2];
rz(-2.3342685) q[2];
sx q[2];
rz(-2.9510483) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30774035) q[1];
sx q[1];
rz(-0.55451143) q[1];
sx q[1];
rz(-2.8608198) q[1];
x q[2];
rz(-2.040287) q[3];
sx q[3];
rz(-1.8113422) q[3];
sx q[3];
rz(2.8194373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0255787) q[2];
sx q[2];
rz(-1.2776926) q[2];
sx q[2];
rz(2.4063827) q[2];
rz(-1.0544581) q[3];
sx q[3];
rz(-2.1078883) q[3];
sx q[3];
rz(-1.2247156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16578199) q[0];
sx q[0];
rz(-2.7175792) q[0];
sx q[0];
rz(1.2531511) q[0];
rz(-1.5755298) q[1];
sx q[1];
rz(-2.0063446) q[1];
sx q[1];
rz(0.53322405) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3778702) q[0];
sx q[0];
rz(-1.3354175) q[0];
sx q[0];
rz(-0.81416582) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8715233) q[2];
sx q[2];
rz(-2.6143619) q[2];
sx q[2];
rz(-0.6747077) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.27665177) q[1];
sx q[1];
rz(-1.7930231) q[1];
sx q[1];
rz(-1.7634311) q[1];
rz(-pi) q[2];
rz(-2.5765057) q[3];
sx q[3];
rz(-1.0990407) q[3];
sx q[3];
rz(3.0155172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6553361) q[2];
sx q[2];
rz(-0.37898263) q[2];
sx q[2];
rz(-2.3113225) q[2];
rz(1.8298979) q[3];
sx q[3];
rz(-1.0887086) q[3];
sx q[3];
rz(2.8970498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61520064) q[0];
sx q[0];
rz(-1.4854234) q[0];
sx q[0];
rz(2.3080589) q[0];
rz(2.1444881) q[1];
sx q[1];
rz(-1.3778967) q[1];
sx q[1];
rz(1.5514143) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98956052) q[0];
sx q[0];
rz(-0.44267926) q[0];
sx q[0];
rz(0.17501207) q[0];
rz(-pi) q[1];
rz(0.35153206) q[2];
sx q[2];
rz(-1.0081292) q[2];
sx q[2];
rz(1.2219791) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4287069) q[1];
sx q[1];
rz(-0.16413153) q[1];
sx q[1];
rz(-2.3219789) q[1];
rz(-pi) q[2];
rz(-1.8206014) q[3];
sx q[3];
rz(-2.0916998) q[3];
sx q[3];
rz(-1.5438207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3414063) q[2];
sx q[2];
rz(-2.1087346) q[2];
sx q[2];
rz(1.7436854) q[2];
rz(0.37311113) q[3];
sx q[3];
rz(-0.92445508) q[3];
sx q[3];
rz(2.1227409) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7031192) q[0];
sx q[0];
rz(-1.4716453) q[0];
sx q[0];
rz(0.31914172) q[0];
rz(-1.9769662) q[1];
sx q[1];
rz(-1.9705557) q[1];
sx q[1];
rz(-3.107531) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24821407) q[0];
sx q[0];
rz(-0.84009087) q[0];
sx q[0];
rz(-2.7284751) q[0];
rz(-pi) q[1];
rz(0.77000846) q[2];
sx q[2];
rz(-1.3427375) q[2];
sx q[2];
rz(0.32021013) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0553375) q[1];
sx q[1];
rz(-3.1172522) q[1];
sx q[1];
rz(1.1175977) q[1];
rz(-pi) q[2];
rz(0.56896992) q[3];
sx q[3];
rz(-0.60980699) q[3];
sx q[3];
rz(-2.7066674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.065319149) q[2];
sx q[2];
rz(-1.3955782) q[2];
sx q[2];
rz(0.070153959) q[2];
rz(1.8326727) q[3];
sx q[3];
rz(-2.2891243) q[3];
sx q[3];
rz(-0.13516983) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081130505) q[0];
sx q[0];
rz(-1.9996996) q[0];
sx q[0];
rz(-0.4749701) q[0];
rz(-1.181107) q[1];
sx q[1];
rz(-2.2468552) q[1];
sx q[1];
rz(-2.2244577) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7118391) q[0];
sx q[0];
rz(-1.8022402) q[0];
sx q[0];
rz(-2.5654456) q[0];
rz(-0.31355942) q[2];
sx q[2];
rz(-0.90679996) q[2];
sx q[2];
rz(-1.0194376) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1247544) q[1];
sx q[1];
rz(-1.4887344) q[1];
sx q[1];
rz(1.7793307) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1876418) q[3];
sx q[3];
rz(-1.0055379) q[3];
sx q[3];
rz(-1.7704028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.26943031) q[2];
sx q[2];
rz(-1.7717125) q[2];
sx q[2];
rz(1.5577462) q[2];
rz(-2.9861084) q[3];
sx q[3];
rz(-1.0289861) q[3];
sx q[3];
rz(-0.99229971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6601722) q[0];
sx q[0];
rz(-1.476113) q[0];
sx q[0];
rz(0.36556622) q[0];
rz(-1.9407326) q[1];
sx q[1];
rz(-1.5301842) q[1];
sx q[1];
rz(-1.0467122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6402046) q[0];
sx q[0];
rz(-1.0602925) q[0];
sx q[0];
rz(-2.1888613) q[0];
x q[1];
rz(-1.0253941) q[2];
sx q[2];
rz(-1.1869245) q[2];
sx q[2];
rz(-1.0486368) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9871708) q[1];
sx q[1];
rz(-0.94564309) q[1];
sx q[1];
rz(2.465017) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0250859) q[3];
sx q[3];
rz(-2.0841597) q[3];
sx q[3];
rz(1.3364163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6863579) q[2];
sx q[2];
rz(-1.1150259) q[2];
sx q[2];
rz(-2.4380747) q[2];
rz(-0.55758682) q[3];
sx q[3];
rz(-0.78123868) q[3];
sx q[3];
rz(-0.47344661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2932417) q[0];
sx q[0];
rz(-2.1108755) q[0];
sx q[0];
rz(-2.1847771) q[0];
rz(-2.1059857) q[1];
sx q[1];
rz(-0.57301141) q[1];
sx q[1];
rz(1.116629) q[1];
rz(2.507092) q[2];
sx q[2];
rz(-0.77526966) q[2];
sx q[2];
rz(-2.1645968) q[2];
rz(1.612929) q[3];
sx q[3];
rz(-1.7880472) q[3];
sx q[3];
rz(0.40732297) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
