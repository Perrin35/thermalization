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
rz(-2.6211205) q[0];
sx q[0];
rz(-1.0360798) q[0];
sx q[0];
rz(0.64089027) q[0];
rz(2.9042397) q[1];
sx q[1];
rz(-0.64164716) q[1];
sx q[1];
rz(0.084903804) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7215243) q[0];
sx q[0];
rz(-0.86948538) q[0];
sx q[0];
rz(-0.12749705) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0831111) q[2];
sx q[2];
rz(-1.6019434) q[2];
sx q[2];
rz(2.8498788) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7218923) q[1];
sx q[1];
rz(-1.2806935) q[1];
sx q[1];
rz(-1.4516524) q[1];
rz(-pi) q[2];
rz(-2.1259948) q[3];
sx q[3];
rz(-2.4732504) q[3];
sx q[3];
rz(2.4618142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.25195965) q[2];
sx q[2];
rz(-1.9388988) q[2];
sx q[2];
rz(-1.6331875) q[2];
rz(0.83186045) q[3];
sx q[3];
rz(-0.32604495) q[3];
sx q[3];
rz(-1.3445725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2107596) q[0];
sx q[0];
rz(-0.63527125) q[0];
sx q[0];
rz(-2.8765836) q[0];
rz(2.3459072) q[1];
sx q[1];
rz(-2.7964451) q[1];
sx q[1];
rz(1.4211242) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20939556) q[0];
sx q[0];
rz(-1.3065814) q[0];
sx q[0];
rz(-2.2351859) q[0];
x q[1];
rz(-0.9515597) q[2];
sx q[2];
rz(-2.4251221) q[2];
sx q[2];
rz(1.8839415) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4231734) q[1];
sx q[1];
rz(-1.8224026) q[1];
sx q[1];
rz(3.0965641) q[1];
rz(3.1116074) q[3];
sx q[3];
rz(-1.4032149) q[3];
sx q[3];
rz(0.94887892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3694094) q[2];
sx q[2];
rz(-2.1464244) q[2];
sx q[2];
rz(0.65917242) q[2];
rz(2.2199421) q[3];
sx q[3];
rz(-2.0989336) q[3];
sx q[3];
rz(-1.5590182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7121861) q[0];
sx q[0];
rz(-2.5262008) q[0];
sx q[0];
rz(1.2834826) q[0];
rz(0.42690024) q[1];
sx q[1];
rz(-2.5476397) q[1];
sx q[1];
rz(-1.1675507) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10159552) q[0];
sx q[0];
rz(-1.9009703) q[0];
sx q[0];
rz(-1.9558681) q[0];
rz(-pi) q[1];
rz(-2.6111772) q[2];
sx q[2];
rz(-1.1258719) q[2];
sx q[2];
rz(0.67131587) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.42347149) q[1];
sx q[1];
rz(-1.2813066) q[1];
sx q[1];
rz(-1.4726588) q[1];
rz(-2.5646788) q[3];
sx q[3];
rz(-0.31539279) q[3];
sx q[3];
rz(-1.9384991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.20649642) q[2];
sx q[2];
rz(-1.4859716) q[2];
sx q[2];
rz(3.0237954) q[2];
rz(-1.3055034) q[3];
sx q[3];
rz(-0.85956231) q[3];
sx q[3];
rz(-0.99087244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8369668) q[0];
sx q[0];
rz(-1.0517629) q[0];
sx q[0];
rz(3.1090609) q[0];
rz(0.98371983) q[1];
sx q[1];
rz(-2.3301221) q[1];
sx q[1];
rz(0.54289114) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1586693) q[0];
sx q[0];
rz(-2.4004333) q[0];
sx q[0];
rz(-0.41038402) q[0];
rz(-pi) q[1];
rz(3.0383598) q[2];
sx q[2];
rz(-1.2599753) q[2];
sx q[2];
rz(-2.9048267) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2192083) q[1];
sx q[1];
rz(-2.3047949) q[1];
sx q[1];
rz(-2.7964785) q[1];
rz(2.3092977) q[3];
sx q[3];
rz(-1.6820388) q[3];
sx q[3];
rz(-2.5865366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5797609) q[2];
sx q[2];
rz(-1.1299955) q[2];
sx q[2];
rz(2.180991) q[2];
rz(-0.97551712) q[3];
sx q[3];
rz(-0.87149182) q[3];
sx q[3];
rz(-2.6196041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(1.8303216) q[0];
sx q[0];
rz(-0.68205849) q[0];
sx q[0];
rz(-2.8090546) q[0];
rz(2.5034816) q[1];
sx q[1];
rz(-0.6256012) q[1];
sx q[1];
rz(0.45613751) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79955937) q[0];
sx q[0];
rz(-1.4385975) q[0];
sx q[0];
rz(3.0782736) q[0];
rz(-0.78807414) q[2];
sx q[2];
rz(-1.9219134) q[2];
sx q[2];
rz(-1.1937564) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.27425452) q[1];
sx q[1];
rz(-0.98877866) q[1];
sx q[1];
rz(0.10669218) q[1];
rz(0.92335574) q[3];
sx q[3];
rz(-1.8643987) q[3];
sx q[3];
rz(1.3577611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0786324) q[2];
sx q[2];
rz(-2.3802064) q[2];
sx q[2];
rz(0.45219335) q[2];
rz(0.26442987) q[3];
sx q[3];
rz(-2.9954438) q[3];
sx q[3];
rz(-1.2368081) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0662769) q[0];
sx q[0];
rz(-2.3105268) q[0];
sx q[0];
rz(2.4237295) q[0];
rz(-1.5526937) q[1];
sx q[1];
rz(-1.6009067) q[1];
sx q[1];
rz(2.9343361) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9167647) q[0];
sx q[0];
rz(-1.5434573) q[0];
sx q[0];
rz(-0.44204373) q[0];
rz(1.4137832) q[2];
sx q[2];
rz(-2.7141389) q[2];
sx q[2];
rz(0.70892109) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0771069) q[1];
sx q[1];
rz(-2.047309) q[1];
sx q[1];
rz(0.48397343) q[1];
rz(-pi) q[2];
rz(-0.10972326) q[3];
sx q[3];
rz(-1.8302813) q[3];
sx q[3];
rz(-0.55486996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3805716) q[2];
sx q[2];
rz(-2.0863159) q[2];
sx q[2];
rz(-2.5804248) q[2];
rz(1.0593972) q[3];
sx q[3];
rz(-1.1166409) q[3];
sx q[3];
rz(1.4580956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2137432) q[0];
sx q[0];
rz(-2.1708467) q[0];
sx q[0];
rz(0.84681502) q[0];
rz(0.98833409) q[1];
sx q[1];
rz(-1.3761995) q[1];
sx q[1];
rz(2.3086937) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0578787) q[0];
sx q[0];
rz(-1.7687609) q[0];
sx q[0];
rz(-1.4651056) q[0];
rz(-0.5055228) q[2];
sx q[2];
rz(-0.16211432) q[2];
sx q[2];
rz(-2.7270728) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8036269) q[1];
sx q[1];
rz(-2.4110893) q[1];
sx q[1];
rz(2.1419163) q[1];
rz(-pi) q[2];
rz(-1.4571413) q[3];
sx q[3];
rz(-1.6196005) q[3];
sx q[3];
rz(-1.5391853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6758468) q[2];
sx q[2];
rz(-2.4387359) q[2];
sx q[2];
rz(-0.77009002) q[2];
rz(3.0105524) q[3];
sx q[3];
rz(-1.4821056) q[3];
sx q[3];
rz(-1.8946764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34762621) q[0];
sx q[0];
rz(-0.93526953) q[0];
sx q[0];
rz(2.6499709) q[0];
rz(2.8288016) q[1];
sx q[1];
rz(-2.8520695) q[1];
sx q[1];
rz(3.0591931) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9461878) q[0];
sx q[0];
rz(-3.051149) q[0];
sx q[0];
rz(-1.7121332) q[0];
rz(-2.267258) q[2];
sx q[2];
rz(-2.0249686) q[2];
sx q[2];
rz(-0.57190013) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.56697733) q[1];
sx q[1];
rz(-1.4094556) q[1];
sx q[1];
rz(1.7477504) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8776603) q[3];
sx q[3];
rz(-2.8147449) q[3];
sx q[3];
rz(-0.89034789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4440492) q[2];
sx q[2];
rz(-1.9335582) q[2];
sx q[2];
rz(0.018012878) q[2];
rz(-0.87120122) q[3];
sx q[3];
rz(-0.33659354) q[3];
sx q[3];
rz(2.474031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1910601) q[0];
sx q[0];
rz(-0.55452269) q[0];
sx q[0];
rz(2.6771123) q[0];
rz(-2.65061) q[1];
sx q[1];
rz(-1.4717439) q[1];
sx q[1];
rz(1.6741265) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53290908) q[0];
sx q[0];
rz(-2.1957046) q[0];
sx q[0];
rz(1.5708357) q[0];
x q[1];
rz(-1.2651132) q[2];
sx q[2];
rz(-1.0846429) q[2];
sx q[2];
rz(0.4403688) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.58337823) q[1];
sx q[1];
rz(-1.631307) q[1];
sx q[1];
rz(-1.8192181) q[1];
x q[2];
rz(2.2676226) q[3];
sx q[3];
rz(-2.2010816) q[3];
sx q[3];
rz(-1.1606248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1269425) q[2];
sx q[2];
rz(-1.6501959) q[2];
sx q[2];
rz(-1.3313741) q[2];
rz(-2.2018382) q[3];
sx q[3];
rz(-1.0569812) q[3];
sx q[3];
rz(-1.1872928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9068271) q[0];
sx q[0];
rz(-0.48521388) q[0];
sx q[0];
rz(1.0377129) q[0];
rz(1.0543793) q[1];
sx q[1];
rz(-1.3155921) q[1];
sx q[1];
rz(1.0494999) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3864256) q[0];
sx q[0];
rz(-1.5706129) q[0];
sx q[0];
rz(1.5210694) q[0];
x q[1];
rz(-2.0590354) q[2];
sx q[2];
rz(-1.0647213) q[2];
sx q[2];
rz(-2.807694) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29925525) q[1];
sx q[1];
rz(-1.6016869) q[1];
sx q[1];
rz(-2.1717291) q[1];
x q[2];
rz(2.353999) q[3];
sx q[3];
rz(-0.75852048) q[3];
sx q[3];
rz(-0.99419981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.99027571) q[2];
sx q[2];
rz(-1.8345202) q[2];
sx q[2];
rz(-1.446373) q[2];
rz(2.9826048) q[3];
sx q[3];
rz(-1.1587326) q[3];
sx q[3];
rz(2.3915763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1405519) q[0];
sx q[0];
rz(-0.97698553) q[0];
sx q[0];
rz(-2.3401596) q[0];
rz(-1.6544381) q[1];
sx q[1];
rz(-2.9947037) q[1];
sx q[1];
rz(-0.53302232) q[1];
rz(-2.2904446) q[2];
sx q[2];
rz(-1.0590886) q[2];
sx q[2];
rz(1.5589489) q[2];
rz(0.63538649) q[3];
sx q[3];
rz(-1.0839331) q[3];
sx q[3];
rz(-1.4270368) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
