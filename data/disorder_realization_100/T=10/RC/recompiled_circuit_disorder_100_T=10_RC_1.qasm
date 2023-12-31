OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93958062) q[0];
sx q[0];
rz(2.7913845) q[0];
sx q[0];
rz(12.933001) q[0];
rz(-2.2740013) q[1];
sx q[1];
rz(-0.64414135) q[1];
sx q[1];
rz(1.6860513) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6018574) q[0];
sx q[0];
rz(-2.5849113) q[0];
sx q[0];
rz(0.88168983) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2725265) q[2];
sx q[2];
rz(-1.367374) q[2];
sx q[2];
rz(-0.8806526) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9283596) q[1];
sx q[1];
rz(-1.8395437) q[1];
sx q[1];
rz(-2.1409722) q[1];
rz(-pi) q[2];
rz(-1.5815758) q[3];
sx q[3];
rz(-1.5228776) q[3];
sx q[3];
rz(2.3673494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.036844) q[2];
sx q[2];
rz(-1.8061183) q[2];
sx q[2];
rz(-2.74995) q[2];
rz(0.13970217) q[3];
sx q[3];
rz(-0.67290664) q[3];
sx q[3];
rz(-3.1203549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7249107) q[0];
sx q[0];
rz(-1.7371438) q[0];
sx q[0];
rz(-1.2765983) q[0];
rz(2.2712767) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(-1.93719) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6635839) q[0];
sx q[0];
rz(-1.1482129) q[0];
sx q[0];
rz(-2.708486) q[0];
rz(-pi) q[1];
rz(2.2764678) q[2];
sx q[2];
rz(-1.428831) q[2];
sx q[2];
rz(2.1324468) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0871036) q[1];
sx q[1];
rz(-0.37885715) q[1];
sx q[1];
rz(-2.8701251) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2435917) q[3];
sx q[3];
rz(-1.3452531) q[3];
sx q[3];
rz(-1.1225835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5045972) q[2];
sx q[2];
rz(-2.6837139) q[2];
sx q[2];
rz(-1.9223928) q[2];
rz(2.7870264) q[3];
sx q[3];
rz(-0.48186007) q[3];
sx q[3];
rz(1.5725117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5511659) q[0];
sx q[0];
rz(-1.5274436) q[0];
sx q[0];
rz(-0.61808008) q[0];
rz(0.016050054) q[1];
sx q[1];
rz(-0.87688223) q[1];
sx q[1];
rz(1.191167) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5601657) q[0];
sx q[0];
rz(-1.7227356) q[0];
sx q[0];
rz(-1.2537969) q[0];
x q[1];
rz(2.8580655) q[2];
sx q[2];
rz(-1.8054188) q[2];
sx q[2];
rz(1.4153751) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9900073) q[1];
sx q[1];
rz(-1.3602435) q[1];
sx q[1];
rz(-1.7819808) q[1];
x q[2];
rz(-0.54179811) q[3];
sx q[3];
rz(-1.9199315) q[3];
sx q[3];
rz(-2.9475714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.96737343) q[2];
sx q[2];
rz(-0.9625532) q[2];
sx q[2];
rz(-2.1276786) q[2];
rz(-1.7845456) q[3];
sx q[3];
rz(-0.8299399) q[3];
sx q[3];
rz(-1.0323662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33092609) q[0];
sx q[0];
rz(-1.4241011) q[0];
sx q[0];
rz(-0.20430918) q[0];
rz(-1.3775685) q[1];
sx q[1];
rz(-1.7267449) q[1];
sx q[1];
rz(0.92299443) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8078634) q[0];
sx q[0];
rz(-0.45818746) q[0];
sx q[0];
rz(-2.3966167) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0976474) q[2];
sx q[2];
rz(-1.5825795) q[2];
sx q[2];
rz(0.95566434) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.877362) q[1];
sx q[1];
rz(-1.1738395) q[1];
sx q[1];
rz(3.0287659) q[1];
x q[2];
rz(-0.022425671) q[3];
sx q[3];
rz(-2.2750345) q[3];
sx q[3];
rz(-1.1426329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4541645) q[2];
sx q[2];
rz(-1.5565846) q[2];
sx q[2];
rz(-2.6320809) q[2];
rz(2.4984958) q[3];
sx q[3];
rz(-2.0431079) q[3];
sx q[3];
rz(-0.96737868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5287857) q[0];
sx q[0];
rz(-1.9354154) q[0];
sx q[0];
rz(0.91941419) q[0];
rz(2.1557504) q[1];
sx q[1];
rz(-1.7835833) q[1];
sx q[1];
rz(-1.7808419) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6960984) q[0];
sx q[0];
rz(-1.7901663) q[0];
sx q[0];
rz(-1.7646837) q[0];
x q[1];
rz(-0.947457) q[2];
sx q[2];
rz(-2.2286378) q[2];
sx q[2];
rz(-0.62589494) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.50966108) q[1];
sx q[1];
rz(-1.7778548) q[1];
sx q[1];
rz(-2.6775371) q[1];
rz(-0.97134437) q[3];
sx q[3];
rz(-1.6347307) q[3];
sx q[3];
rz(-1.5435227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.78648606) q[2];
sx q[2];
rz(-1.7559218) q[2];
sx q[2];
rz(0.11486593) q[2];
rz(-2.6857175) q[3];
sx q[3];
rz(-1.3331648) q[3];
sx q[3];
rz(-1.8615287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.8961261) q[0];
sx q[0];
rz(-1.2446612) q[0];
sx q[0];
rz(-1.2448357) q[0];
rz(-0.70760977) q[1];
sx q[1];
rz(-1.5039624) q[1];
sx q[1];
rz(-0.27522603) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9305785) q[0];
sx q[0];
rz(-1.0761346) q[0];
sx q[0];
rz(2.7633694) q[0];
rz(-2.6456344) q[2];
sx q[2];
rz(-0.99018103) q[2];
sx q[2];
rz(0.95326391) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.59906193) q[1];
sx q[1];
rz(-1.2236569) q[1];
sx q[1];
rz(2.6229726) q[1];
x q[2];
rz(-1.4554548) q[3];
sx q[3];
rz(-2.5177258) q[3];
sx q[3];
rz(-1.3184402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.747921) q[2];
sx q[2];
rz(-1.5966406) q[2];
sx q[2];
rz(0.45864027) q[2];
rz(2.4300872) q[3];
sx q[3];
rz(-0.68370521) q[3];
sx q[3];
rz(-0.59035629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28073072) q[0];
sx q[0];
rz(-1.1870528) q[0];
sx q[0];
rz(1.77805) q[0];
rz(1.4200312) q[1];
sx q[1];
rz(-0.51912156) q[1];
sx q[1];
rz(-2.4218959) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62455432) q[0];
sx q[0];
rz(-1.8554243) q[0];
sx q[0];
rz(2.4863003) q[0];
rz(-pi) q[1];
rz(0.85992203) q[2];
sx q[2];
rz(-2.2051174) q[2];
sx q[2];
rz(0.89564043) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.67932018) q[1];
sx q[1];
rz(-1.6151036) q[1];
sx q[1];
rz(0.54507749) q[1];
rz(1.3004488) q[3];
sx q[3];
rz(-2.2336604) q[3];
sx q[3];
rz(-2.8519252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.65934962) q[2];
sx q[2];
rz(-1.6094001) q[2];
sx q[2];
rz(2.4534295) q[2];
rz(0.041953772) q[3];
sx q[3];
rz(-2.0418906) q[3];
sx q[3];
rz(0.28013128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2280837) q[0];
sx q[0];
rz(-1.9788454) q[0];
sx q[0];
rz(-2.9550609) q[0];
rz(0.60449156) q[1];
sx q[1];
rz(-2.1291321) q[1];
sx q[1];
rz(1.4950745) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.625531) q[0];
sx q[0];
rz(-0.4058668) q[0];
sx q[0];
rz(2.5748475) q[0];
rz(-pi) q[1];
rz(-2.4453156) q[2];
sx q[2];
rz(-1.7211203) q[2];
sx q[2];
rz(2.6339649) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5563158) q[1];
sx q[1];
rz(-1.8299412) q[1];
sx q[1];
rz(-1.6782594) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80963482) q[3];
sx q[3];
rz(-2.0658884) q[3];
sx q[3];
rz(1.9593057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.24215332) q[2];
sx q[2];
rz(-2.184325) q[2];
sx q[2];
rz(3.126826) q[2];
rz(-1.8642558) q[3];
sx q[3];
rz(-1.8321962) q[3];
sx q[3];
rz(-1.7285255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9799141) q[0];
sx q[0];
rz(-2.4665687) q[0];
sx q[0];
rz(-2.263608) q[0];
rz(2.9453078) q[1];
sx q[1];
rz(-1.9613962) q[1];
sx q[1];
rz(1.7810129) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44952794) q[0];
sx q[0];
rz(-0.99943107) q[0];
sx q[0];
rz(0.044762386) q[0];
x q[1];
rz(-0.046710308) q[2];
sx q[2];
rz(-2.197406) q[2];
sx q[2];
rz(-1.7358629) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1189551) q[1];
sx q[1];
rz(-1.2345018) q[1];
sx q[1];
rz(-0.35204661) q[1];
x q[2];
rz(-0.15501539) q[3];
sx q[3];
rz(-1.1559556) q[3];
sx q[3];
rz(-2.993194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.517841) q[2];
sx q[2];
rz(-1.5072284) q[2];
sx q[2];
rz(0.8927792) q[2];
rz(3.1023074) q[3];
sx q[3];
rz(-1.4792484) q[3];
sx q[3];
rz(2.3915496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84353012) q[0];
sx q[0];
rz(-0.003304464) q[0];
sx q[0];
rz(-0.046534006) q[0];
rz(-2.2325366) q[1];
sx q[1];
rz(-2.2687056) q[1];
sx q[1];
rz(0.7199026) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90912949) q[0];
sx q[0];
rz(-1.6502408) q[0];
sx q[0];
rz(-1.3264873) q[0];
x q[1];
rz(3.1245329) q[2];
sx q[2];
rz(-0.3088726) q[2];
sx q[2];
rz(-0.19933867) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1788927) q[1];
sx q[1];
rz(-1.3239412) q[1];
sx q[1];
rz(0.015951338) q[1];
rz(-0.73086892) q[3];
sx q[3];
rz(-2.1248397) q[3];
sx q[3];
rz(2.4027367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7211192) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(-0.026467888) q[2];
rz(-1.8376393) q[3];
sx q[3];
rz(-2.3243258) q[3];
sx q[3];
rz(1.2560237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7883041) q[0];
sx q[0];
rz(-1.7699387) q[0];
sx q[0];
rz(1.8040245) q[0];
rz(2.9251255) q[1];
sx q[1];
rz(-1.6995866) q[1];
sx q[1];
rz(1.235984) q[1];
rz(-0.70710612) q[2];
sx q[2];
rz(-1.8137992) q[2];
sx q[2];
rz(0.9546311) q[2];
rz(-0.032565928) q[3];
sx q[3];
rz(-0.98441549) q[3];
sx q[3];
rz(2.050642) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
