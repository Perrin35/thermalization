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
rz(0.6495629) q[0];
sx q[0];
rz(-0.49512884) q[0];
sx q[0];
rz(0.25547096) q[0];
rz(2.7893692) q[1];
sx q[1];
rz(-1.398634) q[1];
sx q[1];
rz(1.4056828) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3509003) q[0];
sx q[0];
rz(-0.003482799) q[0];
sx q[0];
rz(-3.0538959) q[0];
x q[1];
rz(0.5159401) q[2];
sx q[2];
rz(-1.4778412) q[2];
sx q[2];
rz(-3.1414714) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.60599323) q[1];
sx q[1];
rz(-2.7584425) q[1];
sx q[1];
rz(-1.5616756) q[1];
rz(-pi) q[2];
rz(-3.0902946) q[3];
sx q[3];
rz(-0.56863943) q[3];
sx q[3];
rz(2.7647247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.63456717) q[2];
sx q[2];
rz(-0.031012379) q[2];
sx q[2];
rz(-2.3090889) q[2];
rz(-2.3333874) q[3];
sx q[3];
rz(-3.1263604) q[3];
sx q[3];
rz(0.25850779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38043624) q[0];
sx q[0];
rz(-0.34270898) q[0];
sx q[0];
rz(2.9535182) q[0];
rz(-0.07218083) q[1];
sx q[1];
rz(-1.0344104) q[1];
sx q[1];
rz(1.5123051) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24969026) q[0];
sx q[0];
rz(-1.238126) q[0];
sx q[0];
rz(1.692125) q[0];
rz(-pi) q[1];
x q[1];
rz(2.15835) q[2];
sx q[2];
rz(-3.0863783) q[2];
sx q[2];
rz(3.0172341) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.364671) q[1];
sx q[1];
rz(-1.6699759) q[1];
sx q[1];
rz(-1.6302376) q[1];
rz(1.0076804) q[3];
sx q[3];
rz(-0.84656871) q[3];
sx q[3];
rz(2.7359642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.170681) q[2];
sx q[2];
rz(-0.6450246) q[2];
sx q[2];
rz(1.8485273) q[2];
rz(-0.34717789) q[3];
sx q[3];
rz(-2.8721589) q[3];
sx q[3];
rz(-0.080304317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3229356) q[0];
sx q[0];
rz(-1.8523536) q[0];
sx q[0];
rz(2.6710508) q[0];
rz(1.200354) q[1];
sx q[1];
rz(-2.410694) q[1];
sx q[1];
rz(1.8847195) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6022105) q[0];
sx q[0];
rz(-2.4823144) q[0];
sx q[0];
rz(-2.1677917) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1684709) q[2];
sx q[2];
rz(-0.16301708) q[2];
sx q[2];
rz(2.4578641) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.10373579) q[1];
sx q[1];
rz(-1.5694261) q[1];
sx q[1];
rz(1.5562773) q[1];
x q[2];
rz(1.1670349) q[3];
sx q[3];
rz(-0.65126538) q[3];
sx q[3];
rz(0.854597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54258004) q[2];
sx q[2];
rz(-2.161669) q[2];
sx q[2];
rz(2.2194594) q[2];
rz(-1.7982091) q[3];
sx q[3];
rz(-0.94815367) q[3];
sx q[3];
rz(1.5141727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2682997) q[0];
sx q[0];
rz(-1.0404077) q[0];
sx q[0];
rz(-2.701395) q[0];
rz(1.6104376) q[1];
sx q[1];
rz(-1.6596158) q[1];
sx q[1];
rz(2.8940315) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2092106) q[0];
sx q[0];
rz(-2.4510151) q[0];
sx q[0];
rz(0.85777046) q[0];
rz(-pi) q[1];
rz(0.48657067) q[2];
sx q[2];
rz(-0.19127327) q[2];
sx q[2];
rz(0.74081206) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.3888445) q[1];
sx q[1];
rz(-1.6444863) q[1];
sx q[1];
rz(-1.8703413) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64856802) q[3];
sx q[3];
rz(-1.5952791) q[3];
sx q[3];
rz(-2.690993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89746785) q[2];
sx q[2];
rz(-2.1110057) q[2];
sx q[2];
rz(1.6991276) q[2];
rz(-2.9407732) q[3];
sx q[3];
rz(-1.2186058) q[3];
sx q[3];
rz(2.0012746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9852801) q[0];
sx q[0];
rz(-2.9887178) q[0];
sx q[0];
rz(-0.54541624) q[0];
rz(-3.0975869) q[1];
sx q[1];
rz(-0.017887201) q[1];
sx q[1];
rz(0.47637475) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2891083) q[0];
sx q[0];
rz(-1.327564) q[0];
sx q[0];
rz(-1.6575251) q[0];
rz(-pi) q[1];
rz(2.7353941) q[2];
sx q[2];
rz(-1.9425689) q[2];
sx q[2];
rz(-0.68454725) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0477476) q[1];
sx q[1];
rz(-0.72870164) q[1];
sx q[1];
rz(-1.2959667) q[1];
x q[2];
rz(-3.0073193) q[3];
sx q[3];
rz(-2.0760787) q[3];
sx q[3];
rz(-0.093029417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5215317) q[2];
sx q[2];
rz(-0.69053495) q[2];
sx q[2];
rz(0.36038348) q[2];
rz(2.9195869) q[3];
sx q[3];
rz(-1.6111776) q[3];
sx q[3];
rz(1.7300026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5510657) q[0];
sx q[0];
rz(-1.5532302) q[0];
sx q[0];
rz(1.5565514) q[0];
rz(0.12697728) q[1];
sx q[1];
rz(-1.304909) q[1];
sx q[1];
rz(3.051905) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2253826) q[0];
sx q[0];
rz(-1.0400335) q[0];
sx q[0];
rz(-2.5368774) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96025567) q[2];
sx q[2];
rz(-0.88443236) q[2];
sx q[2];
rz(-2.6853564) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.27881926) q[1];
sx q[1];
rz(-1.3403646) q[1];
sx q[1];
rz(2.6342175) q[1];
rz(-pi) q[2];
rz(-0.60097127) q[3];
sx q[3];
rz(-1.8716164) q[3];
sx q[3];
rz(-1.6114637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.108532) q[2];
sx q[2];
rz(-0.33691275) q[2];
sx q[2];
rz(0.25165558) q[2];
rz(-3.1015977) q[3];
sx q[3];
rz(-0.34188855) q[3];
sx q[3];
rz(2.6778636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43117487) q[0];
sx q[0];
rz(-3.0496821) q[0];
sx q[0];
rz(2.726626) q[0];
rz(1.4893432) q[1];
sx q[1];
rz(-0.057561189) q[1];
sx q[1];
rz(2.8156978) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6221507) q[0];
sx q[0];
rz(-1.2076734) q[0];
sx q[0];
rz(1.0958998) q[0];
rz(-pi) q[1];
rz(-1.4140747) q[2];
sx q[2];
rz(-0.5236519) q[2];
sx q[2];
rz(1.0207506) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.88837543) q[1];
sx q[1];
rz(-1.0723812) q[1];
sx q[1];
rz(1.8640169) q[1];
x q[2];
rz(0.075906673) q[3];
sx q[3];
rz(-1.2682639) q[3];
sx q[3];
rz(2.3940866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9140279) q[2];
sx q[2];
rz(-1.808337) q[2];
sx q[2];
rz(-1.0464767) q[2];
rz(0.21876167) q[3];
sx q[3];
rz(-1.0294139) q[3];
sx q[3];
rz(1.7441162) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47531146) q[0];
sx q[0];
rz(-0.0175716) q[0];
sx q[0];
rz(0.46324357) q[0];
rz(-2.8766368) q[1];
sx q[1];
rz(-3.1400561) q[1];
sx q[1];
rz(1.6418246) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4508598) q[0];
sx q[0];
rz(-1.5895278) q[0];
sx q[0];
rz(2.0035726) q[0];
rz(-1.6258662) q[2];
sx q[2];
rz(-1.1257505) q[2];
sx q[2];
rz(-0.17521706) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38507515) q[1];
sx q[1];
rz(-0.2236872) q[1];
sx q[1];
rz(-1.7023515) q[1];
rz(-1.8231976) q[3];
sx q[3];
rz(-1.6381255) q[3];
sx q[3];
rz(-0.52378262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.90193343) q[2];
sx q[2];
rz(-2.689211) q[2];
sx q[2];
rz(1.3182013) q[2];
rz(0.0031331172) q[3];
sx q[3];
rz(-1.2001218) q[3];
sx q[3];
rz(1.1656632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5636633) q[0];
sx q[0];
rz(-0.85275537) q[0];
sx q[0];
rz(0.71459115) q[0];
rz(0.21358061) q[1];
sx q[1];
rz(-3.112308) q[1];
sx q[1];
rz(-1.9307131) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9506257) q[0];
sx q[0];
rz(-1.4651148) q[0];
sx q[0];
rz(1.8777385) q[0];
x q[1];
rz(1.2229332) q[2];
sx q[2];
rz(-0.90556723) q[2];
sx q[2];
rz(0.72345483) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9309733) q[1];
sx q[1];
rz(-2.5323091) q[1];
sx q[1];
rz(1.0675201) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15809246) q[3];
sx q[3];
rz(-1.7580346) q[3];
sx q[3];
rz(-2.9824389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1750298) q[2];
sx q[2];
rz(-3.1201456) q[2];
sx q[2];
rz(-0.19801298) q[2];
rz(-0.35061947) q[3];
sx q[3];
rz(-1.5712761) q[3];
sx q[3];
rz(-1.143379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.7570033) q[0];
sx q[0];
rz(-0.92246711) q[0];
sx q[0];
rz(-0.61436999) q[0];
rz(3.0820097) q[1];
sx q[1];
rz(-3.0501084) q[1];
sx q[1];
rz(1.7170067) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2587268) q[0];
sx q[0];
rz(-0.75260163) q[0];
sx q[0];
rz(-1.2105788) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0103127) q[2];
sx q[2];
rz(-1.5727709) q[2];
sx q[2];
rz(1.9406089) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9944133) q[1];
sx q[1];
rz(-2.96619) q[1];
sx q[1];
rz(2.0270623) q[1];
rz(-1.9779102) q[3];
sx q[3];
rz(-1.3729248) q[3];
sx q[3];
rz(-2.1382633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0397348) q[2];
sx q[2];
rz(-0.27145824) q[2];
sx q[2];
rz(-0.52379918) q[2];
rz(-2.0959181) q[3];
sx q[3];
rz(-1.8664675) q[3];
sx q[3];
rz(2.6773793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47281784) q[0];
sx q[0];
rz(-1.5325118) q[0];
sx q[0];
rz(-1.4788628) q[0];
rz(-1.7300023) q[1];
sx q[1];
rz(-2.9099606) q[1];
sx q[1];
rz(0.06106613) q[1];
rz(0.642943) q[2];
sx q[2];
rz(-1.8131154) q[2];
sx q[2];
rz(-2.9396069) q[2];
rz(-2.6749776) q[3];
sx q[3];
rz(-1.7468921) q[3];
sx q[3];
rz(0.56260059) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
