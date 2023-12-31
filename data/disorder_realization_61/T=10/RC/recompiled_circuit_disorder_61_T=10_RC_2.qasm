OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7282038) q[0];
sx q[0];
rz(-2.0079186) q[0];
sx q[0];
rz(1.5490305) q[0];
rz(-1.449861) q[1];
sx q[1];
rz(-2.4843042) q[1];
sx q[1];
rz(-0.54418286) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13574164) q[0];
sx q[0];
rz(-0.07677456) q[0];
sx q[0];
rz(0.49850054) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0402096) q[2];
sx q[2];
rz(-1.6147436) q[2];
sx q[2];
rz(-1.6909388) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4197598) q[1];
sx q[1];
rz(-0.83336035) q[1];
sx q[1];
rz(-2.1268197) q[1];
rz(-pi) q[2];
rz(-2.5892026) q[3];
sx q[3];
rz(-3.1079709) q[3];
sx q[3];
rz(-2.606483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0698174) q[2];
sx q[2];
rz(-1.8775619) q[2];
sx q[2];
rz(1.3624181) q[2];
rz(0.028256265) q[3];
sx q[3];
rz(-1.3794206) q[3];
sx q[3];
rz(2.3513667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4966999) q[0];
sx q[0];
rz(-2.4364478) q[0];
sx q[0];
rz(1.171296) q[0];
rz(2.9303739) q[1];
sx q[1];
rz(-2.6995081) q[1];
sx q[1];
rz(-1.7374977) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8453168) q[0];
sx q[0];
rz(-0.8527841) q[0];
sx q[0];
rz(-0.57384558) q[0];
rz(-pi) q[1];
rz(0.30971576) q[2];
sx q[2];
rz(-2.6385912) q[2];
sx q[2];
rz(-2.3878857) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6200871) q[1];
sx q[1];
rz(-1.5477991) q[1];
sx q[1];
rz(-2.8584245) q[1];
x q[2];
rz(-0.25661616) q[3];
sx q[3];
rz(-0.49391541) q[3];
sx q[3];
rz(-1.6625422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8319548) q[2];
sx q[2];
rz(-1.6654623) q[2];
sx q[2];
rz(-2.1976166) q[2];
rz(0.55654636) q[3];
sx q[3];
rz(-2.6440933) q[3];
sx q[3];
rz(-1.8054307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8495162) q[0];
sx q[0];
rz(-2.3777666) q[0];
sx q[0];
rz(1.7180432) q[0];
rz(2.3220093) q[1];
sx q[1];
rz(-0.81033605) q[1];
sx q[1];
rz(-0.56366411) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1492796) q[0];
sx q[0];
rz(-0.57846071) q[0];
sx q[0];
rz(0.1296541) q[0];
rz(-pi) q[1];
rz(-1.7965505) q[2];
sx q[2];
rz(-2.1588237) q[2];
sx q[2];
rz(-0.079344582) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3445661) q[1];
sx q[1];
rz(-2.1246506) q[1];
sx q[1];
rz(-2.7558541) q[1];
rz(-0.62108576) q[3];
sx q[3];
rz(-1.9120145) q[3];
sx q[3];
rz(-0.065545557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6391969) q[2];
sx q[2];
rz(-1.1221308) q[2];
sx q[2];
rz(2.0813265) q[2];
rz(-3.0660196) q[3];
sx q[3];
rz(-0.85791701) q[3];
sx q[3];
rz(-3.0269572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8183427) q[0];
sx q[0];
rz(-0.30775726) q[0];
sx q[0];
rz(-0.19317214) q[0];
rz(1.5441783) q[1];
sx q[1];
rz(-1.1491821) q[1];
sx q[1];
rz(0.65778041) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4260027) q[0];
sx q[0];
rz(-1.7455187) q[0];
sx q[0];
rz(1.5402373) q[0];
rz(-pi) q[1];
rz(1.1409608) q[2];
sx q[2];
rz(-1.2865215) q[2];
sx q[2];
rz(0.62361275) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8448338) q[1];
sx q[1];
rz(-1.8246324) q[1];
sx q[1];
rz(1.0253419) q[1];
rz(1.863402) q[3];
sx q[3];
rz(-1.3232104) q[3];
sx q[3];
rz(2.2729371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0956991) q[2];
sx q[2];
rz(-2.7313488) q[2];
sx q[2];
rz(1.0470541) q[2];
rz(2.0783157) q[3];
sx q[3];
rz(-1.9789109) q[3];
sx q[3];
rz(-1.8803546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.7970153) q[0];
sx q[0];
rz(-0.319096) q[0];
sx q[0];
rz(2.1176594) q[0];
rz(0.78760415) q[1];
sx q[1];
rz(-1.5658295) q[1];
sx q[1];
rz(-2.5147298) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89676566) q[0];
sx q[0];
rz(-1.5875495) q[0];
sx q[0];
rz(-1.5411975) q[0];
rz(2.3847694) q[2];
sx q[2];
rz(-0.42381091) q[2];
sx q[2];
rz(-0.28048453) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31506854) q[1];
sx q[1];
rz(-1.2722004) q[1];
sx q[1];
rz(1.8930757) q[1];
rz(-3.0179126) q[3];
sx q[3];
rz(-1.1808174) q[3];
sx q[3];
rz(-0.19815138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.91288599) q[2];
sx q[2];
rz(-1.9084946) q[2];
sx q[2];
rz(-2.1234925) q[2];
rz(0.61156887) q[3];
sx q[3];
rz(-1.8194149) q[3];
sx q[3];
rz(-0.95156804) q[3];
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
rz(1.6927032) q[0];
sx q[0];
rz(-1.792181) q[0];
sx q[0];
rz(-1.1992136) q[0];
rz(-1.1692283) q[1];
sx q[1];
rz(-2.0904082) q[1];
sx q[1];
rz(2.8170524) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6231461) q[0];
sx q[0];
rz(-2.0707957) q[0];
sx q[0];
rz(2.6536233) q[0];
x q[1];
rz(-0.89914497) q[2];
sx q[2];
rz(-2.2567344) q[2];
sx q[2];
rz(1.1555954) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.68795825) q[1];
sx q[1];
rz(-2.1046241) q[1];
sx q[1];
rz(-0.028793528) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8683897) q[3];
sx q[3];
rz(-1.1614359) q[3];
sx q[3];
rz(-0.62664947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67363182) q[2];
sx q[2];
rz(-1.8362074) q[2];
sx q[2];
rz(1.2754053) q[2];
rz(1.0547137) q[3];
sx q[3];
rz(-0.79189363) q[3];
sx q[3];
rz(-1.595343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4642898) q[0];
sx q[0];
rz(-1.3339366) q[0];
sx q[0];
rz(-1.7328847) q[0];
rz(0.45577058) q[1];
sx q[1];
rz(-0.20142889) q[1];
sx q[1];
rz(-1.2021525) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3287438) q[0];
sx q[0];
rz(-2.3658731) q[0];
sx q[0];
rz(2.1562955) q[0];
rz(-pi) q[1];
rz(0.75041109) q[2];
sx q[2];
rz(-1.8609036) q[2];
sx q[2];
rz(-1.6828705) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89029965) q[1];
sx q[1];
rz(-1.6541462) q[1];
sx q[1];
rz(3.0492196) q[1];
x q[2];
rz(0.2644516) q[3];
sx q[3];
rz(-1.7489986) q[3];
sx q[3];
rz(1.7314535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.000164) q[2];
sx q[2];
rz(-1.2951853) q[2];
sx q[2];
rz(-2.2953575) q[2];
rz(2.6464461) q[3];
sx q[3];
rz(-0.64619243) q[3];
sx q[3];
rz(-1.5406066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0705567) q[0];
sx q[0];
rz(-1.2905916) q[0];
sx q[0];
rz(-2.0284247) q[0];
rz(2.44599) q[1];
sx q[1];
rz(-2.7516987) q[1];
sx q[1];
rz(1.5323458) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0198062) q[0];
sx q[0];
rz(-2.0397423) q[0];
sx q[0];
rz(1.8437587) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3959827) q[2];
sx q[2];
rz(-0.62586212) q[2];
sx q[2];
rz(1.9288043) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6032226) q[1];
sx q[1];
rz(-0.45062989) q[1];
sx q[1];
rz(1.5249114) q[1];
rz(-pi) q[2];
rz(-0.17144449) q[3];
sx q[3];
rz(-0.60895863) q[3];
sx q[3];
rz(-2.3122548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9939076) q[2];
sx q[2];
rz(-2.9116178) q[2];
sx q[2];
rz(-2.2873986) q[2];
rz(-1.2285852) q[3];
sx q[3];
rz(-1.5766141) q[3];
sx q[3];
rz(1.6555697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5937186) q[0];
sx q[0];
rz(-0.49867189) q[0];
sx q[0];
rz(-2.4771931) q[0];
rz(-1.9539072) q[1];
sx q[1];
rz(-0.70078754) q[1];
sx q[1];
rz(-1.857035) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2537456) q[0];
sx q[0];
rz(-2.3652024) q[0];
sx q[0];
rz(1.1573769) q[0];
rz(-pi) q[1];
rz(-0.2662439) q[2];
sx q[2];
rz(-1.2028482) q[2];
sx q[2];
rz(3.0014696) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7582693) q[1];
sx q[1];
rz(-0.82793923) q[1];
sx q[1];
rz(2.5187056) q[1];
x q[2];
rz(-3.0751238) q[3];
sx q[3];
rz(-2.7507938) q[3];
sx q[3];
rz(1.9399411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5294042) q[2];
sx q[2];
rz(-2.2470784) q[2];
sx q[2];
rz(2.5908296) q[2];
rz(-0.70358706) q[3];
sx q[3];
rz(-1.0102605) q[3];
sx q[3];
rz(-1.6350869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4187014) q[0];
sx q[0];
rz(-0.35690618) q[0];
sx q[0];
rz(3.0974467) q[0];
rz(-1.528953) q[1];
sx q[1];
rz(-1.2395369) q[1];
sx q[1];
rz(-0.77967656) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2432602) q[0];
sx q[0];
rz(-1.7120449) q[0];
sx q[0];
rz(2.6000573) q[0];
x q[1];
rz(-1.0829955) q[2];
sx q[2];
rz(-0.81760815) q[2];
sx q[2];
rz(-1.6542733) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.979094) q[1];
sx q[1];
rz(-1.3884123) q[1];
sx q[1];
rz(-0.52796332) q[1];
rz(-pi) q[2];
rz(1.1053106) q[3];
sx q[3];
rz(-2.2959384) q[3];
sx q[3];
rz(-0.5474962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49446517) q[2];
sx q[2];
rz(-2.4202042) q[2];
sx q[2];
rz(1.54281) q[2];
rz(-0.70458448) q[3];
sx q[3];
rz(-1.6274118) q[3];
sx q[3];
rz(-3.1295479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4298532) q[0];
sx q[0];
rz(-1.5562417) q[0];
sx q[0];
rz(-0.65162311) q[0];
rz(1.7977057) q[1];
sx q[1];
rz(-1.4812891) q[1];
sx q[1];
rz(-0.67566009) q[1];
rz(-2.1716739) q[2];
sx q[2];
rz(-0.51545943) q[2];
sx q[2];
rz(-3.1209844) q[2];
rz(0.17776168) q[3];
sx q[3];
rz(-0.64571417) q[3];
sx q[3];
rz(-0.4971102) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
