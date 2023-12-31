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
rz(1.1336741) q[0];
sx q[0];
rz(7.8757474) q[0];
rz(-1.449861) q[1];
sx q[1];
rz(-2.4843042) q[1];
sx q[1];
rz(2.5974098) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.777594) q[0];
sx q[0];
rz(-1.6382123) q[0];
sx q[0];
rz(-1.6075587) q[0];
rz(2.7317023) q[2];
sx q[2];
rz(-0.11046834) q[2];
sx q[2];
rz(-0.2875178) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
x q[2];
rz(-2.5892026) q[3];
sx q[3];
rz(-0.033621764) q[3];
sx q[3];
rz(2.606483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.071775285) q[2];
sx q[2];
rz(-1.8775619) q[2];
sx q[2];
rz(-1.3624181) q[2];
rz(-0.028256265) q[3];
sx q[3];
rz(-1.7621721) q[3];
sx q[3];
rz(-0.79022592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-1.9702966) q[0];
rz(-0.21121875) q[1];
sx q[1];
rz(-2.6995081) q[1];
sx q[1];
rz(-1.7374977) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0728714) q[0];
sx q[0];
rz(-2.255548) q[0];
sx q[0];
rz(-2.1268334) q[0];
rz(-2.6589083) q[2];
sx q[2];
rz(-1.4233372) q[2];
sx q[2];
rz(0.54373103) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1120226) q[1];
sx q[1];
rz(-0.28407541) q[1];
sx q[1];
rz(-0.082138852) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6614463) q[3];
sx q[3];
rz(-1.4501791) q[3];
sx q[3];
rz(0.13531019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8319548) q[2];
sx q[2];
rz(-1.6654623) q[2];
sx q[2];
rz(2.1976166) q[2];
rz(-2.5850463) q[3];
sx q[3];
rz(-2.6440933) q[3];
sx q[3];
rz(-1.8054307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8495162) q[0];
sx q[0];
rz(-0.76382604) q[0];
sx q[0];
rz(1.4235494) q[0];
rz(-0.81958333) q[1];
sx q[1];
rz(-0.81033605) q[1];
sx q[1];
rz(-0.56366411) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99479988) q[0];
sx q[0];
rz(-2.1437862) q[0];
sx q[0];
rz(-1.6550199) q[0];
rz(-pi) q[1];
rz(2.5416449) q[2];
sx q[2];
rz(-1.7581345) q[2];
sx q[2];
rz(-1.3647321) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1395531) q[1];
sx q[1];
rz(-2.4783652) q[1];
sx q[1];
rz(-1.0242277) q[1];
rz(-1.1590957) q[3];
sx q[3];
rz(-2.1512096) q[3];
sx q[3];
rz(1.2702277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6391969) q[2];
sx q[2];
rz(-1.1221308) q[2];
sx q[2];
rz(-2.0813265) q[2];
rz(-3.0660196) q[3];
sx q[3];
rz(-0.85791701) q[3];
sx q[3];
rz(-3.0269572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8183427) q[0];
sx q[0];
rz(-0.30775726) q[0];
sx q[0];
rz(-2.9484205) q[0];
rz(1.5441783) q[1];
sx q[1];
rz(-1.1491821) q[1];
sx q[1];
rz(0.65778041) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13947978) q[0];
sx q[0];
rz(-1.6008899) q[0];
sx q[0];
rz(2.9667903) q[0];
rz(0.9592922) q[2];
sx q[2];
rz(-2.6311953) q[2];
sx q[2];
rz(1.4962326) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2967589) q[1];
sx q[1];
rz(-1.3169603) q[1];
sx q[1];
rz(1.0253419) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2781906) q[3];
sx q[3];
rz(-1.3232104) q[3];
sx q[3];
rz(0.86865559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0956991) q[2];
sx q[2];
rz(-0.4102439) q[2];
sx q[2];
rz(1.0470541) q[2];
rz(-2.0783157) q[3];
sx q[3];
rz(-1.1626817) q[3];
sx q[3];
rz(1.261238) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3445774) q[0];
sx q[0];
rz(-0.319096) q[0];
sx q[0];
rz(2.1176594) q[0];
rz(-0.78760415) q[1];
sx q[1];
rz(-1.5757631) q[1];
sx q[1];
rz(-2.5147298) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89676566) q[0];
sx q[0];
rz(-1.5875495) q[0];
sx q[0];
rz(-1.5411975) q[0];
x q[1];
rz(0.75682326) q[2];
sx q[2];
rz(-2.7177817) q[2];
sx q[2];
rz(2.8611081) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.53351952) q[1];
sx q[1];
rz(-2.7058209) q[1];
sx q[1];
rz(-0.79969745) q[1];
rz(-pi) q[2];
rz(-1.1781138) q[3];
sx q[3];
rz(-1.4564449) q[3];
sx q[3];
rz(-1.3254195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2287067) q[2];
sx q[2];
rz(-1.9084946) q[2];
sx q[2];
rz(2.1234925) q[2];
rz(2.5300238) q[3];
sx q[3];
rz(-1.8194149) q[3];
sx q[3];
rz(0.95156804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.4488895) q[0];
sx q[0];
rz(-1.792181) q[0];
sx q[0];
rz(1.1992136) q[0];
rz(1.1692283) q[1];
sx q[1];
rz(-2.0904082) q[1];
sx q[1];
rz(0.32454023) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6231461) q[0];
sx q[0];
rz(-1.070797) q[0];
sx q[0];
rz(0.48796939) q[0];
rz(-pi) q[1];
rz(-0.65002273) q[2];
sx q[2];
rz(-2.2215002) q[2];
sx q[2];
rz(-1.0879773) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3970916) q[1];
sx q[1];
rz(-2.607064) q[1];
sx q[1];
rz(1.6194653) q[1];
rz(-pi) q[2];
rz(1.8683897) q[3];
sx q[3];
rz(-1.9801567) q[3];
sx q[3];
rz(0.62664947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.67363182) q[2];
sx q[2];
rz(-1.8362074) q[2];
sx q[2];
rz(1.2754053) q[2];
rz(-1.0547137) q[3];
sx q[3];
rz(-2.349699) q[3];
sx q[3];
rz(1.5462497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.6773029) q[0];
sx q[0];
rz(-1.3339366) q[0];
sx q[0];
rz(1.7328847) q[0];
rz(-2.6858221) q[1];
sx q[1];
rz(-0.20142889) q[1];
sx q[1];
rz(-1.2021525) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5802655) q[0];
sx q[0];
rz(-0.94764493) q[0];
sx q[0];
rz(-0.49669637) q[0];
rz(-pi) q[1];
rz(2.7289594) q[2];
sx q[2];
rz(-2.3473783) q[2];
sx q[2];
rz(-2.7318294) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.68820885) q[1];
sx q[1];
rz(-1.6628478) q[1];
sx q[1];
rz(-1.6545014) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60336242) q[3];
sx q[3];
rz(-2.8238736) q[3];
sx q[3];
rz(0.41894223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.000164) q[2];
sx q[2];
rz(-1.8464073) q[2];
sx q[2];
rz(2.2953575) q[2];
rz(-2.6464461) q[3];
sx q[3];
rz(-2.4954002) q[3];
sx q[3];
rz(1.600986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0705567) q[0];
sx q[0];
rz(-1.8510011) q[0];
sx q[0];
rz(-2.0284247) q[0];
rz(-2.44599) q[1];
sx q[1];
rz(-2.7516987) q[1];
sx q[1];
rz(-1.5323458) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42513645) q[0];
sx q[0];
rz(-1.3279337) q[0];
sx q[0];
rz(-0.48432414) q[0];
x q[1];
rz(2.0267018) q[2];
sx q[2];
rz(-1.1258944) q[2];
sx q[2];
rz(2.0632495) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6541877) q[1];
sx q[1];
rz(-1.1206756) q[1];
sx q[1];
rz(0.022189157) q[1];
x q[2];
rz(-0.17144449) q[3];
sx q[3];
rz(-2.532634) q[3];
sx q[3];
rz(2.3122548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9939076) q[2];
sx q[2];
rz(-2.9116178) q[2];
sx q[2];
rz(2.2873986) q[2];
rz(-1.9130075) q[3];
sx q[3];
rz(-1.5649786) q[3];
sx q[3];
rz(1.6555697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5937186) q[0];
sx q[0];
rz(-0.49867189) q[0];
sx q[0];
rz(-2.4771931) q[0];
rz(1.9539072) q[1];
sx q[1];
rz(-2.4408051) q[1];
sx q[1];
rz(1.2845576) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62038319) q[0];
sx q[0];
rz(-1.8561583) q[0];
sx q[0];
rz(0.83831212) q[0];
x q[1];
rz(-2.8753488) q[2];
sx q[2];
rz(-1.9387445) q[2];
sx q[2];
rz(3.0014696) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38332332) q[1];
sx q[1];
rz(-2.3136534) q[1];
sx q[1];
rz(2.5187056) q[1];
x q[2];
rz(-2.7515718) q[3];
sx q[3];
rz(-1.5961002) q[3];
sx q[3];
rz(0.43061531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5294042) q[2];
sx q[2];
rz(-2.2470784) q[2];
sx q[2];
rz(2.5908296) q[2];
rz(2.4380056) q[3];
sx q[3];
rz(-1.0102605) q[3];
sx q[3];
rz(1.5065058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4187014) q[0];
sx q[0];
rz(-0.35690618) q[0];
sx q[0];
rz(-3.0974467) q[0];
rz(1.6126397) q[1];
sx q[1];
rz(-1.2395369) q[1];
sx q[1];
rz(-0.77967656) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3846489) q[0];
sx q[0];
rz(-2.1063519) q[0];
sx q[0];
rz(1.4063565) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0585971) q[2];
sx q[2];
rz(-2.3239845) q[2];
sx q[2];
rz(1.6542733) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.30291468) q[1];
sx q[1];
rz(-2.089114) q[1];
sx q[1];
rz(-1.360449) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46882792) q[3];
sx q[3];
rz(-2.3033597) q[3];
sx q[3];
rz(-1.9459141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49446517) q[2];
sx q[2];
rz(-0.72138849) q[2];
sx q[2];
rz(1.54281) q[2];
rz(-2.4370082) q[3];
sx q[3];
rz(-1.5141809) q[3];
sx q[3];
rz(-3.1295479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4298532) q[0];
sx q[0];
rz(-1.5562417) q[0];
sx q[0];
rz(-0.65162311) q[0];
rz(-1.7977057) q[1];
sx q[1];
rz(-1.6603036) q[1];
sx q[1];
rz(2.4659326) q[1];
rz(0.30998183) q[2];
sx q[2];
rz(-1.1520755) q[2];
sx q[2];
rz(0.68785695) q[2];
rz(1.4383437) q[3];
sx q[3];
rz(-0.9369029) q[3];
sx q[3];
rz(2.423219) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
