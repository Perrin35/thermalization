OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4274347) q[0];
sx q[0];
rz(-0.43570575) q[0];
sx q[0];
rz(0.92619196) q[0];
rz(1.9594833) q[1];
sx q[1];
rz(-0.73298454) q[1];
sx q[1];
rz(0.37252537) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6719138) q[0];
sx q[0];
rz(-0.9979453) q[0];
sx q[0];
rz(-1.8125305) q[0];
rz(-2.9973642) q[2];
sx q[2];
rz(-0.63604522) q[2];
sx q[2];
rz(-0.47338212) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8718611) q[1];
sx q[1];
rz(-1.4565804) q[1];
sx q[1];
rz(-1.7608587) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11301179) q[3];
sx q[3];
rz(-1.3556619) q[3];
sx q[3];
rz(0.61258951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7131876) q[2];
sx q[2];
rz(-1.5593854) q[2];
sx q[2];
rz(-0.92450809) q[2];
rz(-1.6690147) q[3];
sx q[3];
rz(-1.2481097) q[3];
sx q[3];
rz(1.6424461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9532042) q[0];
sx q[0];
rz(-0.062347118) q[0];
sx q[0];
rz(1.6054608) q[0];
rz(2.9470782) q[1];
sx q[1];
rz(-1.8201927) q[1];
sx q[1];
rz(-3.0867192) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40614906) q[0];
sx q[0];
rz(-2.9950954) q[0];
sx q[0];
rz(-2.2551401) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5218711) q[2];
sx q[2];
rz(-2.0267068) q[2];
sx q[2];
rz(-2.1330619) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.730719) q[1];
sx q[1];
rz(-1.5748071) q[1];
sx q[1];
rz(1.2282759) q[1];
x q[2];
rz(-0.67316405) q[3];
sx q[3];
rz(-1.737397) q[3];
sx q[3];
rz(1.3413606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.43869552) q[2];
sx q[2];
rz(-0.56151152) q[2];
sx q[2];
rz(1.6595718) q[2];
rz(-2.1510018) q[3];
sx q[3];
rz(-1.4086658) q[3];
sx q[3];
rz(-1.3668485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
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
rz(2.7822587) q[0];
sx q[0];
rz(-2.6222836) q[0];
sx q[0];
rz(-2.7666132) q[0];
rz(-2.8938876) q[1];
sx q[1];
rz(-1.9539555) q[1];
sx q[1];
rz(-1.3365655) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1776745) q[0];
sx q[0];
rz(-1.5128711) q[0];
sx q[0];
rz(-1.792568) q[0];
x q[1];
rz(-2.5518718) q[2];
sx q[2];
rz(-1.1955185) q[2];
sx q[2];
rz(-2.7207665) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.8001576) q[1];
sx q[1];
rz(-2.2812124) q[1];
sx q[1];
rz(-2.7374949) q[1];
rz(1.1683447) q[3];
sx q[3];
rz(-0.84512049) q[3];
sx q[3];
rz(2.890051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1027362) q[2];
sx q[2];
rz(-0.83013022) q[2];
sx q[2];
rz(-2.1170199) q[2];
rz(1.864795) q[3];
sx q[3];
rz(-1.5921311) q[3];
sx q[3];
rz(-1.5156486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(2.9643726) q[0];
sx q[0];
rz(-1.6765046) q[0];
sx q[0];
rz(0.033551034) q[0];
rz(-1.9112446) q[1];
sx q[1];
rz(-0.76554811) q[1];
sx q[1];
rz(-1.4039325) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0008154) q[0];
sx q[0];
rz(-2.7404832) q[0];
sx q[0];
rz(-0.51638575) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83547445) q[2];
sx q[2];
rz(-0.88047853) q[2];
sx q[2];
rz(-0.26963216) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0701323) q[1];
sx q[1];
rz(-1.4590108) q[1];
sx q[1];
rz(-1.8624572) q[1];
rz(-pi) q[2];
x q[2];
rz(0.085539354) q[3];
sx q[3];
rz(-0.88118689) q[3];
sx q[3];
rz(-2.6181108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6491062) q[2];
sx q[2];
rz(-2.2968569) q[2];
sx q[2];
rz(0.21326324) q[2];
rz(-0.02031859) q[3];
sx q[3];
rz(-1.820194) q[3];
sx q[3];
rz(0.4030574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7810818) q[0];
sx q[0];
rz(-1.1163982) q[0];
sx q[0];
rz(2.1257341) q[0];
rz(-1.5083183) q[1];
sx q[1];
rz(-2.1834178) q[1];
sx q[1];
rz(0.034084592) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7406697) q[0];
sx q[0];
rz(-1.4105083) q[0];
sx q[0];
rz(-2.7435061) q[0];
rz(-1.3047406) q[2];
sx q[2];
rz(-1.012946) q[2];
sx q[2];
rz(1.7709874) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9220229) q[1];
sx q[1];
rz(-1.8551738) q[1];
sx q[1];
rz(0.0071830458) q[1];
rz(-pi) q[2];
rz(0.94200763) q[3];
sx q[3];
rz(-1.994547) q[3];
sx q[3];
rz(2.8429192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.74026996) q[2];
sx q[2];
rz(-2.5677742) q[2];
sx q[2];
rz(-1.1711228) q[2];
rz(-1.3327538) q[3];
sx q[3];
rz(-1.7141902) q[3];
sx q[3];
rz(-0.12935054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68833441) q[0];
sx q[0];
rz(-1.9535221) q[0];
sx q[0];
rz(-0.43584287) q[0];
rz(-1.1054989) q[1];
sx q[1];
rz(-1.8445797) q[1];
sx q[1];
rz(-1.2058535) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5642731) q[0];
sx q[0];
rz(-2.3015129) q[0];
sx q[0];
rz(-2.8217836) q[0];
rz(-pi) q[1];
rz(-3.1255683) q[2];
sx q[2];
rz(-1.4756225) q[2];
sx q[2];
rz(-1.3118088) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0059433) q[1];
sx q[1];
rz(-1.597153) q[1];
sx q[1];
rz(0.51075682) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.033127012) q[3];
sx q[3];
rz(-1.1713542) q[3];
sx q[3];
rz(1.9314194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.897052) q[2];
sx q[2];
rz(-0.62770939) q[2];
sx q[2];
rz(-0.051606027) q[2];
rz(2.7339325) q[3];
sx q[3];
rz(-1.1838653) q[3];
sx q[3];
rz(-2.6962962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62149858) q[0];
sx q[0];
rz(-2.5243653) q[0];
sx q[0];
rz(2.4682585) q[0];
rz(0.78397059) q[1];
sx q[1];
rz(-1.4173123) q[1];
sx q[1];
rz(2.6046682) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1896742) q[0];
sx q[0];
rz(-3.0351743) q[0];
sx q[0];
rz(-0.90344067) q[0];
rz(1.0579254) q[2];
sx q[2];
rz(-0.29209902) q[2];
sx q[2];
rz(-0.7823173) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.19837025) q[1];
sx q[1];
rz(-2.3343625) q[1];
sx q[1];
rz(3.0109349) q[1];
x q[2];
rz(2.5506006) q[3];
sx q[3];
rz(-1.6984852) q[3];
sx q[3];
rz(-1.5203116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.15988222) q[2];
sx q[2];
rz(-1.8843702) q[2];
sx q[2];
rz(-2.9505777) q[2];
rz(-2.8602709) q[3];
sx q[3];
rz(-1.1912991) q[3];
sx q[3];
rz(2.7991926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0328338) q[0];
sx q[0];
rz(-2.7712951) q[0];
sx q[0];
rz(-0.38761815) q[0];
rz(-3.0265813) q[1];
sx q[1];
rz(-1.7128046) q[1];
sx q[1];
rz(-2.8040335) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.176929) q[0];
sx q[0];
rz(-1.3725855) q[0];
sx q[0];
rz(2.6508509) q[0];
rz(-0.93162025) q[2];
sx q[2];
rz(-1.6281623) q[2];
sx q[2];
rz(0.66040874) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2355289) q[1];
sx q[1];
rz(-1.4396832) q[1];
sx q[1];
rz(-0.99793418) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0739325) q[3];
sx q[3];
rz(-2.3849871) q[3];
sx q[3];
rz(-0.17186952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.18743029) q[2];
sx q[2];
rz(-0.53415161) q[2];
sx q[2];
rz(1.8743275) q[2];
rz(-2.3665442) q[3];
sx q[3];
rz(-1.6036443) q[3];
sx q[3];
rz(-2.3601941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9545814) q[0];
sx q[0];
rz(-1.0382074) q[0];
sx q[0];
rz(0.22098456) q[0];
rz(-2.2194608) q[1];
sx q[1];
rz(-1.8716967) q[1];
sx q[1];
rz(-2.1386713) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.636165) q[0];
sx q[0];
rz(-1.3190077) q[0];
sx q[0];
rz(-1.4279143) q[0];
x q[1];
rz(-2.2641364) q[2];
sx q[2];
rz(-2.4306731) q[2];
sx q[2];
rz(3.1117709) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3993966) q[1];
sx q[1];
rz(-2.4376215) q[1];
sx q[1];
rz(2.3854371) q[1];
x q[2];
rz(-0.27243154) q[3];
sx q[3];
rz(-2.2386132) q[3];
sx q[3];
rz(-1.9200793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6691436) q[2];
sx q[2];
rz(-2.3995235) q[2];
sx q[2];
rz(-2.9830902) q[2];
rz(2.6878099) q[3];
sx q[3];
rz(-0.78275371) q[3];
sx q[3];
rz(-2.2973072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52699387) q[0];
sx q[0];
rz(-2.232593) q[0];
sx q[0];
rz(-0.19113834) q[0];
rz(-0.29516164) q[1];
sx q[1];
rz(-0.8894397) q[1];
sx q[1];
rz(-0.89231649) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79217171) q[0];
sx q[0];
rz(-2.3755709) q[0];
sx q[0];
rz(-1.8269405) q[0];
rz(3.0478165) q[2];
sx q[2];
rz(-2.2095592) q[2];
sx q[2];
rz(-0.057657777) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.031242) q[1];
sx q[1];
rz(-0.24734766) q[1];
sx q[1];
rz(2.8691643) q[1];
rz(0.43182208) q[3];
sx q[3];
rz(-1.338827) q[3];
sx q[3];
rz(-0.98917978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3698547) q[2];
sx q[2];
rz(-0.83738804) q[2];
sx q[2];
rz(-2.2820293) q[2];
rz(-1.9178948) q[3];
sx q[3];
rz(-2.2151291) q[3];
sx q[3];
rz(-2.3971476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1214462) q[0];
sx q[0];
rz(-0.93562026) q[0];
sx q[0];
rz(1.390441) q[0];
rz(-1.4795115) q[1];
sx q[1];
rz(-2.6585487) q[1];
sx q[1];
rz(-1.9225635) q[1];
rz(-1.2791469) q[2];
sx q[2];
rz(-1.5445865) q[2];
sx q[2];
rz(-0.18538743) q[2];
rz(-1.5020694) q[3];
sx q[3];
rz(-1.6106265) q[3];
sx q[3];
rz(2.0911218) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
