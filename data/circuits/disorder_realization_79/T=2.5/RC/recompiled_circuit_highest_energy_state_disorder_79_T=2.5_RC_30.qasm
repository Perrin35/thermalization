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
rz(2.6464638) q[0];
sx q[0];
rz(9.169307) q[0];
rz(2.7893692) q[1];
sx q[1];
rz(-1.398634) q[1];
sx q[1];
rz(-1.7359098) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2737925) q[0];
sx q[0];
rz(-1.5711014) q[0];
sx q[0];
rz(0.0034694151) q[0];
rz(-pi) q[1];
x q[1];
rz(2.954835) q[2];
sx q[2];
rz(-0.52350145) q[2];
sx q[2];
rz(1.40846) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5454332) q[1];
sx q[1];
rz(-1.1876629) q[1];
sx q[1];
rz(3.1379164) q[1];
rz(-pi) q[2];
rz(2.5735503) q[3];
sx q[3];
rz(-1.5431817) q[3];
sx q[3];
rz(1.9909007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.63456717) q[2];
sx q[2];
rz(-3.1105803) q[2];
sx q[2];
rz(2.3090889) q[2];
rz(-0.80820525) q[3];
sx q[3];
rz(-0.015232239) q[3];
sx q[3];
rz(-2.8830849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.38043624) q[0];
sx q[0];
rz(-2.7988837) q[0];
sx q[0];
rz(2.9535182) q[0];
rz(3.0694118) q[1];
sx q[1];
rz(-1.0344104) q[1];
sx q[1];
rz(1.5123051) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3609027) q[0];
sx q[0];
rz(-1.6854428) q[0];
sx q[0];
rz(0.33495442) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1109643) q[2];
sx q[2];
rz(-1.616744) q[2];
sx q[2];
rz(-0.46389889) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.20023274) q[1];
sx q[1];
rz(-1.5116475) q[1];
sx q[1];
rz(0.099353921) q[1];
rz(2.1339122) q[3];
sx q[3];
rz(-2.2950239) q[3];
sx q[3];
rz(-0.40562848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9709116) q[2];
sx q[2];
rz(-0.6450246) q[2];
sx q[2];
rz(-1.8485273) q[2];
rz(2.7944148) q[3];
sx q[3];
rz(-0.2694338) q[3];
sx q[3];
rz(0.080304317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8186571) q[0];
sx q[0];
rz(-1.8523536) q[0];
sx q[0];
rz(-2.6710508) q[0];
rz(-1.200354) q[1];
sx q[1];
rz(-2.410694) q[1];
sx q[1];
rz(1.2568731) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82913933) q[0];
sx q[0];
rz(-1.0395673) q[0];
sx q[0];
rz(-0.4108528) q[0];
x q[1];
rz(0.09229163) q[2];
sx q[2];
rz(-1.4362291) q[2];
sx q[2];
rz(1.2876266) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5804435) q[1];
sx q[1];
rz(-3.1270091) q[1];
sx q[1];
rz(-1.4766978) q[1];
rz(2.8506365) q[3];
sx q[3];
rz(-2.1621063) q[3];
sx q[3];
rz(1.3475498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5990126) q[2];
sx q[2];
rz(-2.161669) q[2];
sx q[2];
rz(-0.92213321) q[2];
rz(-1.7982091) q[3];
sx q[3];
rz(-2.193439) q[3];
sx q[3];
rz(1.62742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.2682997) q[0];
sx q[0];
rz(-2.101185) q[0];
sx q[0];
rz(2.701395) q[0];
rz(-1.6104376) q[1];
sx q[1];
rz(-1.4819769) q[1];
sx q[1];
rz(2.8940315) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9151814) q[0];
sx q[0];
rz(-1.1410211) q[0];
sx q[0];
rz(-2.129401) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4804968) q[2];
sx q[2];
rz(-1.7396428) q[2];
sx q[2];
rz(0.24659469) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9369071) q[1];
sx q[1];
rz(-1.2720894) q[1];
sx q[1];
rz(-3.0644817) q[1];
rz(1.5400793) q[3];
sx q[3];
rz(-0.92245543) q[3];
sx q[3];
rz(1.101644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2441248) q[2];
sx q[2];
rz(-1.030587) q[2];
sx q[2];
rz(-1.6991276) q[2];
rz(2.9407732) q[3];
sx q[3];
rz(-1.9229869) q[3];
sx q[3];
rz(-1.1403181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9852801) q[0];
sx q[0];
rz(-0.15287481) q[0];
sx q[0];
rz(0.54541624) q[0];
rz(-0.044005752) q[1];
sx q[1];
rz(-3.1237055) q[1];
sx q[1];
rz(0.47637475) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8808419) q[0];
sx q[0];
rz(-1.6549661) q[0];
sx q[0];
rz(0.24411402) q[0];
rz(-pi) q[1];
rz(-2.3628391) q[2];
sx q[2];
rz(-2.5980332) q[2];
sx q[2];
rz(-1.5875848) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.68639437) q[1];
sx q[1];
rz(-2.2664811) q[1];
sx q[1];
rz(-0.2376539) q[1];
x q[2];
rz(2.0799177) q[3];
sx q[3];
rz(-1.4533852) q[3];
sx q[3];
rz(-1.5985296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5215317) q[2];
sx q[2];
rz(-2.4510577) q[2];
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
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5905269) q[0];
sx q[0];
rz(-1.5883625) q[0];
sx q[0];
rz(1.5565514) q[0];
rz(0.12697728) q[1];
sx q[1];
rz(-1.8366837) q[1];
sx q[1];
rz(0.089687673) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68192775) q[0];
sx q[0];
rz(-1.0583504) q[0];
sx q[0];
rz(-0.95109032) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96025567) q[2];
sx q[2];
rz(-0.88443236) q[2];
sx q[2];
rz(-2.6853564) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6818004) q[1];
sx q[1];
rz(-2.5885317) q[1];
sx q[1];
rz(-0.44981594) q[1];
x q[2];
rz(0.50181194) q[3];
sx q[3];
rz(-2.4779421) q[3];
sx q[3];
rz(2.7743031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.108532) q[2];
sx q[2];
rz(-2.8046799) q[2];
sx q[2];
rz(2.8899371) q[2];
rz(-0.039994914) q[3];
sx q[3];
rz(-2.7997041) q[3];
sx q[3];
rz(2.6778636) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43117487) q[0];
sx q[0];
rz(-0.091910563) q[0];
sx q[0];
rz(-0.41496667) q[0];
rz(1.4893432) q[1];
sx q[1];
rz(-3.0840315) q[1];
sx q[1];
rz(-2.8156978) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.519442) q[0];
sx q[0];
rz(-1.2076734) q[0];
sx q[0];
rz(-1.0958998) q[0];
rz(-pi) q[1];
rz(-1.7275179) q[2];
sx q[2];
rz(-0.5236519) q[2];
sx q[2];
rz(2.120842) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6025117) q[1];
sx q[1];
rz(-1.314114) q[1];
sx q[1];
rz(-0.51694444) q[1];
rz(-pi) q[2];
rz(1.8741499) q[3];
sx q[3];
rz(-1.6432495) q[3];
sx q[3];
rz(0.8006351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9140279) q[2];
sx q[2];
rz(-1.3332557) q[2];
sx q[2];
rz(1.0464767) q[2];
rz(-2.922831) q[3];
sx q[3];
rz(-2.1121787) q[3];
sx q[3];
rz(-1.7441162) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47531146) q[0];
sx q[0];
rz(-3.1240211) q[0];
sx q[0];
rz(-0.46324357) q[0];
rz(0.26495588) q[1];
sx q[1];
rz(-0.0015365096) q[1];
sx q[1];
rz(-1.6418246) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87141052) q[0];
sx q[0];
rz(-2.0034916) q[0];
sx q[0];
rz(3.1209594) q[0];
rz(-pi) q[1];
rz(-2.6959571) q[2];
sx q[2];
rz(-1.6204973) q[2];
sx q[2];
rz(-1.4193064) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8913938) q[1];
sx q[1];
rz(-1.7925182) q[1];
sx q[1];
rz(3.1117597) q[1];
x q[2];
rz(-1.8231976) q[3];
sx q[3];
rz(-1.5034672) q[3];
sx q[3];
rz(-2.61781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2396592) q[2];
sx q[2];
rz(-0.45238164) q[2];
sx q[2];
rz(1.3182013) q[2];
rz(3.1384595) q[3];
sx q[3];
rz(-1.2001218) q[3];
sx q[3];
rz(-1.1656632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5779293) q[0];
sx q[0];
rz(-2.2888373) q[0];
sx q[0];
rz(0.71459115) q[0];
rz(-0.21358061) q[1];
sx q[1];
rz(-0.029284632) q[1];
sx q[1];
rz(-1.9307131) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1909669) q[0];
sx q[0];
rz(-1.4651148) q[0];
sx q[0];
rz(-1.2638541) q[0];
rz(-1.9186594) q[2];
sx q[2];
rz(-2.2360254) q[2];
sx q[2];
rz(-0.72345483) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.2106194) q[1];
sx q[1];
rz(-0.60928357) q[1];
sx q[1];
rz(2.0740725) q[1];
rz(2.2641422) q[3];
sx q[3];
rz(-2.8971379) q[3];
sx q[3];
rz(-0.86737421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.96656281) q[2];
sx q[2];
rz(-0.021447072) q[2];
sx q[2];
rz(2.9435797) q[2];
rz(-2.7909732) q[3];
sx q[3];
rz(-1.5712761) q[3];
sx q[3];
rz(-1.9982136) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38458934) q[0];
sx q[0];
rz(-2.2191255) q[0];
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
rz(0.40648288) q[0];
sx q[0];
rz(-2.2648659) q[0];
sx q[0];
rz(-2.8227692) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5727881) q[2];
sx q[2];
rz(-1.4395166) q[2];
sx q[2];
rz(-0.37007331) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1471794) q[1];
sx q[1];
rz(-0.1754027) q[1];
sx q[1];
rz(2.0270623) q[1];
rz(1.1636825) q[3];
sx q[3];
rz(-1.3729248) q[3];
sx q[3];
rz(1.0033293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1018579) q[2];
sx q[2];
rz(-0.27145824) q[2];
sx q[2];
rz(-2.6177935) q[2];
rz(1.0456746) q[3];
sx q[3];
rz(-1.8664675) q[3];
sx q[3];
rz(-0.4642134) q[3];
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
rz(-2.4986497) q[2];
sx q[2];
rz(-1.8131154) q[2];
sx q[2];
rz(-2.9396069) q[2];
rz(1.3741334) q[3];
sx q[3];
rz(-2.0296367) q[3];
sx q[3];
rz(2.2214132) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
