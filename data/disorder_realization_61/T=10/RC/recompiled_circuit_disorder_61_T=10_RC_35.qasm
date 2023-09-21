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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13574164) q[0];
sx q[0];
rz(-0.07677456) q[0];
sx q[0];
rz(-2.6430921) q[0];
x q[1];
rz(-1.6149701) q[2];
sx q[2];
rz(-1.4695115) q[2];
sx q[2];
rz(3.0169808) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.45317263) q[1];
sx q[1];
rz(-1.9721713) q[1];
sx q[1];
rz(-0.81897264) q[1];
rz(2.5892026) q[3];
sx q[3];
rz(-0.033621764) q[3];
sx q[3];
rz(0.53510964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.071775285) q[2];
sx q[2];
rz(-1.8775619) q[2];
sx q[2];
rz(1.3624181) q[2];
rz(-3.1133364) q[3];
sx q[3];
rz(-1.3794206) q[3];
sx q[3];
rz(2.3513667) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4966999) q[0];
sx q[0];
rz(-2.4364478) q[0];
sx q[0];
rz(1.171296) q[0];
rz(-2.9303739) q[1];
sx q[1];
rz(-0.44208458) q[1];
sx q[1];
rz(1.404095) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1275741) q[0];
sx q[0];
rz(-1.1496135) q[0];
sx q[0];
rz(-0.76571) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48268433) q[2];
sx q[2];
rz(-1.4233372) q[2];
sx q[2];
rz(-0.54373103) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.055981936) q[1];
sx q[1];
rz(-1.2877052) q[1];
sx q[1];
rz(-1.5468456) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7066129) q[3];
sx q[3];
rz(-2.0471626) q[3];
sx q[3];
rz(1.3729031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8319548) q[2];
sx q[2];
rz(-1.4761304) q[2];
sx q[2];
rz(-2.1976166) q[2];
rz(2.5850463) q[3];
sx q[3];
rz(-0.49749938) q[3];
sx q[3];
rz(-1.8054307) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8495162) q[0];
sx q[0];
rz(-0.76382604) q[0];
sx q[0];
rz(1.7180432) q[0];
rz(-0.81958333) q[1];
sx q[1];
rz(-2.3312566) q[1];
sx q[1];
rz(0.56366411) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99479988) q[0];
sx q[0];
rz(-2.1437862) q[0];
sx q[0];
rz(-1.6550199) q[0];
rz(-1.7965505) q[2];
sx q[2];
rz(-0.98276897) q[2];
sx q[2];
rz(0.079344582) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7970265) q[1];
sx q[1];
rz(-1.016942) q[1];
sx q[1];
rz(2.7558541) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5205069) q[3];
sx q[3];
rz(-1.2295782) q[3];
sx q[3];
rz(-3.0760471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.50239572) q[2];
sx q[2];
rz(-2.0194619) q[2];
sx q[2];
rz(-1.0602661) q[2];
rz(-0.075573102) q[3];
sx q[3];
rz(-0.85791701) q[3];
sx q[3];
rz(3.0269572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.32325) q[0];
sx q[0];
rz(-2.8338354) q[0];
sx q[0];
rz(0.19317214) q[0];
rz(-1.5441783) q[1];
sx q[1];
rz(-1.9924106) q[1];
sx q[1];
rz(-2.4838122) q[1];
rz(pi/2) q[2];
sx q[2];
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
rz(-1.64536) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.12280497) q[1];
sx q[1];
rz(-2.0969166) q[1];
sx q[1];
rz(2.8469574) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8834881) q[3];
sx q[3];
rz(-1.8542284) q[3];
sx q[3];
rz(-2.5131445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0956991) q[2];
sx q[2];
rz(-0.4102439) q[2];
sx q[2];
rz(-1.0470541) q[2];
rz(-1.0632769) q[3];
sx q[3];
rz(-1.9789109) q[3];
sx q[3];
rz(1.261238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3445774) q[0];
sx q[0];
rz(-0.319096) q[0];
sx q[0];
rz(1.0239333) q[0];
rz(-2.3539885) q[1];
sx q[1];
rz(-1.5658295) q[1];
sx q[1];
rz(-2.5147298) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1889362) q[0];
sx q[0];
rz(-0.034010012) q[0];
sx q[0];
rz(2.0859499) q[0];
rz(-pi) q[1];
rz(2.3847694) q[2];
sx q[2];
rz(-0.42381091) q[2];
sx q[2];
rz(2.8611081) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.353646) q[1];
sx q[1];
rz(-1.8783356) q[1];
sx q[1];
rz(0.31378444) q[1];
rz(-1.8623779) q[3];
sx q[3];
rz(-2.7334308) q[3];
sx q[3];
rz(0.51418958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.91288599) q[2];
sx q[2];
rz(-1.2330981) q[2];
sx q[2];
rz(-2.1234925) q[2];
rz(-0.61156887) q[3];
sx q[3];
rz(-1.8194149) q[3];
sx q[3];
rz(-2.1900246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488895) q[0];
sx q[0];
rz(-1.792181) q[0];
sx q[0];
rz(-1.1992136) q[0];
rz(-1.9723643) q[1];
sx q[1];
rz(-1.0511845) q[1];
sx q[1];
rz(-0.32454023) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3015315) q[0];
sx q[0];
rz(-1.1467629) q[0];
sx q[0];
rz(1.016894) q[0];
rz(-pi) q[1];
rz(-0.89914497) q[2];
sx q[2];
rz(-2.2567344) q[2];
sx q[2];
rz(-1.9859973) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3970916) q[1];
sx q[1];
rz(-2.607064) q[1];
sx q[1];
rz(1.5221273) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42602496) q[3];
sx q[3];
rz(-1.8431438) q[3];
sx q[3];
rz(-0.8226738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67363182) q[2];
sx q[2];
rz(-1.3053852) q[2];
sx q[2];
rz(1.2754053) q[2];
rz(-1.0547137) q[3];
sx q[3];
rz(-0.79189363) q[3];
sx q[3];
rz(1.595343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4642898) q[0];
sx q[0];
rz(-1.3339366) q[0];
sx q[0];
rz(-1.4087079) q[0];
rz(-2.6858221) q[1];
sx q[1];
rz(-2.9401638) q[1];
sx q[1];
rz(-1.9394402) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81284886) q[0];
sx q[0];
rz(-0.7757196) q[0];
sx q[0];
rz(-2.1562955) q[0];
rz(-1.9583148) q[2];
sx q[2];
rz(-2.2829208) q[2];
sx q[2];
rz(0.14856635) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.68820885) q[1];
sx q[1];
rz(-1.4787448) q[1];
sx q[1];
rz(1.6545014) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7552745) q[3];
sx q[3];
rz(-1.8309621) q[3];
sx q[3];
rz(0.20862143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.000164) q[2];
sx q[2];
rz(-1.8464073) q[2];
sx q[2];
rz(-0.84623519) q[2];
rz(-0.49514654) q[3];
sx q[3];
rz(-2.4954002) q[3];
sx q[3];
rz(1.5406066) q[3];
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
rz(-3.0705567) q[0];
sx q[0];
rz(-1.2905916) q[0];
sx q[0];
rz(-2.0284247) q[0];
rz(2.44599) q[1];
sx q[1];
rz(-2.7516987) q[1];
sx q[1];
rz(-1.6092469) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1217864) q[0];
sx q[0];
rz(-1.1018503) q[0];
sx q[0];
rz(1.297834) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3959827) q[2];
sx q[2];
rz(-2.5157305) q[2];
sx q[2];
rz(-1.9288043) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.073736139) q[1];
sx q[1];
rz(-1.590775) q[1];
sx q[1];
rz(-2.0210135) q[1];
rz(1.6892151) q[3];
sx q[3];
rz(-0.97202557) q[3];
sx q[3];
rz(-2.1042202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1476851) q[2];
sx q[2];
rz(-2.9116178) q[2];
sx q[2];
rz(-0.85419401) q[2];
rz(-1.2285852) q[3];
sx q[3];
rz(-1.5649786) q[3];
sx q[3];
rz(1.4860229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5478741) q[0];
sx q[0];
rz(-0.49867189) q[0];
sx q[0];
rz(2.4771931) q[0];
rz(1.1876855) q[1];
sx q[1];
rz(-0.70078754) q[1];
sx q[1];
rz(1.2845576) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70243835) q[0];
sx q[0];
rz(-0.87411532) q[0];
sx q[0];
rz(0.37581635) q[0];
rz(0.9719073) q[2];
sx q[2];
rz(-0.45058695) q[2];
sx q[2];
rz(-0.50859261) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.43209546) q[1];
sx q[1];
rz(-0.92952432) q[1];
sx q[1];
rz(1.0048559) q[1];
rz(-3.0751238) q[3];
sx q[3];
rz(-0.39079881) q[3];
sx q[3];
rz(-1.9399411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5294042) q[2];
sx q[2];
rz(-0.89451423) q[2];
sx q[2];
rz(-0.55076304) q[2];
rz(-2.4380056) q[3];
sx q[3];
rz(-1.0102605) q[3];
sx q[3];
rz(-1.5065058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7228912) q[0];
sx q[0];
rz(-2.7846865) q[0];
sx q[0];
rz(-0.044145949) q[0];
rz(1.6126397) q[1];
sx q[1];
rz(-1.9020558) q[1];
sx q[1];
rz(-2.3619161) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89833242) q[0];
sx q[0];
rz(-1.7120449) q[0];
sx q[0];
rz(0.54153533) q[0];
x q[1];
rz(-2.3264364) q[2];
sx q[2];
rz(-1.9197459) q[2];
sx q[2];
rz(-0.43158917) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.70982198) q[1];
sx q[1];
rz(-2.5858472) q[1];
sx q[1];
rz(-0.35094786) q[1];
rz(-pi) q[2];
rz(0.78124222) q[3];
sx q[3];
rz(-1.2282073) q[3];
sx q[3];
rz(2.4398746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.49446517) q[2];
sx q[2];
rz(-2.4202042) q[2];
sx q[2];
rz(-1.54281) q[2];
rz(-2.4370082) q[3];
sx q[3];
rz(-1.5141809) q[3];
sx q[3];
rz(-3.1295479) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71173944) q[0];
sx q[0];
rz(-1.585351) q[0];
sx q[0];
rz(2.4899695) q[0];
rz(1.7977057) q[1];
sx q[1];
rz(-1.4812891) q[1];
sx q[1];
rz(-0.67566009) q[1];
rz(1.1336397) q[2];
sx q[2];
rz(-1.2883678) q[2];
sx q[2];
rz(-1.012445) q[2];
rz(-1.4383437) q[3];
sx q[3];
rz(-2.2046897) q[3];
sx q[3];
rz(-0.71837367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
