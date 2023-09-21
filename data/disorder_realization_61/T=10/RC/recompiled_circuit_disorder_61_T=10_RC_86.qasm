OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4133889) q[0];
sx q[0];
rz(-1.1336741) q[0];
sx q[0];
rz(1.5925621) q[0];
rz(1.6917317) q[1];
sx q[1];
rz(5.6258968) q[1];
sx q[1];
rz(13.110553) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9323174) q[0];
sx q[0];
rz(-1.6074751) q[0];
sx q[0];
rz(-0.067461405) q[0];
rz(2.7317023) q[2];
sx q[2];
rz(-0.11046834) q[2];
sx q[2];
rz(-0.2875178) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4197598) q[1];
sx q[1];
rz(-2.3082323) q[1];
sx q[1];
rz(1.014773) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1129684) q[3];
sx q[3];
rz(-1.588436) q[3];
sx q[3];
rz(-2.6580435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.071775285) q[2];
sx q[2];
rz(-1.8775619) q[2];
sx q[2];
rz(1.7791746) q[2];
rz(3.1133364) q[3];
sx q[3];
rz(-1.3794206) q[3];
sx q[3];
rz(-2.3513667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.64489275) q[0];
sx q[0];
rz(-0.70514482) q[0];
sx q[0];
rz(-1.171296) q[0];
rz(2.9303739) q[1];
sx q[1];
rz(-2.6995081) q[1];
sx q[1];
rz(1.404095) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0140186) q[0];
sx q[0];
rz(-1.1496135) q[0];
sx q[0];
rz(0.76571) q[0];
x q[1];
rz(0.30971576) q[2];
sx q[2];
rz(-0.5030015) q[2];
sx q[2];
rz(-0.75370698) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6200871) q[1];
sx q[1];
rz(-1.5477991) q[1];
sx q[1];
rz(-2.8584245) q[1];
rz(2.6614463) q[3];
sx q[3];
rz(-1.4501791) q[3];
sx q[3];
rz(-3.0062825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8319548) q[2];
sx q[2];
rz(-1.4761304) q[2];
sx q[2];
rz(-0.94397604) q[2];
rz(-2.5850463) q[3];
sx q[3];
rz(-0.49749938) q[3];
sx q[3];
rz(-1.336162) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8495162) q[0];
sx q[0];
rz(-0.76382604) q[0];
sx q[0];
rz(-1.4235494) q[0];
rz(2.3220093) q[1];
sx q[1];
rz(-0.81033605) q[1];
sx q[1];
rz(2.5779285) q[1];
rz(pi/2) q[2];
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
rz(-1.3450422) q[2];
sx q[2];
rz(-0.98276897) q[2];
sx q[2];
rz(3.0622481) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0020395) q[1];
sx q[1];
rz(-0.66322749) q[1];
sx q[1];
rz(2.1173649) q[1];
x q[2];
rz(2.5936801) q[3];
sx q[3];
rz(-0.69763819) q[3];
sx q[3];
rz(-1.9426949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6391969) q[2];
sx q[2];
rz(-1.1221308) q[2];
sx q[2];
rz(-1.0602661) q[2];
rz(-0.075573102) q[3];
sx q[3];
rz(-0.85791701) q[3];
sx q[3];
rz(-0.11463541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.32325) q[0];
sx q[0];
rz(-2.8338354) q[0];
sx q[0];
rz(2.9484205) q[0];
rz(1.5974143) q[1];
sx q[1];
rz(-1.9924106) q[1];
sx q[1];
rz(0.65778041) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6000711) q[0];
sx q[0];
rz(-0.17734781) q[0];
sx q[0];
rz(2.970201) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8305956) q[2];
sx q[2];
rz(-1.9823091) q[2];
sx q[2];
rz(-2.3222773) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8448338) q[1];
sx q[1];
rz(-1.3169603) q[1];
sx q[1];
rz(-1.0253419) q[1];
x q[2];
rz(-2.8834881) q[3];
sx q[3];
rz(-1.2873642) q[3];
sx q[3];
rz(2.5131445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0956991) q[2];
sx q[2];
rz(-2.7313488) q[2];
sx q[2];
rz(2.0945385) q[2];
rz(-2.0783157) q[3];
sx q[3];
rz(-1.9789109) q[3];
sx q[3];
rz(-1.261238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7970153) q[0];
sx q[0];
rz(-2.8224967) q[0];
sx q[0];
rz(-2.1176594) q[0];
rz(0.78760415) q[1];
sx q[1];
rz(-1.5658295) q[1];
sx q[1];
rz(-2.5147298) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67353467) q[0];
sx q[0];
rz(-1.600391) q[0];
sx q[0];
rz(-3.1248321) q[0];
rz(-2.3847694) q[2];
sx q[2];
rz(-2.7177817) q[2];
sx q[2];
rz(-0.28048453) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7879466) q[1];
sx q[1];
rz(-1.8783356) q[1];
sx q[1];
rz(-2.8278082) q[1];
x q[2];
rz(-3.0179126) q[3];
sx q[3];
rz(-1.9607753) q[3];
sx q[3];
rz(0.19815138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2287067) q[2];
sx q[2];
rz(-1.9084946) q[2];
sx q[2];
rz(1.0181001) q[2];
rz(2.5300238) q[3];
sx q[3];
rz(-1.3221778) q[3];
sx q[3];
rz(2.1900246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6927032) q[0];
sx q[0];
rz(-1.3494116) q[0];
sx q[0];
rz(-1.1992136) q[0];
rz(-1.9723643) q[1];
sx q[1];
rz(-2.0904082) q[1];
sx q[1];
rz(0.32454023) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31774662) q[0];
sx q[0];
rz(-2.4577603) q[0];
sx q[0];
rz(-2.2800287) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2424477) q[2];
sx q[2];
rz(-0.88485826) q[2];
sx q[2];
rz(-1.9859973) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4536344) q[1];
sx q[1];
rz(-2.1046241) q[1];
sx q[1];
rz(0.028793528) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5472774) q[3];
sx q[3];
rz(-0.50110498) q[3];
sx q[3];
rz(-1.8584115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6773029) q[0];
sx q[0];
rz(-1.807656) q[0];
sx q[0];
rz(-1.4087079) q[0];
rz(0.45577058) q[1];
sx q[1];
rz(-2.9401638) q[1];
sx q[1];
rz(-1.9394402) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31583187) q[0];
sx q[0];
rz(-1.1734661) q[0];
sx q[0];
rz(0.8855008) q[0];
x q[1];
rz(-2.7289594) q[2];
sx q[2];
rz(-0.79421439) q[2];
sx q[2];
rz(0.40976322) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.251293) q[1];
sx q[1];
rz(-1.4874465) q[1];
sx q[1];
rz(3.0492196) q[1];
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
rz(-pi/2) q[1];
rz(1.000164) q[2];
sx q[2];
rz(-1.2951853) q[2];
sx q[2];
rz(0.84623519) q[2];
rz(2.6464461) q[3];
sx q[3];
rz(-2.4954002) q[3];
sx q[3];
rz(1.5406066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-3.0705567) q[0];
sx q[0];
rz(-1.2905916) q[0];
sx q[0];
rz(-2.0284247) q[0];
rz(0.69560266) q[1];
sx q[1];
rz(-0.38989392) q[1];
sx q[1];
rz(-1.6092469) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42513645) q[0];
sx q[0];
rz(-1.3279337) q[0];
sx q[0];
rz(0.48432414) q[0];
rz(0.74560994) q[2];
sx q[2];
rz(-2.5157305) q[2];
sx q[2];
rz(1.2127884) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6541877) q[1];
sx q[1];
rz(-1.1206756) q[1];
sx q[1];
rz(-3.1194035) q[1];
x q[2];
rz(1.6892151) q[3];
sx q[3];
rz(-0.97202557) q[3];
sx q[3];
rz(-2.1042202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9939076) q[2];
sx q[2];
rz(-0.2299749) q[2];
sx q[2];
rz(2.2873986) q[2];
rz(1.9130075) q[3];
sx q[3];
rz(-1.5649786) q[3];
sx q[3];
rz(1.4860229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5478741) q[0];
sx q[0];
rz(-0.49867189) q[0];
sx q[0];
rz(-0.6643995) q[0];
rz(-1.9539072) q[1];
sx q[1];
rz(-0.70078754) q[1];
sx q[1];
rz(-1.857035) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4391543) q[0];
sx q[0];
rz(-0.87411532) q[0];
sx q[0];
rz(-2.7657763) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9719073) q[2];
sx q[2];
rz(-0.45058695) q[2];
sx q[2];
rz(2.633) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.43209546) q[1];
sx q[1];
rz(-0.92952432) q[1];
sx q[1];
rz(2.1367367) q[1];
rz(1.5434389) q[3];
sx q[3];
rz(-1.1809071) q[3];
sx q[3];
rz(-2.011812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5294042) q[2];
sx q[2];
rz(-0.89451423) q[2];
sx q[2];
rz(2.5908296) q[2];
rz(0.70358706) q[3];
sx q[3];
rz(-1.0102605) q[3];
sx q[3];
rz(1.6350869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4187014) q[0];
sx q[0];
rz(-2.7846865) q[0];
sx q[0];
rz(-0.044145949) q[0];
rz(1.528953) q[1];
sx q[1];
rz(-1.2395369) q[1];
sx q[1];
rz(-2.3619161) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6990307) q[0];
sx q[0];
rz(-2.583722) q[0];
sx q[0];
rz(-2.872422) q[0];
rz(-0.463562) q[2];
sx q[2];
rz(-2.2710685) q[2];
sx q[2];
rz(2.314032) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.30291468) q[1];
sx q[1];
rz(-2.089114) q[1];
sx q[1];
rz(1.7811437) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78124222) q[3];
sx q[3];
rz(-1.2282073) q[3];
sx q[3];
rz(-0.70171802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.49446517) q[2];
sx q[2];
rz(-2.4202042) q[2];
sx q[2];
rz(1.5987827) q[2];
rz(0.70458448) q[3];
sx q[3];
rz(-1.5141809) q[3];
sx q[3];
rz(0.012044756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
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
rz(1.7977057) q[1];
sx q[1];
rz(-1.4812891) q[1];
sx q[1];
rz(-0.67566009) q[1];
rz(0.30998183) q[2];
sx q[2];
rz(-1.1520755) q[2];
sx q[2];
rz(0.68785695) q[2];
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