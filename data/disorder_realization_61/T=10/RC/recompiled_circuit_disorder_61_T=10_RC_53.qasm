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
rz(2.5974098) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.777594) q[0];
sx q[0];
rz(-1.5033804) q[0];
sx q[0];
rz(1.534034) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0402096) q[2];
sx q[2];
rz(-1.5268491) q[2];
sx q[2];
rz(-1.6909388) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6737353) q[1];
sx q[1];
rz(-0.89092548) q[1];
sx q[1];
rz(-2.6152337) q[1];
rz(-pi) q[2];
rz(-0.028624264) q[3];
sx q[3];
rz(-1.5531566) q[3];
sx q[3];
rz(2.6580435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.071775285) q[2];
sx q[2];
rz(-1.2640307) q[2];
sx q[2];
rz(1.7791746) q[2];
rz(-3.1133364) q[3];
sx q[3];
rz(-1.7621721) q[3];
sx q[3];
rz(0.79022592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4966999) q[0];
sx q[0];
rz(-0.70514482) q[0];
sx q[0];
rz(1.171296) q[0];
rz(2.9303739) q[1];
sx q[1];
rz(-0.44208458) q[1];
sx q[1];
rz(-1.404095) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0140186) q[0];
sx q[0];
rz(-1.9919792) q[0];
sx q[0];
rz(-0.76571) q[0];
rz(1.7369466) q[2];
sx q[2];
rz(-2.0478021) q[2];
sx q[2];
rz(-1.1039066) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.055981936) q[1];
sx q[1];
rz(-1.2877052) q[1];
sx q[1];
rz(-1.594747) q[1];
x q[2];
rz(-2.6614463) q[3];
sx q[3];
rz(-1.4501791) q[3];
sx q[3];
rz(3.0062825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.30963787) q[2];
sx q[2];
rz(-1.4761304) q[2];
sx q[2];
rz(-0.94397604) q[2];
rz(2.5850463) q[3];
sx q[3];
rz(-2.6440933) q[3];
sx q[3];
rz(1.8054307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29207644) q[0];
sx q[0];
rz(-0.76382604) q[0];
sx q[0];
rz(1.4235494) q[0];
rz(-2.3220093) q[1];
sx q[1];
rz(-0.81033605) q[1];
sx q[1];
rz(0.56366411) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6113341) q[0];
sx q[0];
rz(-1.5000492) q[0];
sx q[0];
rz(0.57460873) q[0];
rz(-0.32391874) q[2];
sx q[2];
rz(-0.62506667) q[2];
sx q[2];
rz(2.6697086) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0020395) q[1];
sx q[1];
rz(-0.66322749) q[1];
sx q[1];
rz(1.0242277) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54791252) q[3];
sx q[3];
rz(-0.69763819) q[3];
sx q[3];
rz(-1.1988977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.50239572) q[2];
sx q[2];
rz(-2.0194619) q[2];
sx q[2];
rz(2.0813265) q[2];
rz(0.075573102) q[3];
sx q[3];
rz(-0.85791701) q[3];
sx q[3];
rz(0.11463541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.32325) q[0];
sx q[0];
rz(-2.8338354) q[0];
sx q[0];
rz(-2.9484205) q[0];
rz(1.5974143) q[1];
sx q[1];
rz(-1.9924106) q[1];
sx q[1];
rz(0.65778041) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4260027) q[0];
sx q[0];
rz(-1.7455187) q[0];
sx q[0];
rz(1.6013553) q[0];
x q[1];
rz(1.1409608) q[2];
sx q[2];
rz(-1.2865215) q[2];
sx q[2];
rz(0.62361275) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.12280497) q[1];
sx q[1];
rz(-2.0969166) q[1];
sx q[1];
rz(-0.29463525) q[1];
rz(1.2781906) q[3];
sx q[3];
rz(-1.3232104) q[3];
sx q[3];
rz(0.86865559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0956991) q[2];
sx q[2];
rz(-2.7313488) q[2];
sx q[2];
rz(2.0945385) q[2];
rz(1.0632769) q[3];
sx q[3];
rz(-1.1626817) q[3];
sx q[3];
rz(1.261238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7970153) q[0];
sx q[0];
rz(-0.319096) q[0];
sx q[0];
rz(-1.0239333) q[0];
rz(0.78760415) q[1];
sx q[1];
rz(-1.5757631) q[1];
sx q[1];
rz(-0.62686282) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89676566) q[0];
sx q[0];
rz(-1.5875495) q[0];
sx q[0];
rz(-1.6003952) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75682326) q[2];
sx q[2];
rz(-2.7177817) q[2];
sx q[2];
rz(2.8611081) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.353646) q[1];
sx q[1];
rz(-1.8783356) q[1];
sx q[1];
rz(-2.8278082) q[1];
x q[2];
rz(1.8623779) q[3];
sx q[3];
rz(-0.40816187) q[3];
sx q[3];
rz(0.51418958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.91288599) q[2];
sx q[2];
rz(-1.2330981) q[2];
sx q[2];
rz(-1.0181001) q[2];
rz(2.5300238) q[3];
sx q[3];
rz(-1.3221778) q[3];
sx q[3];
rz(-0.95156804) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6927032) q[0];
sx q[0];
rz(-1.3494116) q[0];
sx q[0];
rz(-1.942379) q[0];
rz(-1.9723643) q[1];
sx q[1];
rz(-1.0511845) q[1];
sx q[1];
rz(-0.32454023) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51844653) q[0];
sx q[0];
rz(-1.070797) q[0];
sx q[0];
rz(2.6536233) q[0];
x q[1];
rz(-0.65002273) q[2];
sx q[2];
rz(-2.2215002) q[2];
sx q[2];
rz(-1.0879773) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.74450103) q[1];
sx q[1];
rz(-2.607064) q[1];
sx q[1];
rz(1.6194653) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2732029) q[3];
sx q[3];
rz(-1.9801567) q[3];
sx q[3];
rz(-2.5149432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.67363182) q[2];
sx q[2];
rz(-1.8362074) q[2];
sx q[2];
rz(-1.2754053) q[2];
rz(-1.0547137) q[3];
sx q[3];
rz(-2.349699) q[3];
sx q[3];
rz(1.5462497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4642898) q[0];
sx q[0];
rz(-1.3339366) q[0];
sx q[0];
rz(-1.4087079) q[0];
rz(2.6858221) q[1];
sx q[1];
rz(-2.9401638) q[1];
sx q[1];
rz(1.9394402) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31583187) q[0];
sx q[0];
rz(-1.1734661) q[0];
sx q[0];
rz(2.2560918) q[0];
rz(0.75041109) q[2];
sx q[2];
rz(-1.8609036) q[2];
sx q[2];
rz(1.4587221) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.251293) q[1];
sx q[1];
rz(-1.6541462) q[1];
sx q[1];
rz(-0.092373089) q[1];
x q[2];
rz(-0.2644516) q[3];
sx q[3];
rz(-1.7489986) q[3];
sx q[3];
rz(1.4101392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1414286) q[2];
sx q[2];
rz(-1.2951853) q[2];
sx q[2];
rz(-0.84623519) q[2];
rz(0.49514654) q[3];
sx q[3];
rz(-0.64619243) q[3];
sx q[3];
rz(1.5406066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0705567) q[0];
sx q[0];
rz(-1.2905916) q[0];
sx q[0];
rz(2.0284247) q[0];
rz(-2.44599) q[1];
sx q[1];
rz(-2.7516987) q[1];
sx q[1];
rz(-1.5323458) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1217864) q[0];
sx q[0];
rz(-2.0397423) q[0];
sx q[0];
rz(-1.297834) q[0];
rz(-pi) q[1];
rz(-2.3959827) q[2];
sx q[2];
rz(-2.5157305) q[2];
sx q[2];
rz(1.2127884) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6032226) q[1];
sx q[1];
rz(-2.6909628) q[1];
sx q[1];
rz(-1.6166812) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4523776) q[3];
sx q[3];
rz(-2.1695671) q[3];
sx q[3];
rz(-2.1042202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1476851) q[2];
sx q[2];
rz(-2.9116178) q[2];
sx q[2];
rz(-2.2873986) q[2];
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
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5937186) q[0];
sx q[0];
rz(-2.6429208) q[0];
sx q[0];
rz(2.4771931) q[0];
rz(-1.9539072) q[1];
sx q[1];
rz(-0.70078754) q[1];
sx q[1];
rz(-1.857035) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5212095) q[0];
sx q[0];
rz(-1.8561583) q[0];
sx q[0];
rz(-0.83831212) q[0];
x q[1];
rz(1.9509435) q[2];
sx q[2];
rz(-1.8188307) q[2];
sx q[2];
rz(-1.5284577) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7582693) q[1];
sx q[1];
rz(-2.3136534) q[1];
sx q[1];
rz(-2.5187056) q[1];
rz(3.0751238) q[3];
sx q[3];
rz(-2.7507938) q[3];
sx q[3];
rz(-1.9399411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.61218843) q[2];
sx q[2];
rz(-2.2470784) q[2];
sx q[2];
rz(0.55076304) q[2];
rz(-2.4380056) q[3];
sx q[3];
rz(-2.1313322) q[3];
sx q[3];
rz(-1.6350869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
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
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6990307) q[0];
sx q[0];
rz(-0.5578707) q[0];
sx q[0];
rz(-0.2691707) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3264364) q[2];
sx q[2];
rz(-1.9197459) q[2];
sx q[2];
rz(-2.7100035) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.30291468) q[1];
sx q[1];
rz(-2.089114) q[1];
sx q[1];
rz(1.7811437) q[1];
rz(-0.46882792) q[3];
sx q[3];
rz(-2.3033597) q[3];
sx q[3];
rz(-1.9459141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49446517) q[2];
sx q[2];
rz(-0.72138849) q[2];
sx q[2];
rz(1.5987827) q[2];
rz(-2.4370082) q[3];
sx q[3];
rz(-1.6274118) q[3];
sx q[3];
rz(3.1295479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.96991878) q[2];
sx q[2];
rz(-2.6261332) q[2];
sx q[2];
rz(0.020608227) q[2];
rz(-2.963831) q[3];
sx q[3];
rz(-0.64571417) q[3];
sx q[3];
rz(-0.4971102) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
