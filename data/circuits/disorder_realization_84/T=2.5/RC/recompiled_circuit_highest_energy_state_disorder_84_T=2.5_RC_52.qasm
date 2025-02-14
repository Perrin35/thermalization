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
rz(2.5166729) q[0];
sx q[0];
rz(-1.807037) q[0];
sx q[0];
rz(0.053243756) q[0];
rz(-2.6553395) q[1];
sx q[1];
rz(-3.0142205) q[1];
sx q[1];
rz(-1.706634) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4006179) q[0];
sx q[0];
rz(-1.8549998) q[0];
sx q[0];
rz(1.7114729) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4515013) q[2];
sx q[2];
rz(-1.3053467) q[2];
sx q[2];
rz(2.1448898) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8088805) q[1];
sx q[1];
rz(-1.1774027) q[1];
sx q[1];
rz(-1.62864) q[1];
rz(-pi) q[2];
rz(0.49519914) q[3];
sx q[3];
rz(-0.99118587) q[3];
sx q[3];
rz(-0.51900348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3367553) q[2];
sx q[2];
rz(-1.7982322) q[2];
sx q[2];
rz(-0.27377823) q[2];
rz(1.9233507) q[3];
sx q[3];
rz(-0.67353407) q[3];
sx q[3];
rz(0.28606733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1186721) q[0];
sx q[0];
rz(-2.3781222) q[0];
sx q[0];
rz(0.77962312) q[0];
rz(-0.66462213) q[1];
sx q[1];
rz(-1.9326991) q[1];
sx q[1];
rz(-1.2453311) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0459291) q[0];
sx q[0];
rz(-1.6324658) q[0];
sx q[0];
rz(1.685919) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0491291) q[2];
sx q[2];
rz(-0.86148724) q[2];
sx q[2];
rz(-0.48395983) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7088931) q[1];
sx q[1];
rz(-1.782976) q[1];
sx q[1];
rz(3.0794705) q[1];
rz(-pi) q[2];
rz(-0.20540463) q[3];
sx q[3];
rz(-1.0213791) q[3];
sx q[3];
rz(-0.37976625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7094884) q[2];
sx q[2];
rz(-0.69807845) q[2];
sx q[2];
rz(-3.091605) q[2];
rz(1.6677808) q[3];
sx q[3];
rz(-1.7260909) q[3];
sx q[3];
rz(-2.2059691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98899984) q[0];
sx q[0];
rz(-2.279156) q[0];
sx q[0];
rz(1.9631901) q[0];
rz(1.9940469) q[1];
sx q[1];
rz(-0.18745628) q[1];
sx q[1];
rz(2.5386834) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1499304) q[0];
sx q[0];
rz(-1.7310744) q[0];
sx q[0];
rz(-0.10326881) q[0];
rz(-2.1772867) q[2];
sx q[2];
rz(-2.3539001) q[2];
sx q[2];
rz(0.76157969) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.29290043) q[1];
sx q[1];
rz(-0.91182263) q[1];
sx q[1];
rz(-2.1007358) q[1];
rz(-pi) q[2];
rz(-2.6125978) q[3];
sx q[3];
rz(-2.1784665) q[3];
sx q[3];
rz(-1.120847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9511562) q[2];
sx q[2];
rz(-0.78593212) q[2];
sx q[2];
rz(-0.3832761) q[2];
rz(-1.8303309) q[3];
sx q[3];
rz(-2.0537328) q[3];
sx q[3];
rz(-3.0200628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.424778) q[0];
sx q[0];
rz(-0.99262339) q[0];
sx q[0];
rz(0.92856652) q[0];
rz(2.3021452) q[1];
sx q[1];
rz(-1.3176368) q[1];
sx q[1];
rz(-2.3323257) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8438551) q[0];
sx q[0];
rz(-1.017573) q[0];
sx q[0];
rz(-2.7064302) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4014612) q[2];
sx q[2];
rz(-1.1462829) q[2];
sx q[2];
rz(-2.9237539) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.88621186) q[1];
sx q[1];
rz(-1.6396697) q[1];
sx q[1];
rz(1.3248454) q[1];
rz(-pi) q[2];
rz(2.2146642) q[3];
sx q[3];
rz(-1.2371105) q[3];
sx q[3];
rz(2.9179171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3889918) q[2];
sx q[2];
rz(-1.5316803) q[2];
sx q[2];
rz(1.4468225) q[2];
rz(-1.1270479) q[3];
sx q[3];
rz(-2.4067252) q[3];
sx q[3];
rz(-2.5126422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22555722) q[0];
sx q[0];
rz(-1.1981523) q[0];
sx q[0];
rz(-1.7115364) q[0];
rz(-2.9043708) q[1];
sx q[1];
rz(-2.1784541) q[1];
sx q[1];
rz(2.846834) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4859516) q[0];
sx q[0];
rz(-1.427934) q[0];
sx q[0];
rz(-0.44012196) q[0];
rz(0.37173231) q[2];
sx q[2];
rz(-2.1991538) q[2];
sx q[2];
rz(-2.6184751) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7511616) q[1];
sx q[1];
rz(-2.0195578) q[1];
sx q[1];
rz(2.6780891) q[1];
rz(-pi) q[2];
rz(0.2278403) q[3];
sx q[3];
rz(-0.58594847) q[3];
sx q[3];
rz(-0.70825746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1416867) q[2];
sx q[2];
rz(-0.20076951) q[2];
sx q[2];
rz(-3.0086009) q[2];
rz(2.3717132) q[3];
sx q[3];
rz(-1.2917638) q[3];
sx q[3];
rz(0.31908116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28984508) q[0];
sx q[0];
rz(-0.3388437) q[0];
sx q[0];
rz(1.4122562) q[0];
rz(3.0899561) q[1];
sx q[1];
rz(-2.5187571) q[1];
sx q[1];
rz(0.19270611) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7190921) q[0];
sx q[0];
rz(-2.0845045) q[0];
sx q[0];
rz(1.5925608) q[0];
rz(-pi) q[1];
rz(0.73130812) q[2];
sx q[2];
rz(-2.0530862) q[2];
sx q[2];
rz(2.363229) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.440408) q[1];
sx q[1];
rz(-1.1478979) q[1];
sx q[1];
rz(-0.50114034) q[1];
x q[2];
rz(-2.5000743) q[3];
sx q[3];
rz(-1.8810442) q[3];
sx q[3];
rz(-1.9812685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0685588) q[2];
sx q[2];
rz(-1.2898338) q[2];
sx q[2];
rz(1.7603091) q[2];
rz(-2.6265788) q[3];
sx q[3];
rz(-2.2541855) q[3];
sx q[3];
rz(-1.7611586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16159049) q[0];
sx q[0];
rz(-1.9965633) q[0];
sx q[0];
rz(-0.34647754) q[0];
rz(1.0193635) q[1];
sx q[1];
rz(-0.66277021) q[1];
sx q[1];
rz(1.4541218) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2812735) q[0];
sx q[0];
rz(-0.77206445) q[0];
sx q[0];
rz(-0.66769974) q[0];
rz(-2.2367661) q[2];
sx q[2];
rz(-0.77205333) q[2];
sx q[2];
rz(1.1386629) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7123966) q[1];
sx q[1];
rz(-2.0919261) q[1];
sx q[1];
rz(-0.31983917) q[1];
rz(2.1460974) q[3];
sx q[3];
rz(-0.18333215) q[3];
sx q[3];
rz(-3.0005279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.60160294) q[2];
sx q[2];
rz(-1.8113965) q[2];
sx q[2];
rz(0.23923624) q[2];
rz(0.38419497) q[3];
sx q[3];
rz(-2.4592063) q[3];
sx q[3];
rz(0.95170963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20188986) q[0];
sx q[0];
rz(-1.6738142) q[0];
sx q[0];
rz(-1.7643167) q[0];
rz(-1.2205623) q[1];
sx q[1];
rz(-2.2790597) q[1];
sx q[1];
rz(-2.9753704) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8140504) q[0];
sx q[0];
rz(-0.53314994) q[0];
sx q[0];
rz(2.0943805) q[0];
x q[1];
rz(-0.73428369) q[2];
sx q[2];
rz(-1.3314651) q[2];
sx q[2];
rz(0.83877968) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.080648153) q[1];
sx q[1];
rz(-1.540144) q[1];
sx q[1];
rz(2.9303307) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1020426) q[3];
sx q[3];
rz(-0.47787468) q[3];
sx q[3];
rz(-2.9210764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.347747) q[2];
sx q[2];
rz(-0.85024992) q[2];
sx q[2];
rz(-2.094685) q[2];
rz(-1.8868123) q[3];
sx q[3];
rz(-2.051765) q[3];
sx q[3];
rz(-1.2966398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5407402) q[0];
sx q[0];
rz(-2.7678601) q[0];
sx q[0];
rz(-1.1131713) q[0];
rz(-0.81740776) q[1];
sx q[1];
rz(-1.4459041) q[1];
sx q[1];
rz(-2.7730952) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8857074) q[0];
sx q[0];
rz(-1.2987381) q[0];
sx q[0];
rz(-1.2124802) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0037874) q[2];
sx q[2];
rz(-2.3944582) q[2];
sx q[2];
rz(1.8835619) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8702183) q[1];
sx q[1];
rz(-1.8370359) q[1];
sx q[1];
rz(-2.4143747) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0115252) q[3];
sx q[3];
rz(-0.91010053) q[3];
sx q[3];
rz(-1.6370156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.14785279) q[2];
sx q[2];
rz(-0.61823121) q[2];
sx q[2];
rz(-1.3573793) q[2];
rz(-2.3944858) q[3];
sx q[3];
rz(-2.1785469) q[3];
sx q[3];
rz(-0.67556206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6152182) q[0];
sx q[0];
rz(-0.45658657) q[0];
sx q[0];
rz(2.0988462) q[0];
rz(0.77049795) q[1];
sx q[1];
rz(-1.0105402) q[1];
sx q[1];
rz(-2.3811293) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7211518) q[0];
sx q[0];
rz(-1.8095224) q[0];
sx q[0];
rz(-1.9107642) q[0];
x q[1];
rz(1.3622401) q[2];
sx q[2];
rz(-2.957666) q[2];
sx q[2];
rz(0.49801258) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9925966) q[1];
sx q[1];
rz(-1.1482052) q[1];
sx q[1];
rz(-1.6572957) q[1];
x q[2];
rz(-2.4861927) q[3];
sx q[3];
rz(-1.8701347) q[3];
sx q[3];
rz(1.4531524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9383135) q[2];
sx q[2];
rz(-0.61050296) q[2];
sx q[2];
rz(-2.7321613) q[2];
rz(-1.5832541) q[3];
sx q[3];
rz(-2.6810472) q[3];
sx q[3];
rz(0.62622768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4938477) q[0];
sx q[0];
rz(-2.0997601) q[0];
sx q[0];
rz(0.44152015) q[0];
rz(-1.582675) q[1];
sx q[1];
rz(-1.9428923) q[1];
sx q[1];
rz(2.518242) q[1];
rz(1.5123488) q[2];
sx q[2];
rz(-1.4679906) q[2];
sx q[2];
rz(1.0124258) q[2];
rz(-0.49420301) q[3];
sx q[3];
rz(-0.27402607) q[3];
sx q[3];
rz(-1.5370697) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
