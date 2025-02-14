OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6899941) q[0];
sx q[0];
rz(-2.8328083) q[0];
sx q[0];
rz(0.2398332) q[0];
rz(2.580515) q[1];
sx q[1];
rz(-0.74220389) q[1];
sx q[1];
rz(1.5129369) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5133249) q[0];
sx q[0];
rz(-1.6986934) q[0];
sx q[0];
rz(1.9403752) q[0];
rz(1.0983019) q[2];
sx q[2];
rz(-1.3368946) q[2];
sx q[2];
rz(-0.59532524) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.35084942) q[1];
sx q[1];
rz(-0.51209699) q[1];
sx q[1];
rz(-2.8382906) q[1];
x q[2];
rz(0.058006279) q[3];
sx q[3];
rz(-1.4027624) q[3];
sx q[3];
rz(-1.9498273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.98770398) q[2];
sx q[2];
rz(-2.14553) q[2];
sx q[2];
rz(1.055701) q[2];
rz(-2.740247) q[3];
sx q[3];
rz(-1.6258806) q[3];
sx q[3];
rz(-1.5041941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5098679) q[0];
sx q[0];
rz(-0.26917502) q[0];
sx q[0];
rz(-1.7357695) q[0];
rz(2.3954605) q[1];
sx q[1];
rz(-1.1765307) q[1];
sx q[1];
rz(-0.064780898) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21532741) q[0];
sx q[0];
rz(-2.522922) q[0];
sx q[0];
rz(-2.4149405) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4006813) q[2];
sx q[2];
rz(-0.89808849) q[2];
sx q[2];
rz(1.9667039) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.83103115) q[1];
sx q[1];
rz(-1.6509202) q[1];
sx q[1];
rz(1.7145654) q[1];
rz(1.359932) q[3];
sx q[3];
rz(-1.5393889) q[3];
sx q[3];
rz(-1.5109911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.259321) q[2];
sx q[2];
rz(-0.52053014) q[2];
sx q[2];
rz(0.041291324) q[2];
rz(-2.0222372) q[3];
sx q[3];
rz(-1.3675523) q[3];
sx q[3];
rz(-1.1411427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3248046) q[0];
sx q[0];
rz(-0.4275221) q[0];
sx q[0];
rz(0.69931716) q[0];
rz(-0.62272561) q[1];
sx q[1];
rz(-2.3284349) q[1];
sx q[1];
rz(2.8831388) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7197882) q[0];
sx q[0];
rz(-2.5076619) q[0];
sx q[0];
rz(3.1209206) q[0];
x q[1];
rz(-1.4213561) q[2];
sx q[2];
rz(-0.65776134) q[2];
sx q[2];
rz(0.8873111) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.093100637) q[1];
sx q[1];
rz(-2.6922467) q[1];
sx q[1];
rz(-3.0952342) q[1];
x q[2];
rz(-1.4914037) q[3];
sx q[3];
rz(-2.7104122) q[3];
sx q[3];
rz(3.0357547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2059325) q[2];
sx q[2];
rz(-1.460133) q[2];
sx q[2];
rz(2.4350731) q[2];
rz(2.5126854) q[3];
sx q[3];
rz(-0.8258515) q[3];
sx q[3];
rz(1.1227054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0141456) q[0];
sx q[0];
rz(-1.7147467) q[0];
sx q[0];
rz(-2.1674147) q[0];
rz(1.4840508) q[1];
sx q[1];
rz(-1.0933135) q[1];
sx q[1];
rz(-2.1407703) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38353048) q[0];
sx q[0];
rz(-2.5208726) q[0];
sx q[0];
rz(-1.3002943) q[0];
x q[1];
rz(0.77299574) q[2];
sx q[2];
rz(-2.2064035) q[2];
sx q[2];
rz(2.5781812) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.82881935) q[1];
sx q[1];
rz(-2.3760894) q[1];
sx q[1];
rz(2.440052) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6203247) q[3];
sx q[3];
rz(-1.2504745) q[3];
sx q[3];
rz(-1.76521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.68690825) q[2];
sx q[2];
rz(-1.3409706) q[2];
sx q[2];
rz(-2.3243813) q[2];
rz(2.5456083) q[3];
sx q[3];
rz(-1.417336) q[3];
sx q[3];
rz(-1.7433085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3572094) q[0];
sx q[0];
rz(-1.4056118) q[0];
sx q[0];
rz(2.0696409) q[0];
rz(1.1373854) q[1];
sx q[1];
rz(-1.6268566) q[1];
sx q[1];
rz(-0.17328182) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3629775) q[0];
sx q[0];
rz(-2.5234875) q[0];
sx q[0];
rz(-0.72160665) q[0];
rz(1.0801804) q[2];
sx q[2];
rz(-1.2254493) q[2];
sx q[2];
rz(2.5290348) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.73663354) q[1];
sx q[1];
rz(-1.6425321) q[1];
sx q[1];
rz(-0.81291764) q[1];
x q[2];
rz(1.4777484) q[3];
sx q[3];
rz(-1.4085839) q[3];
sx q[3];
rz(1.4402267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5002354) q[2];
sx q[2];
rz(-2.6884029) q[2];
sx q[2];
rz(2.267061) q[2];
rz(-0.66679653) q[3];
sx q[3];
rz(-1.4614481) q[3];
sx q[3];
rz(-0.82505208) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2905529) q[0];
sx q[0];
rz(-2.9930826) q[0];
sx q[0];
rz(-2.9669951) q[0];
rz(-2.5945276) q[1];
sx q[1];
rz(-2.6262296) q[1];
sx q[1];
rz(-1.863716) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39659835) q[0];
sx q[0];
rz(-1.5088668) q[0];
sx q[0];
rz(-1.8577736) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1013979) q[2];
sx q[2];
rz(-2.2661327) q[2];
sx q[2];
rz(2.4027195) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4365579) q[1];
sx q[1];
rz(-0.83022699) q[1];
sx q[1];
rz(0.48506801) q[1];
rz(-2.296769) q[3];
sx q[3];
rz(-0.81171821) q[3];
sx q[3];
rz(-2.7561848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7586907) q[2];
sx q[2];
rz(-2.8768657) q[2];
sx q[2];
rz(-0.36941377) q[2];
rz(0.82768011) q[3];
sx q[3];
rz(-1.8281432) q[3];
sx q[3];
rz(2.9874492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5740042) q[0];
sx q[0];
rz(-0.91667691) q[0];
sx q[0];
rz(-0.23707238) q[0];
rz(1.392662) q[1];
sx q[1];
rz(-1.9475513) q[1];
sx q[1];
rz(1.0038092) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061351731) q[0];
sx q[0];
rz(-2.6727242) q[0];
sx q[0];
rz(-0.37567015) q[0];
rz(-1.1091484) q[2];
sx q[2];
rz(-0.91141111) q[2];
sx q[2];
rz(-0.75379363) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8449515) q[1];
sx q[1];
rz(-1.2403204) q[1];
sx q[1];
rz(0.56092324) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0521734) q[3];
sx q[3];
rz(-2.0767835) q[3];
sx q[3];
rz(-0.86694366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6134593) q[2];
sx q[2];
rz(-1.5789092) q[2];
sx q[2];
rz(-0.51631874) q[2];
rz(-0.83827072) q[3];
sx q[3];
rz(-1.4811938) q[3];
sx q[3];
rz(2.9576438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2876005) q[0];
sx q[0];
rz(-0.28165278) q[0];
sx q[0];
rz(0.91424346) q[0];
rz(2.5573348) q[1];
sx q[1];
rz(-1.4062358) q[1];
sx q[1];
rz(-0.90726888) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6833008) q[0];
sx q[0];
rz(-2.2193546) q[0];
sx q[0];
rz(0.022819937) q[0];
rz(-pi) q[1];
rz(-1.7444939) q[2];
sx q[2];
rz(-0.96250421) q[2];
sx q[2];
rz(2.4243958) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4248391) q[1];
sx q[1];
rz(-2.3100634) q[1];
sx q[1];
rz(-2.5118973) q[1];
rz(-pi) q[2];
rz(2.1694555) q[3];
sx q[3];
rz(-1.3479509) q[3];
sx q[3];
rz(-0.11892648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3153136) q[2];
sx q[2];
rz(-1.4297012) q[2];
sx q[2];
rz(-0.86177525) q[2];
rz(-1.2365384) q[3];
sx q[3];
rz(-3.0573513) q[3];
sx q[3];
rz(-1.2110565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5009907) q[0];
sx q[0];
rz(-0.7826829) q[0];
sx q[0];
rz(-1.030141) q[0];
rz(0.31632272) q[1];
sx q[1];
rz(-1.373469) q[1];
sx q[1];
rz(-2.8533459) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.41053) q[0];
sx q[0];
rz(-1.1708492) q[0];
sx q[0];
rz(-1.4131141) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4924303) q[2];
sx q[2];
rz(-1.1466951) q[2];
sx q[2];
rz(2.8555388) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1128149) q[1];
sx q[1];
rz(-2.3773686) q[1];
sx q[1];
rz(-2.1995596) q[1];
rz(2.8602198) q[3];
sx q[3];
rz(-2.4503539) q[3];
sx q[3];
rz(-0.465525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2890702) q[2];
sx q[2];
rz(-1.4129637) q[2];
sx q[2];
rz(-2.9150325) q[2];
rz(1.4551) q[3];
sx q[3];
rz(-2.2454567) q[3];
sx q[3];
rz(-2.4494825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9914472) q[0];
sx q[0];
rz(-1.2757855) q[0];
sx q[0];
rz(-0.62514296) q[0];
rz(-1.217968) q[1];
sx q[1];
rz(-0.70911276) q[1];
sx q[1];
rz(-1.9295173) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3071412) q[0];
sx q[0];
rz(-2.0491776) q[0];
sx q[0];
rz(-0.68104736) q[0];
rz(2.8202989) q[2];
sx q[2];
rz(-0.69707631) q[2];
sx q[2];
rz(-1.4980121) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9245476) q[1];
sx q[1];
rz(-0.44483652) q[1];
sx q[1];
rz(1.43631) q[1];
x q[2];
rz(-2.3732568) q[3];
sx q[3];
rz(-0.4076805) q[3];
sx q[3];
rz(2.3967495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7071699) q[2];
sx q[2];
rz(-1.3753128) q[2];
sx q[2];
rz(-1.0506857) q[2];
rz(0.55189842) q[3];
sx q[3];
rz(-2.4514908) q[3];
sx q[3];
rz(1.8570073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62591775) q[0];
sx q[0];
rz(-0.82763012) q[0];
sx q[0];
rz(-0.81830842) q[0];
rz(-0.7863518) q[1];
sx q[1];
rz(-0.66023371) q[1];
sx q[1];
rz(0.22088851) q[1];
rz(-0.037842265) q[2];
sx q[2];
rz(-1.2076245) q[2];
sx q[2];
rz(-2.5943499) q[2];
rz(1.4825357) q[3];
sx q[3];
rz(-0.88659989) q[3];
sx q[3];
rz(-2.551589) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
