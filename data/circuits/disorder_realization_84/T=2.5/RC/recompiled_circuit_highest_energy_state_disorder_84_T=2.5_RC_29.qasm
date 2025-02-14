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
rz(-0.62491971) q[0];
sx q[0];
rz(4.9486296) q[0];
sx q[0];
rz(9.3715342) q[0];
rz(0.48625311) q[1];
sx q[1];
rz(-0.12737218) q[1];
sx q[1];
rz(1.706634) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.931728) q[0];
sx q[0];
rz(-1.7057944) q[0];
sx q[0];
rz(-0.28688669) q[0];
x q[1];
rz(1.6900914) q[2];
sx q[2];
rz(-1.3053467) q[2];
sx q[2];
rz(-0.99670289) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9588152) q[1];
sx q[1];
rz(-0.39740409) q[1];
sx q[1];
rz(-0.13840492) q[1];
rz(-2.1987183) q[3];
sx q[3];
rz(-2.3981895) q[3];
sx q[3];
rz(-1.8441895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3367553) q[2];
sx q[2];
rz(-1.3433604) q[2];
sx q[2];
rz(-0.27377823) q[2];
rz(-1.9233507) q[3];
sx q[3];
rz(-2.4680586) q[3];
sx q[3];
rz(0.28606733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022920595) q[0];
sx q[0];
rz(-2.3781222) q[0];
sx q[0];
rz(-0.77962312) q[0];
rz(2.4769705) q[1];
sx q[1];
rz(-1.9326991) q[1];
sx q[1];
rz(1.8962616) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1561403) q[0];
sx q[0];
rz(-0.13053556) q[0];
sx q[0];
rz(-1.077561) q[0];
rz(-pi) q[1];
rz(-1.4636366) q[2];
sx q[2];
rz(-0.71427155) q[2];
sx q[2];
rz(-0.34253866) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9960963) q[1];
sx q[1];
rz(-2.9206373) q[1];
sx q[1];
rz(-1.851382) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8923442) q[3];
sx q[3];
rz(-2.5587497) q[3];
sx q[3];
rz(3.1414523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.43210426) q[2];
sx q[2];
rz(-2.4435142) q[2];
sx q[2];
rz(-3.091605) q[2];
rz(-1.4738119) q[3];
sx q[3];
rz(-1.4155017) q[3];
sx q[3];
rz(-0.93562359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98899984) q[0];
sx q[0];
rz(-0.86243668) q[0];
sx q[0];
rz(-1.9631901) q[0];
rz(1.9940469) q[1];
sx q[1];
rz(-0.18745628) q[1];
sx q[1];
rz(2.5386834) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1499304) q[0];
sx q[0];
rz(-1.4105182) q[0];
sx q[0];
rz(3.0383238) q[0];
rz(-pi) q[1];
rz(0.52003543) q[2];
sx q[2];
rz(-0.94920659) q[2];
sx q[2];
rz(0.015403143) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5192831) q[1];
sx q[1];
rz(-1.9819248) q[1];
sx q[1];
rz(-2.4101188) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6125978) q[3];
sx q[3];
rz(-0.96312614) q[3];
sx q[3];
rz(-1.120847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19043645) q[2];
sx q[2];
rz(-2.3556605) q[2];
sx q[2];
rz(-2.7583165) q[2];
rz(1.3112618) q[3];
sx q[3];
rz(-2.0537328) q[3];
sx q[3];
rz(-3.0200628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.424778) q[0];
sx q[0];
rz(-0.99262339) q[0];
sx q[0];
rz(-2.2130261) q[0];
rz(2.3021452) q[1];
sx q[1];
rz(-1.3176368) q[1];
sx q[1];
rz(-2.3323257) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1195212) q[0];
sx q[0];
rz(-0.68953994) q[0];
sx q[0];
rz(0.9718231) q[0];
x q[1];
rz(-1.7401314) q[2];
sx q[2];
rz(-1.9953097) q[2];
sx q[2];
rz(-2.9237539) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7245771) q[1];
sx q[1];
rz(-2.8863686) q[1];
sx q[1];
rz(-1.8468813) q[1];
rz(-pi) q[2];
rz(1.0471116) q[3];
sx q[3];
rz(-0.71412702) q[3];
sx q[3];
rz(-1.7585825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3889918) q[2];
sx q[2];
rz(-1.5316803) q[2];
sx q[2];
rz(1.6947702) q[2];
rz(-2.0145448) q[3];
sx q[3];
rz(-0.73486745) q[3];
sx q[3];
rz(-2.5126422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9160354) q[0];
sx q[0];
rz(-1.1981523) q[0];
sx q[0];
rz(1.7115364) q[0];
rz(-0.23722181) q[1];
sx q[1];
rz(-0.96313852) q[1];
sx q[1];
rz(-0.29475862) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98210456) q[0];
sx q[0];
rz(-2.0061261) q[0];
sx q[0];
rz(1.728471) q[0];
rz(-pi) q[1];
rz(2.2331401) q[2];
sx q[2];
rz(-1.2725432) q[2];
sx q[2];
rz(-0.82240381) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.393906) q[1];
sx q[1];
rz(-1.1562041) q[1];
sx q[1];
rz(1.0769597) q[1];
rz(-pi) q[2];
rz(0.57394694) q[3];
sx q[3];
rz(-1.4455631) q[3];
sx q[3];
rz(0.67170152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99990591) q[2];
sx q[2];
rz(-0.20076951) q[2];
sx q[2];
rz(-3.0086009) q[2];
rz(-0.76987949) q[3];
sx q[3];
rz(-1.8498288) q[3];
sx q[3];
rz(2.8225115) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8517476) q[0];
sx q[0];
rz(-0.3388437) q[0];
sx q[0];
rz(-1.4122562) q[0];
rz(0.051636592) q[1];
sx q[1];
rz(-2.5187571) q[1];
sx q[1];
rz(2.9488865) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4225005) q[0];
sx q[0];
rz(-2.0845045) q[0];
sx q[0];
rz(-1.5490319) q[0];
rz(-pi) q[1];
rz(0.66484837) q[2];
sx q[2];
rz(-2.2908205) q[2];
sx q[2];
rz(-2.8262422) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.440408) q[1];
sx q[1];
rz(-1.1478979) q[1];
sx q[1];
rz(-2.6404523) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9514378) q[3];
sx q[3];
rz(-2.1771113) q[3];
sx q[3];
rz(-2.9553138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0685588) q[2];
sx q[2];
rz(-1.2898338) q[2];
sx q[2];
rz(1.3812836) q[2];
rz(-0.51501385) q[3];
sx q[3];
rz(-0.88740715) q[3];
sx q[3];
rz(1.380434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16159049) q[0];
sx q[0];
rz(-1.1450293) q[0];
sx q[0];
rz(0.34647754) q[0];
rz(1.0193635) q[1];
sx q[1];
rz(-0.66277021) q[1];
sx q[1];
rz(-1.6874708) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2812735) q[0];
sx q[0];
rz(-2.3695282) q[0];
sx q[0];
rz(-0.66769974) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90482651) q[2];
sx q[2];
rz(-2.3695393) q[2];
sx q[2];
rz(-2.0029298) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2994308) q[1];
sx q[1];
rz(-2.5379532) q[1];
sx q[1];
rz(1.0697212) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4164657) q[3];
sx q[3];
rz(-1.6701506) q[3];
sx q[3];
rz(-1.1442483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5399897) q[2];
sx q[2];
rz(-1.8113965) q[2];
sx q[2];
rz(0.23923624) q[2];
rz(2.7573977) q[3];
sx q[3];
rz(-0.68238634) q[3];
sx q[3];
rz(0.95170963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9397028) q[0];
sx q[0];
rz(-1.6738142) q[0];
sx q[0];
rz(-1.3772759) q[0];
rz(1.2205623) q[1];
sx q[1];
rz(-2.2790597) q[1];
sx q[1];
rz(2.9753704) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2234874) q[0];
sx q[0];
rz(-2.0265731) q[0];
sx q[0];
rz(-2.8546643) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.407309) q[2];
sx q[2];
rz(-1.8101276) q[2];
sx q[2];
rz(-2.302813) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6580171) q[1];
sx q[1];
rz(-1.7819575) q[1];
sx q[1];
rz(-1.5394475) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0395501) q[3];
sx q[3];
rz(-2.663718) q[3];
sx q[3];
rz(-0.22051624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.347747) q[2];
sx q[2];
rz(-2.2913427) q[2];
sx q[2];
rz(-2.094685) q[2];
rz(-1.2547803) q[3];
sx q[3];
rz(-1.0898277) q[3];
sx q[3];
rz(-1.2966398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5407402) q[0];
sx q[0];
rz(-0.37373251) q[0];
sx q[0];
rz(-2.0284213) q[0];
rz(0.81740776) q[1];
sx q[1];
rz(-1.4459041) q[1];
sx q[1];
rz(2.7730952) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2042646) q[0];
sx q[0];
rz(-2.6952792) q[0];
sx q[0];
rz(0.89881368) q[0];
x q[1];
rz(3.0037874) q[2];
sx q[2];
rz(-0.74713444) q[2];
sx q[2];
rz(-1.2580308) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6121713) q[1];
sx q[1];
rz(-2.2670806) q[1];
sx q[1];
rz(1.9208292) q[1];
rz(1.0115252) q[3];
sx q[3];
rz(-2.2314921) q[3];
sx q[3];
rz(1.504577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.14785279) q[2];
sx q[2];
rz(-2.5233614) q[2];
sx q[2];
rz(1.7842133) q[2];
rz(0.74710685) q[3];
sx q[3];
rz(-0.96304572) q[3];
sx q[3];
rz(0.67556206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6152182) q[0];
sx q[0];
rz(-2.6850061) q[0];
sx q[0];
rz(-1.0427465) q[0];
rz(0.77049795) q[1];
sx q[1];
rz(-2.1310525) q[1];
sx q[1];
rz(-0.76046336) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5805394) q[0];
sx q[0];
rz(-2.7288611) q[0];
sx q[0];
rz(-0.94032918) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7508239) q[2];
sx q[2];
rz(-1.6086726) q[2];
sx q[2];
rz(1.2779209) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.14899602) q[1];
sx q[1];
rz(-1.9933874) q[1];
sx q[1];
rz(1.6572957) q[1];
rz(-pi) q[2];
rz(-0.46871878) q[3];
sx q[3];
rz(-0.7112452) q[3];
sx q[3];
rz(2.6576692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-2.515365) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4938477) q[0];
sx q[0];
rz(-2.0997601) q[0];
sx q[0];
rz(0.44152015) q[0];
rz(-1.5589177) q[1];
sx q[1];
rz(-1.1987004) q[1];
sx q[1];
rz(-0.62335062) q[1];
rz(2.6263993) q[2];
sx q[2];
rz(-3.0233848) q[2];
sx q[2];
rz(1.530627) q[2];
rz(1.7033475) q[3];
sx q[3];
rz(-1.3302531) q[3];
sx q[3];
rz(2.1147685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
