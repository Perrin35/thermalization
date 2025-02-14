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
rz(-0.15859088) q[0];
sx q[0];
rz(-0.51704419) q[0];
sx q[0];
rz(1.3102732) q[0];
rz(-0.34215555) q[1];
sx q[1];
rz(-1.181239) q[1];
sx q[1];
rz(0.96460834) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2268596) q[0];
sx q[0];
rz(-0.64382416) q[0];
sx q[0];
rz(-2.9152246) q[0];
x q[1];
rz(-0.50268051) q[2];
sx q[2];
rz(-2.7516904) q[2];
sx q[2];
rz(-1.0165868) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1312932) q[1];
sx q[1];
rz(-2.2939766) q[1];
sx q[1];
rz(-2.7406969) q[1];
rz(0.89972605) q[3];
sx q[3];
rz(-1.6613071) q[3];
sx q[3];
rz(-2.8495726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1356807) q[2];
sx q[2];
rz(-2.9475806) q[2];
sx q[2];
rz(-1.4646336) q[2];
rz(-1.8320734) q[3];
sx q[3];
rz(-0.94683164) q[3];
sx q[3];
rz(0.32002282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.2784073) q[0];
sx q[0];
rz(-2.0018556) q[0];
sx q[0];
rz(-1.7742668) q[0];
rz(1.6525846) q[1];
sx q[1];
rz(-0.79843489) q[1];
sx q[1];
rz(2.8489825) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3037138) q[0];
sx q[0];
rz(-2.7927164) q[0];
sx q[0];
rz(2.8512591) q[0];
rz(-1.4433447) q[2];
sx q[2];
rz(-1.410219) q[2];
sx q[2];
rz(1.0739971) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.039174883) q[1];
sx q[1];
rz(-2.8264649) q[1];
sx q[1];
rz(-0.71690317) q[1];
x q[2];
rz(2.8752453) q[3];
sx q[3];
rz(-1.6345336) q[3];
sx q[3];
rz(-1.6132465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1613529) q[2];
sx q[2];
rz(-0.62675256) q[2];
sx q[2];
rz(-0.18388595) q[2];
rz(0.45267496) q[3];
sx q[3];
rz(-1.5225182) q[3];
sx q[3];
rz(-3.0116853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1995131) q[0];
sx q[0];
rz(-0.71490723) q[0];
sx q[0];
rz(-2.3364501) q[0];
rz(1.0792271) q[1];
sx q[1];
rz(-2.3715623) q[1];
sx q[1];
rz(1.8260746) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.091324) q[0];
sx q[0];
rz(-1.6299695) q[0];
sx q[0];
rz(-2.4135804) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82683021) q[2];
sx q[2];
rz(-1.2203802) q[2];
sx q[2];
rz(2.0598799) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.47655573) q[1];
sx q[1];
rz(-1.2893234) q[1];
sx q[1];
rz(-0.095515619) q[1];
x q[2];
rz(2.838921) q[3];
sx q[3];
rz(-1.4883452) q[3];
sx q[3];
rz(-0.77471212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9341854) q[2];
sx q[2];
rz(-1.4855874) q[2];
sx q[2];
rz(2.5918813) q[2];
rz(3.1304729) q[3];
sx q[3];
rz(-2.34237) q[3];
sx q[3];
rz(-1.3219249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84871197) q[0];
sx q[0];
rz(-0.95893812) q[0];
sx q[0];
rz(2.8380561) q[0];
rz(-2.983298) q[1];
sx q[1];
rz(-1.0946495) q[1];
sx q[1];
rz(2.1319481) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5490259) q[0];
sx q[0];
rz(-2.0011304) q[0];
sx q[0];
rz(0.31220147) q[0];
x q[1];
rz(-2.7042806) q[2];
sx q[2];
rz(-0.6491937) q[2];
sx q[2];
rz(-1.6216261) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9480966) q[1];
sx q[1];
rz(-1.0973012) q[1];
sx q[1];
rz(-1.7481453) q[1];
rz(0.03831717) q[3];
sx q[3];
rz(-2.6434757) q[3];
sx q[3];
rz(1.2850873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2046854) q[2];
sx q[2];
rz(-2.3202809) q[2];
sx q[2];
rz(-0.26910195) q[2];
rz(1.6875632) q[3];
sx q[3];
rz(-1.5498091) q[3];
sx q[3];
rz(1.5627741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.27075416) q[0];
sx q[0];
rz(-1.0601059) q[0];
sx q[0];
rz(-2.7794072) q[0];
rz(-0.65131342) q[1];
sx q[1];
rz(-1.4910699) q[1];
sx q[1];
rz(-1.6168894) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0391738) q[0];
sx q[0];
rz(-0.96511594) q[0];
sx q[0];
rz(-1.7282979) q[0];
rz(-pi) q[1];
rz(-0.082090898) q[2];
sx q[2];
rz(-2.1920648) q[2];
sx q[2];
rz(2.4022849) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.32423577) q[1];
sx q[1];
rz(-0.65297665) q[1];
sx q[1];
rz(-0.95545902) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7157341) q[3];
sx q[3];
rz(-0.79090624) q[3];
sx q[3];
rz(0.87825852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0321956) q[2];
sx q[2];
rz(-0.95462644) q[2];
sx q[2];
rz(-1.7677914) q[2];
rz(0.43921709) q[3];
sx q[3];
rz(-2.4924811) q[3];
sx q[3];
rz(-0.60905987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1010901) q[0];
sx q[0];
rz(-0.73223615) q[0];
sx q[0];
rz(-0.26443431) q[0];
rz(2.4624372) q[1];
sx q[1];
rz(-0.87398386) q[1];
sx q[1];
rz(2.1566379) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0593582) q[0];
sx q[0];
rz(-1.2684039) q[0];
sx q[0];
rz(-1.9134269) q[0];
x q[1];
rz(-0.28552766) q[2];
sx q[2];
rz(-2.0047024) q[2];
sx q[2];
rz(-1.7800231) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86414908) q[1];
sx q[1];
rz(-1.207749) q[1];
sx q[1];
rz(-2.6226642) q[1];
x q[2];
rz(1.9846647) q[3];
sx q[3];
rz(-0.47734161) q[3];
sx q[3];
rz(-1.1787924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.039845) q[2];
sx q[2];
rz(-2.9797649) q[2];
sx q[2];
rz(-0.43231535) q[2];
rz(2.1281706) q[3];
sx q[3];
rz(-1.6067959) q[3];
sx q[3];
rz(-2.5920946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6559615) q[0];
sx q[0];
rz(-0.39325842) q[0];
sx q[0];
rz(-1.1095169) q[0];
rz(-2.3484777) q[1];
sx q[1];
rz(-1.7879281) q[1];
sx q[1];
rz(-1.5464334) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4684179) q[0];
sx q[0];
rz(-2.1270848) q[0];
sx q[0];
rz(-2.2493717) q[0];
x q[1];
rz(2.4199615) q[2];
sx q[2];
rz(-1.7296556) q[2];
sx q[2];
rz(1.5312486) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2364051) q[1];
sx q[1];
rz(-1.3705472) q[1];
sx q[1];
rz(0.88242857) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71641201) q[3];
sx q[3];
rz(-1.4906814) q[3];
sx q[3];
rz(0.83601213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.82211295) q[2];
sx q[2];
rz(-0.18222465) q[2];
sx q[2];
rz(-1.7721843) q[2];
rz(2.2111514) q[3];
sx q[3];
rz(-1.7934099) q[3];
sx q[3];
rz(-1.5038917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1996138) q[0];
sx q[0];
rz(-1.1811341) q[0];
sx q[0];
rz(0.30366316) q[0];
rz(-0.97967255) q[1];
sx q[1];
rz(-0.62398463) q[1];
sx q[1];
rz(-2.7365275) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9054149) q[0];
sx q[0];
rz(-1.9472031) q[0];
sx q[0];
rz(-1.3628528) q[0];
rz(-pi) q[1];
rz(2.9689661) q[2];
sx q[2];
rz(-1.1571047) q[2];
sx q[2];
rz(2.911441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.67928708) q[1];
sx q[1];
rz(-2.132173) q[1];
sx q[1];
rz(1.7803468) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9977244) q[3];
sx q[3];
rz(-1.3865304) q[3];
sx q[3];
rz(-1.1729405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27641174) q[2];
sx q[2];
rz(-2.1750735) q[2];
sx q[2];
rz(0.69918862) q[2];
rz(0.80896038) q[3];
sx q[3];
rz(-1.8374846) q[3];
sx q[3];
rz(-2.8185524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-2.9998099) q[0];
sx q[0];
rz(-1.1572105) q[0];
sx q[0];
rz(3.004177) q[0];
rz(0.049086463) q[1];
sx q[1];
rz(-1.1156491) q[1];
sx q[1];
rz(0.5079937) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79530966) q[0];
sx q[0];
rz(-2.6457743) q[0];
sx q[0];
rz(-3.0853372) q[0];
rz(-2.9692978) q[2];
sx q[2];
rz(-1.5827547) q[2];
sx q[2];
rz(-1.1520141) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0145871) q[1];
sx q[1];
rz(-1.0795322) q[1];
sx q[1];
rz(0.36775695) q[1];
rz(-pi) q[2];
rz(-1.6925595) q[3];
sx q[3];
rz(-1.6596438) q[3];
sx q[3];
rz(2.4112687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5892443) q[2];
sx q[2];
rz(-1.2184703) q[2];
sx q[2];
rz(-0.98769665) q[2];
rz(-2.6622631) q[3];
sx q[3];
rz(-1.5440116) q[3];
sx q[3];
rz(2.2492669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.669303) q[0];
sx q[0];
rz(-2.9030114) q[0];
sx q[0];
rz(-0.14078374) q[0];
rz(0.36551481) q[1];
sx q[1];
rz(-2.3194158) q[1];
sx q[1];
rz(0.64868322) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9539584) q[0];
sx q[0];
rz(-1.2231069) q[0];
sx q[0];
rz(2.4566922) q[0];
rz(2.3769605) q[2];
sx q[2];
rz(-1.4864385) q[2];
sx q[2];
rz(0.81071883) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2461697) q[1];
sx q[1];
rz(-1.6270346) q[1];
sx q[1];
rz(0.59631056) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8148531) q[3];
sx q[3];
rz(-1.3413197) q[3];
sx q[3];
rz(-0.47577259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6792128) q[2];
sx q[2];
rz(-2.1361735) q[2];
sx q[2];
rz(3.0688378) q[2];
rz(-2.4788729) q[3];
sx q[3];
rz(-1.8279671) q[3];
sx q[3];
rz(-0.16547468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0872021) q[0];
sx q[0];
rz(-1.5329755) q[0];
sx q[0];
rz(-1.6736915) q[0];
rz(-1.6064593) q[1];
sx q[1];
rz(-0.88773334) q[1];
sx q[1];
rz(1.6233374) q[1];
rz(0.22850484) q[2];
sx q[2];
rz(-1.6863556) q[2];
sx q[2];
rz(-0.40876331) q[2];
rz(0.40409186) q[3];
sx q[3];
rz(-0.93032144) q[3];
sx q[3];
rz(-3.0138187) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
