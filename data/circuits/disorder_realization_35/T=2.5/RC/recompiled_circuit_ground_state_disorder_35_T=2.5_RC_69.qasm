OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93249455) q[0];
sx q[0];
rz(-0.068109186) q[0];
sx q[0];
rz(0.77686247) q[0];
rz(1.8344954) q[1];
sx q[1];
rz(-2.1721462) q[1];
sx q[1];
rz(-1.3219272) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37047077) q[0];
sx q[0];
rz(-1.4323455) q[0];
sx q[0];
rz(0.8531424) q[0];
rz(2.5547682) q[2];
sx q[2];
rz(-0.57673645) q[2];
sx q[2];
rz(1.5962034) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.81404954) q[1];
sx q[1];
rz(-2.7316878) q[1];
sx q[1];
rz(1.4973212) q[1];
x q[2];
rz(-2.1573651) q[3];
sx q[3];
rz(-2.7968458) q[3];
sx q[3];
rz(-1.00351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0778568) q[2];
sx q[2];
rz(-1.7813762) q[2];
sx q[2];
rz(0.35749164) q[2];
rz(2.5743971) q[3];
sx q[3];
rz(-1.7287858) q[3];
sx q[3];
rz(2.7020057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1883063) q[0];
sx q[0];
rz(-0.28306857) q[0];
sx q[0];
rz(1.333492) q[0];
rz(-0.4370583) q[1];
sx q[1];
rz(-0.60085618) q[1];
sx q[1];
rz(-2.3016047) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14535689) q[0];
sx q[0];
rz(-1.0635492) q[0];
sx q[0];
rz(-0.17542393) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8905021) q[2];
sx q[2];
rz(-0.60949627) q[2];
sx q[2];
rz(1.9025546) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31038302) q[1];
sx q[1];
rz(-1.0070756) q[1];
sx q[1];
rz(-0.68117627) q[1];
rz(-pi) q[2];
rz(-1.6343104) q[3];
sx q[3];
rz(-0.5507142) q[3];
sx q[3];
rz(-2.7515819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4240894) q[2];
sx q[2];
rz(-1.4428029) q[2];
sx q[2];
rz(1.7355512) q[2];
rz(-2.9794335) q[3];
sx q[3];
rz(-1.5806961) q[3];
sx q[3];
rz(0.46419188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9185987) q[0];
sx q[0];
rz(-0.080568947) q[0];
sx q[0];
rz(-0.42174569) q[0];
rz(-2.2987507) q[1];
sx q[1];
rz(-1.755736) q[1];
sx q[1];
rz(0.27580321) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0378135) q[0];
sx q[0];
rz(-2.7386754) q[0];
sx q[0];
rz(-0.62312868) q[0];
rz(-2.7245443) q[2];
sx q[2];
rz(-0.50784238) q[2];
sx q[2];
rz(1.5323973) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7843102) q[1];
sx q[1];
rz(-1.5187335) q[1];
sx q[1];
rz(-0.32084685) q[1];
x q[2];
rz(0.098297997) q[3];
sx q[3];
rz(-0.87140981) q[3];
sx q[3];
rz(0.15546945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80321035) q[2];
sx q[2];
rz(-2.3440177) q[2];
sx q[2];
rz(-0.69022834) q[2];
rz(2.7817182) q[3];
sx q[3];
rz(-1.7706324) q[3];
sx q[3];
rz(-2.7580822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2279219) q[0];
sx q[0];
rz(-2.0991195) q[0];
sx q[0];
rz(0.87164718) q[0];
rz(-2.7323885) q[1];
sx q[1];
rz(-2.6458461) q[1];
sx q[1];
rz(-0.31235487) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6808514) q[0];
sx q[0];
rz(-0.4609403) q[0];
sx q[0];
rz(-2.8747706) q[0];
x q[1];
rz(1.4710452) q[2];
sx q[2];
rz(-1.3653697) q[2];
sx q[2];
rz(-0.86459407) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4620381) q[1];
sx q[1];
rz(-1.7741388) q[1];
sx q[1];
rz(0.36851369) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7658224) q[3];
sx q[3];
rz(-0.93149501) q[3];
sx q[3];
rz(2.2488058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.13545869) q[2];
sx q[2];
rz(-1.6929071) q[2];
sx q[2];
rz(2.3166166) q[2];
rz(-1.4247591) q[3];
sx q[3];
rz(-0.32130876) q[3];
sx q[3];
rz(-2.7062374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7216126) q[0];
sx q[0];
rz(-1.551832) q[0];
sx q[0];
rz(2.2671674) q[0];
rz(0.32886109) q[1];
sx q[1];
rz(-1.7560274) q[1];
sx q[1];
rz(1.0985451) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4916353) q[0];
sx q[0];
rz(-1.1960317) q[0];
sx q[0];
rz(1.6894421) q[0];
rz(-0.036089049) q[2];
sx q[2];
rz(-0.83613013) q[2];
sx q[2];
rz(0.058908894) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.85531536) q[1];
sx q[1];
rz(-2.554727) q[1];
sx q[1];
rz(-1.0081968) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1608814) q[3];
sx q[3];
rz(-1.3385337) q[3];
sx q[3];
rz(1.4651608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.99981368) q[2];
sx q[2];
rz(-1.9580611) q[2];
sx q[2];
rz(2.4675274) q[2];
rz(1.4874124) q[3];
sx q[3];
rz(-1.4451278) q[3];
sx q[3];
rz(-2.8580247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0312626) q[0];
sx q[0];
rz(-2.1188348) q[0];
sx q[0];
rz(-1.7559825) q[0];
rz(-1.1194057) q[1];
sx q[1];
rz(-1.4675354) q[1];
sx q[1];
rz(-1.3853692) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0273756) q[0];
sx q[0];
rz(-1.3189335) q[0];
sx q[0];
rz(0.92433874) q[0];
rz(-pi) q[1];
rz(-2.9846956) q[2];
sx q[2];
rz(-0.50373915) q[2];
sx q[2];
rz(0.7157601) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0138356) q[1];
sx q[1];
rz(-2.7633939) q[1];
sx q[1];
rz(2.5950123) q[1];
rz(-2.9117111) q[3];
sx q[3];
rz(-2.982989) q[3];
sx q[3];
rz(1.6386709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0783656) q[2];
sx q[2];
rz(-2.207022) q[2];
sx q[2];
rz(2.7654977) q[2];
rz(-2.0070576) q[3];
sx q[3];
rz(-0.36342707) q[3];
sx q[3];
rz(-2.334972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1217839) q[0];
sx q[0];
rz(-0.60096318) q[0];
sx q[0];
rz(3.0390749) q[0];
rz(-2.5566697) q[1];
sx q[1];
rz(-0.94766098) q[1];
sx q[1];
rz(-1.1835416) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69307704) q[0];
sx q[0];
rz(-1.428837) q[0];
sx q[0];
rz(-1.7008812) q[0];
rz(-pi) q[1];
rz(-1.905422) q[2];
sx q[2];
rz(-2.2058704) q[2];
sx q[2];
rz(-2.2872567) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.592748) q[1];
sx q[1];
rz(-1.3164489) q[1];
sx q[1];
rz(-0.91491048) q[1];
rz(-pi) q[2];
rz(0.4966708) q[3];
sx q[3];
rz(-2.6210945) q[3];
sx q[3];
rz(-0.48540533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79503235) q[2];
sx q[2];
rz(-0.9298032) q[2];
sx q[2];
rz(-2.9417876) q[2];
rz(-2.0810769) q[3];
sx q[3];
rz(-0.54148713) q[3];
sx q[3];
rz(0.04960355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.878433) q[0];
sx q[0];
rz(-1.4819772) q[0];
sx q[0];
rz(-0.31295452) q[0];
rz(2.1921659) q[1];
sx q[1];
rz(-1.3525617) q[1];
sx q[1];
rz(-1.4535646) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2672103) q[0];
sx q[0];
rz(-1.4802762) q[0];
sx q[0];
rz(-0.1874013) q[0];
x q[1];
rz(-1.6498927) q[2];
sx q[2];
rz(-2.1888615) q[2];
sx q[2];
rz(-1.6343509) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1968699) q[1];
sx q[1];
rz(-2.2579282) q[1];
sx q[1];
rz(-1.8623167) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9096776) q[3];
sx q[3];
rz(-2.0032855) q[3];
sx q[3];
rz(2.5755455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.35385418) q[2];
sx q[2];
rz(-2.165803) q[2];
sx q[2];
rz(-1.3346416) q[2];
rz(1.8169962) q[3];
sx q[3];
rz(-2.4748804) q[3];
sx q[3];
rz(1.9273531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.036309328) q[0];
sx q[0];
rz(-2.5258625) q[0];
sx q[0];
rz(0.69806725) q[0];
rz(2.2659194) q[1];
sx q[1];
rz(-2.0798637) q[1];
sx q[1];
rz(2.7811513) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1669641) q[0];
sx q[0];
rz(-0.95050838) q[0];
sx q[0];
rz(-1.0247158) q[0];
x q[1];
rz(-1.5064042) q[2];
sx q[2];
rz(-1.1887822) q[2];
sx q[2];
rz(2.7329993) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0086509) q[1];
sx q[1];
rz(-1.1906149) q[1];
sx q[1];
rz(-0.93312414) q[1];
rz(1.4029986) q[3];
sx q[3];
rz(-1.0281963) q[3];
sx q[3];
rz(0.29796539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2299819) q[2];
sx q[2];
rz(-1.4260099) q[2];
sx q[2];
rz(-0.81725517) q[2];
rz(-2.3115555) q[3];
sx q[3];
rz(-0.54098141) q[3];
sx q[3];
rz(2.9221007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5068186) q[0];
sx q[0];
rz(-2.2932678) q[0];
sx q[0];
rz(3.1095374) q[0];
rz(1.3196779) q[1];
sx q[1];
rz(-1.3751043) q[1];
sx q[1];
rz(2.4874036) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.74943) q[0];
sx q[0];
rz(-1.6902958) q[0];
sx q[0];
rz(0.197535) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0243353) q[2];
sx q[2];
rz(-0.83570601) q[2];
sx q[2];
rz(2.105956) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3820104) q[1];
sx q[1];
rz(-2.0339662) q[1];
sx q[1];
rz(2.5968371) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1919674) q[3];
sx q[3];
rz(-2.4740268) q[3];
sx q[3];
rz(2.6258056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0111982) q[2];
sx q[2];
rz(-2.0989213) q[2];
sx q[2];
rz(-1.0065669) q[2];
rz(-2.448163) q[3];
sx q[3];
rz(-1.3207685) q[3];
sx q[3];
rz(-1.8600195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7158159) q[0];
sx q[0];
rz(-2.4632813) q[0];
sx q[0];
rz(3.0671469) q[0];
rz(2.2609932) q[1];
sx q[1];
rz(-1.1136628) q[1];
sx q[1];
rz(-2.640092) q[1];
rz(0.975561) q[2];
sx q[2];
rz(-0.79339334) q[2];
sx q[2];
rz(-1.7421772) q[2];
rz(1.1893336) q[3];
sx q[3];
rz(-0.91100024) q[3];
sx q[3];
rz(-0.80445214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
