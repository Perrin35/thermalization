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
rz(1.2020741) q[0];
sx q[0];
rz(-2.095686) q[0];
sx q[0];
rz(2.3910971) q[0];
rz(-2.9208288) q[1];
sx q[1];
rz(-1.060744) q[1];
sx q[1];
rz(-1.9579252) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1738244) q[0];
sx q[0];
rz(-0.90961752) q[0];
sx q[0];
rz(-2.1029766) q[0];
rz(1.4040825) q[2];
sx q[2];
rz(-1.445747) q[2];
sx q[2];
rz(0.84593433) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4422801) q[1];
sx q[1];
rz(-0.89095014) q[1];
sx q[1];
rz(-1.555843) q[1];
x q[2];
rz(1.2624083) q[3];
sx q[3];
rz(-1.3779791) q[3];
sx q[3];
rz(-2.8578025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3262647) q[2];
sx q[2];
rz(-1.1417737) q[2];
sx q[2];
rz(0.092078837) q[2];
rz(2.0103256) q[3];
sx q[3];
rz(-1.0695894) q[3];
sx q[3];
rz(-2.2856975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6446514) q[0];
sx q[0];
rz(-2.6583789) q[0];
sx q[0];
rz(-0.53571969) q[0];
rz(-0.016853111) q[1];
sx q[1];
rz(-0.87937513) q[1];
sx q[1];
rz(-3.1244315) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6494228) q[0];
sx q[0];
rz(-1.0809521) q[0];
sx q[0];
rz(0.36618284) q[0];
x q[1];
rz(0.13249915) q[2];
sx q[2];
rz(-2.6713058) q[2];
sx q[2];
rz(-0.4343701) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.90777724) q[1];
sx q[1];
rz(-1.3388947) q[1];
sx q[1];
rz(1.7323937) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0073418) q[3];
sx q[3];
rz(-1.6608616) q[3];
sx q[3];
rz(1.4081189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0422334) q[2];
sx q[2];
rz(-1.6681654) q[2];
sx q[2];
rz(-0.24822203) q[2];
rz(0.92784268) q[3];
sx q[3];
rz(-0.78301269) q[3];
sx q[3];
rz(-0.46970126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27702734) q[0];
sx q[0];
rz(-1.1772573) q[0];
sx q[0];
rz(-2.9789341) q[0];
rz(-2.3795369) q[1];
sx q[1];
rz(-0.22480741) q[1];
sx q[1];
rz(0.83388296) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8441761) q[0];
sx q[0];
rz(-2.5582097) q[0];
sx q[0];
rz(-2.9177505) q[0];
rz(-0.35807407) q[2];
sx q[2];
rz(-0.34293567) q[2];
sx q[2];
rz(-1.6274522) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0790339) q[1];
sx q[1];
rz(-1.9818881) q[1];
sx q[1];
rz(2.316019) q[1];
rz(0.44077222) q[3];
sx q[3];
rz(-1.5197872) q[3];
sx q[3];
rz(0.15404242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0791846) q[2];
sx q[2];
rz(-2.0024313) q[2];
sx q[2];
rz(2.8774234) q[2];
rz(-2.0832113) q[3];
sx q[3];
rz(-0.94401413) q[3];
sx q[3];
rz(-3.1112166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77085483) q[0];
sx q[0];
rz(-0.83468947) q[0];
sx q[0];
rz(-1.7133065) q[0];
rz(2.3772073) q[1];
sx q[1];
rz(-1.1420219) q[1];
sx q[1];
rz(1.3216602) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3300312) q[0];
sx q[0];
rz(-1.9178277) q[0];
sx q[0];
rz(2.7156732) q[0];
x q[1];
rz(3.0899449) q[2];
sx q[2];
rz(-1.2394636) q[2];
sx q[2];
rz(0.29690642) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9191967) q[1];
sx q[1];
rz(-0.42342227) q[1];
sx q[1];
rz(-2.3906796) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0760078) q[3];
sx q[3];
rz(-2.3632999) q[3];
sx q[3];
rz(-0.28631223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4798212) q[2];
sx q[2];
rz(-1.4134553) q[2];
sx q[2];
rz(-2.856355) q[2];
rz(-1.8789004) q[3];
sx q[3];
rz(-1.217239) q[3];
sx q[3];
rz(-1.4234022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30094639) q[0];
sx q[0];
rz(-0.74746376) q[0];
sx q[0];
rz(-0.88229156) q[0];
rz(0.67059416) q[1];
sx q[1];
rz(-2.0457485) q[1];
sx q[1];
rz(-2.1569599) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0142954) q[0];
sx q[0];
rz(-1.6335575) q[0];
sx q[0];
rz(1.5393637) q[0];
x q[1];
rz(-3.0987344) q[2];
sx q[2];
rz(-2.1120666) q[2];
sx q[2];
rz(1.0554016) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.24422503) q[1];
sx q[1];
rz(-2.3363718) q[1];
sx q[1];
rz(2.3697788) q[1];
rz(-pi) q[2];
x q[2];
rz(2.601131) q[3];
sx q[3];
rz(-1.988388) q[3];
sx q[3];
rz(2.3706974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.275445) q[2];
sx q[2];
rz(-0.53469849) q[2];
sx q[2];
rz(-0.61694413) q[2];
rz(-2.3894737) q[3];
sx q[3];
rz(-1.813846) q[3];
sx q[3];
rz(-0.45929685) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98944703) q[0];
sx q[0];
rz(-2.4247657) q[0];
sx q[0];
rz(-0.77051198) q[0];
rz(2.41467) q[1];
sx q[1];
rz(-2.9194071) q[1];
sx q[1];
rz(0.70473421) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7677697) q[0];
sx q[0];
rz(-0.85568586) q[0];
sx q[0];
rz(0.82775791) q[0];
rz(-pi) q[1];
rz(2.0915085) q[2];
sx q[2];
rz(-1.8936704) q[2];
sx q[2];
rz(2.8276557) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.875346) q[1];
sx q[1];
rz(-2.7901712) q[1];
sx q[1];
rz(-3.0446192) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1416311) q[3];
sx q[3];
rz(-1.3558767) q[3];
sx q[3];
rz(0.31610148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3720234) q[2];
sx q[2];
rz(-1.4906887) q[2];
sx q[2];
rz(0.43323576) q[2];
rz(0.20089928) q[3];
sx q[3];
rz(-1.9379987) q[3];
sx q[3];
rz(0.5350298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(2.3055426) q[0];
sx q[0];
rz(-1.1490281) q[0];
sx q[0];
rz(-2.7591163) q[0];
rz(1.0410694) q[1];
sx q[1];
rz(-1.46547) q[1];
sx q[1];
rz(-0.44630757) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.324821) q[0];
sx q[0];
rz(-1.3123625) q[0];
sx q[0];
rz(-1.517061) q[0];
rz(-1.3417753) q[2];
sx q[2];
rz(-1.9441368) q[2];
sx q[2];
rz(2.7315706) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.016231) q[1];
sx q[1];
rz(-0.70763677) q[1];
sx q[1];
rz(1.0098004) q[1];
rz(-0.84844037) q[3];
sx q[3];
rz(-2.1506243) q[3];
sx q[3];
rz(2.8715796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.81991601) q[2];
sx q[2];
rz(-1.1175464) q[2];
sx q[2];
rz(-2.4134911) q[2];
rz(0.32650945) q[3];
sx q[3];
rz(-1.510334) q[3];
sx q[3];
rz(-2.2842893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2364872) q[0];
sx q[0];
rz(-0.40920722) q[0];
sx q[0];
rz(-1.8439199) q[0];
rz(2.1406651) q[1];
sx q[1];
rz(-0.74556723) q[1];
sx q[1];
rz(0.40850684) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4057874) q[0];
sx q[0];
rz(-1.5716432) q[0];
sx q[0];
rz(-1.5849831) q[0];
rz(-pi) q[1];
rz(-0.82168545) q[2];
sx q[2];
rz(-0.61640451) q[2];
sx q[2];
rz(-2.6954755) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0169819) q[1];
sx q[1];
rz(-1.4370185) q[1];
sx q[1];
rz(-3.0264808) q[1];
rz(-0.85902416) q[3];
sx q[3];
rz(-0.54104303) q[3];
sx q[3];
rz(-2.4699901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0268176) q[2];
sx q[2];
rz(-1.6208181) q[2];
sx q[2];
rz(0.82359037) q[2];
rz(2.1987727) q[3];
sx q[3];
rz(-1.7190245) q[3];
sx q[3];
rz(-1.6387117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59458643) q[0];
sx q[0];
rz(-2.6972045) q[0];
sx q[0];
rz(-1.2521) q[0];
rz(-1.3638672) q[1];
sx q[1];
rz(-1.6412647) q[1];
sx q[1];
rz(1.3346671) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20864381) q[0];
sx q[0];
rz(-1.874039) q[0];
sx q[0];
rz(-1.6088845) q[0];
rz(-pi) q[1];
rz(-0.64200114) q[2];
sx q[2];
rz(-1.4307012) q[2];
sx q[2];
rz(2.7336655) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6872944) q[1];
sx q[1];
rz(-2.0210938) q[1];
sx q[1];
rz(1.4895093) q[1];
rz(-pi) q[2];
rz(1.077721) q[3];
sx q[3];
rz(-2.6209313) q[3];
sx q[3];
rz(-0.23948224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1198662) q[2];
sx q[2];
rz(-2.848337) q[2];
sx q[2];
rz(0.29652706) q[2];
rz(-1.6622539) q[3];
sx q[3];
rz(-1.5115503) q[3];
sx q[3];
rz(-2.9843946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6327561) q[0];
sx q[0];
rz(-2.4907676) q[0];
sx q[0];
rz(0.79488361) q[0];
rz(2.5859313) q[1];
sx q[1];
rz(-1.2763005) q[1];
sx q[1];
rz(1.6937675) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2771196) q[0];
sx q[0];
rz(-2.2896447) q[0];
sx q[0];
rz(2.6299631) q[0];
rz(-pi) q[1];
rz(-0.49782355) q[2];
sx q[2];
rz(-1.6256672) q[2];
sx q[2];
rz(0.16615088) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.70778479) q[1];
sx q[1];
rz(-2.6205739) q[1];
sx q[1];
rz(1.1624813) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1673739) q[3];
sx q[3];
rz(-2.0250626) q[3];
sx q[3];
rz(-2.5429986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29622233) q[2];
sx q[2];
rz(-1.6310548) q[2];
sx q[2];
rz(-1.4307865) q[2];
rz(2.3885662) q[3];
sx q[3];
rz(-2.3236661) q[3];
sx q[3];
rz(-0.16451612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0059148) q[0];
sx q[0];
rz(-2.4527241) q[0];
sx q[0];
rz(1.2149568) q[0];
rz(-1.6750492) q[1];
sx q[1];
rz(-0.92862447) q[1];
sx q[1];
rz(-2.5359572) q[1];
rz(-2.5765924) q[2];
sx q[2];
rz(-1.0641268) q[2];
sx q[2];
rz(0.63882154) q[2];
rz(1.4004272) q[3];
sx q[3];
rz(-1.6765249) q[3];
sx q[3];
rz(-0.75679814) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
