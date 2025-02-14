OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72579223) q[0];
sx q[0];
rz(-1.1644726) q[0];
sx q[0];
rz(1.8902984) q[0];
rz(2.8407821) q[1];
sx q[1];
rz(-1.7164813) q[1];
sx q[1];
rz(-3.0256174) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.011256) q[0];
sx q[0];
rz(-2.7948871) q[0];
sx q[0];
rz(-1.1954855) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1484199) q[2];
sx q[2];
rz(-2.6184888) q[2];
sx q[2];
rz(2.1577699) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8311316) q[1];
sx q[1];
rz(-2.3678105) q[1];
sx q[1];
rz(1.6602064) q[1];
rz(0.54268022) q[3];
sx q[3];
rz(-1.451406) q[3];
sx q[3];
rz(-0.38633864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3731602) q[2];
sx q[2];
rz(-0.065040437) q[2];
sx q[2];
rz(-1.6981286) q[2];
rz(0.69624919) q[3];
sx q[3];
rz(-1.6257476) q[3];
sx q[3];
rz(-0.35370091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26831216) q[0];
sx q[0];
rz(-3.066635) q[0];
sx q[0];
rz(-2.7798376) q[0];
rz(-1.7573645) q[1];
sx q[1];
rz(-0.39159602) q[1];
sx q[1];
rz(1.0272383) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7536094) q[0];
sx q[0];
rz(-1.7513203) q[0];
sx q[0];
rz(0.87708824) q[0];
x q[1];
rz(2.9774819) q[2];
sx q[2];
rz(-2.3375953) q[2];
sx q[2];
rz(1.6156593) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9578182) q[1];
sx q[1];
rz(-1.8190067) q[1];
sx q[1];
rz(0.91285984) q[1];
rz(-pi) q[2];
rz(-1.6546685) q[3];
sx q[3];
rz(-1.7954651) q[3];
sx q[3];
rz(-0.39473907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1343753) q[2];
sx q[2];
rz(-0.5642887) q[2];
sx q[2];
rz(-1.0312274) q[2];
rz(2.4550896) q[3];
sx q[3];
rz(-1.406823) q[3];
sx q[3];
rz(-2.8197207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4903851) q[0];
sx q[0];
rz(-1.729916) q[0];
sx q[0];
rz(1.7899293) q[0];
rz(1.5216113) q[1];
sx q[1];
rz(-0.23623513) q[1];
sx q[1];
rz(-1.1865541) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31165115) q[0];
sx q[0];
rz(-1.9970987) q[0];
sx q[0];
rz(-1.9894129) q[0];
rz(-pi) q[1];
rz(1.4954689) q[2];
sx q[2];
rz(-1.1161303) q[2];
sx q[2];
rz(2.1476434) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.12315519) q[1];
sx q[1];
rz(-1.5644524) q[1];
sx q[1];
rz(-3.1325587) q[1];
rz(2.260798) q[3];
sx q[3];
rz(-1.8024416) q[3];
sx q[3];
rz(0.85154545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.97570193) q[2];
sx q[2];
rz(-1.8791135) q[2];
sx q[2];
rz(-0.28124896) q[2];
rz(0.37505546) q[3];
sx q[3];
rz(-2.227759) q[3];
sx q[3];
rz(-0.32947549) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81958812) q[0];
sx q[0];
rz(-2.3276038) q[0];
sx q[0];
rz(2.8170012) q[0];
rz(1.5931256) q[1];
sx q[1];
rz(-0.32061583) q[1];
sx q[1];
rz(-2.2494907) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4504334) q[0];
sx q[0];
rz(-1.5884134) q[0];
sx q[0];
rz(-1.5595647) q[0];
x q[1];
rz(0.11406819) q[2];
sx q[2];
rz(-2.4270396) q[2];
sx q[2];
rz(-2.8482506) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0405691) q[1];
sx q[1];
rz(-1.8547579) q[1];
sx q[1];
rz(-0.26974399) q[1];
x q[2];
rz(3.1079477) q[3];
sx q[3];
rz(-0.97760289) q[3];
sx q[3];
rz(0.9481687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0804245) q[2];
sx q[2];
rz(-1.3865546) q[2];
sx q[2];
rz(-2.5648153) q[2];
rz(-3.0171677) q[3];
sx q[3];
rz(-0.7500698) q[3];
sx q[3];
rz(-0.09593825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4633923) q[0];
sx q[0];
rz(-0.41825774) q[0];
sx q[0];
rz(2.8611355) q[0];
rz(-2.8434143) q[1];
sx q[1];
rz(-3.0715946) q[1];
sx q[1];
rz(-0.15204522) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9579118) q[0];
sx q[0];
rz(-3.0820345) q[0];
sx q[0];
rz(-2.0516711) q[0];
rz(-pi) q[1];
rz(1.399081) q[2];
sx q[2];
rz(-2.6125312) q[2];
sx q[2];
rz(2.3429639) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.070134461) q[1];
sx q[1];
rz(-1.9862174) q[1];
sx q[1];
rz(-0.43674119) q[1];
rz(1.9535311) q[3];
sx q[3];
rz(-1.5627341) q[3];
sx q[3];
rz(0.950056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8449479) q[2];
sx q[2];
rz(-1.351202) q[2];
sx q[2];
rz(-0.039999261) q[2];
rz(0.1376888) q[3];
sx q[3];
rz(-0.44860336) q[3];
sx q[3];
rz(2.3562446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8765151) q[0];
sx q[0];
rz(-0.10049937) q[0];
sx q[0];
rz(-2.5608089) q[0];
rz(0.68897092) q[1];
sx q[1];
rz(-0.13801485) q[1];
sx q[1];
rz(1.2977915) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9823124) q[0];
sx q[0];
rz(-1.1508688) q[0];
sx q[0];
rz(-0.4585271) q[0];
rz(-pi) q[1];
x q[1];
rz(3.077438) q[2];
sx q[2];
rz(-1.0279845) q[2];
sx q[2];
rz(-2.9911957) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4124085) q[1];
sx q[1];
rz(-1.3505591) q[1];
sx q[1];
rz(-2.7413648) q[1];
x q[2];
rz(-2.5291555) q[3];
sx q[3];
rz(-2.4531657) q[3];
sx q[3];
rz(1.5189309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40655228) q[2];
sx q[2];
rz(-1.8031969) q[2];
sx q[2];
rz(-1.5470541) q[2];
rz(2.7100587) q[3];
sx q[3];
rz(-2.0412692) q[3];
sx q[3];
rz(-2.4908454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65289098) q[0];
sx q[0];
rz(-0.015691375) q[0];
sx q[0];
rz(2.5456862) q[0];
rz(-1.3504922) q[1];
sx q[1];
rz(-0.48521438) q[1];
sx q[1];
rz(-0.57292604) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6011621) q[0];
sx q[0];
rz(-1.5613544) q[0];
sx q[0];
rz(-1.5736395) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7155034) q[2];
sx q[2];
rz(-1.3182148) q[2];
sx q[2];
rz(-1.9468228) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0966037) q[1];
sx q[1];
rz(-1.8879688) q[1];
sx q[1];
rz(-0.60896341) q[1];
x q[2];
rz(2.7468729) q[3];
sx q[3];
rz(-0.79666797) q[3];
sx q[3];
rz(-1.7834237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4493745) q[2];
sx q[2];
rz(-2.2276679) q[2];
sx q[2];
rz(-0.91478983) q[2];
rz(-1.5335013) q[3];
sx q[3];
rz(-1.1398818) q[3];
sx q[3];
rz(-2.1515414) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6221301) q[0];
sx q[0];
rz(-1.7161481) q[0];
sx q[0];
rz(-0.064432681) q[0];
rz(0.26647767) q[1];
sx q[1];
rz(-3.0173512) q[1];
sx q[1];
rz(-1.2640094) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5504996) q[0];
sx q[0];
rz(-1.8944709) q[0];
sx q[0];
rz(-1.8514015) q[0];
rz(-pi) q[1];
rz(2.6904628) q[2];
sx q[2];
rz(-2.3887803) q[2];
sx q[2];
rz(0.058016404) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.83710421) q[1];
sx q[1];
rz(-1.7715369) q[1];
sx q[1];
rz(0.31977579) q[1];
x q[2];
rz(-2.1355992) q[3];
sx q[3];
rz(-2.0801615) q[3];
sx q[3];
rz(2.6107174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6182575) q[2];
sx q[2];
rz(-1.9674415) q[2];
sx q[2];
rz(2.0972283) q[2];
rz(0.70762819) q[3];
sx q[3];
rz(-0.39053598) q[3];
sx q[3];
rz(2.0948758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7758961) q[0];
sx q[0];
rz(-1.5713659) q[0];
sx q[0];
rz(-0.43029341) q[0];
rz(2.4039092) q[1];
sx q[1];
rz(-0.084641181) q[1];
sx q[1];
rz(-2.784909) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1708313) q[0];
sx q[0];
rz(-2.8342842) q[0];
sx q[0];
rz(-0.31913646) q[0];
rz(-pi) q[1];
rz(-1.6820774) q[2];
sx q[2];
rz(-0.26646895) q[2];
sx q[2];
rz(-0.9916208) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.99706471) q[1];
sx q[1];
rz(-0.72880077) q[1];
sx q[1];
rz(-0.65741922) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55583007) q[3];
sx q[3];
rz(-0.14066108) q[3];
sx q[3];
rz(2.6859796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3714527) q[2];
sx q[2];
rz(-0.74718237) q[2];
sx q[2];
rz(2.746197) q[2];
rz(-1.4622408) q[3];
sx q[3];
rz(-1.3926927) q[3];
sx q[3];
rz(-3.1126378) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28703654) q[0];
sx q[0];
rz(-0.77571464) q[0];
sx q[0];
rz(-2.4391644) q[0];
rz(-0.17829819) q[1];
sx q[1];
rz(-2.4874004) q[1];
sx q[1];
rz(-1.5231232) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1078774) q[0];
sx q[0];
rz(-1.807567) q[0];
sx q[0];
rz(-1.9997755) q[0];
rz(1.4728925) q[2];
sx q[2];
rz(-1.5852829) q[2];
sx q[2];
rz(1.9928428) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2158015) q[1];
sx q[1];
rz(-0.58432441) q[1];
sx q[1];
rz(-0.25694848) q[1];
x q[2];
rz(-0.50242063) q[3];
sx q[3];
rz(-1.9137041) q[3];
sx q[3];
rz(-0.08758647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6452597) q[2];
sx q[2];
rz(-2.1767904) q[2];
sx q[2];
rz(-2.1920835) q[2];
rz(2.005714) q[3];
sx q[3];
rz(-0.94616008) q[3];
sx q[3];
rz(0.56376636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5180494) q[0];
sx q[0];
rz(-1.8907056) q[0];
sx q[0];
rz(-0.84778669) q[0];
rz(-1.7060948) q[1];
sx q[1];
rz(-1.2789627) q[1];
sx q[1];
rz(0.36437558) q[1];
rz(-3.0355787) q[2];
sx q[2];
rz(-2.2681375) q[2];
sx q[2];
rz(-2.9974147) q[2];
rz(-1.3798716) q[3];
sx q[3];
rz(-2.8287201) q[3];
sx q[3];
rz(-0.74486275) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
