OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(2.2979484) q[0];
sx q[0];
rz(9.2568682) q[0];
rz(-1.9703938) q[1];
sx q[1];
rz(-0.29532239) q[1];
sx q[1];
rz(3.0854316) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18626285) q[0];
sx q[0];
rz(-1.3583399) q[0];
sx q[0];
rz(2.0583378) q[0];
rz(2.216823) q[2];
sx q[2];
rz(-1.3999108) q[2];
sx q[2];
rz(2.8303763) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.15966283) q[1];
sx q[1];
rz(-0.61239457) q[1];
sx q[1];
rz(-2.0484522) q[1];
rz(-0.25986259) q[3];
sx q[3];
rz(-1.5279603) q[3];
sx q[3];
rz(-0.44997893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7636259) q[2];
sx q[2];
rz(-2.8597735) q[2];
sx q[2];
rz(-0.4326694) q[2];
rz(-1.9487322) q[3];
sx q[3];
rz(-1.2377219) q[3];
sx q[3];
rz(0.38309923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52779657) q[0];
sx q[0];
rz(-0.48848099) q[0];
sx q[0];
rz(-1.312785) q[0];
rz(0.20547543) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(1.1516494) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5702471) q[0];
sx q[0];
rz(-2.3803108) q[0];
sx q[0];
rz(1.135457) q[0];
rz(-pi) q[1];
rz(1.0863016) q[2];
sx q[2];
rz(-1.0420348) q[2];
sx q[2];
rz(-0.27258401) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0484867) q[1];
sx q[1];
rz(-1.9890607) q[1];
sx q[1];
rz(0.47936819) q[1];
x q[2];
rz(-2.6395256) q[3];
sx q[3];
rz(-1.2289398) q[3];
sx q[3];
rz(-1.348192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1318704) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(-2.9197664) q[2];
rz(-0.37718537) q[3];
sx q[3];
rz(-0.42755869) q[3];
sx q[3];
rz(2.3691573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31056988) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(0.036852766) q[0];
rz(-2.316078) q[1];
sx q[1];
rz(-1.8258391) q[1];
sx q[1];
rz(0.056578606) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.409257) q[0];
sx q[0];
rz(-2.0154672) q[0];
sx q[0];
rz(-2.1952941) q[0];
rz(-1.4372196) q[2];
sx q[2];
rz(-1.3165511) q[2];
sx q[2];
rz(-0.1453407) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.024836) q[1];
sx q[1];
rz(-2.7783238) q[1];
sx q[1];
rz(2.5519752) q[1];
rz(-pi) q[2];
rz(1.9583086) q[3];
sx q[3];
rz(-0.96252493) q[3];
sx q[3];
rz(2.0351978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8686707) q[2];
sx q[2];
rz(-1.4977095) q[2];
sx q[2];
rz(2.2154714) q[2];
rz(-2.5849294) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(-1.042897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91519231) q[0];
sx q[0];
rz(-0.72015786) q[0];
sx q[0];
rz(-0.74209374) q[0];
rz(1.1391976) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(2.6779968) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47953654) q[0];
sx q[0];
rz(-2.2502796) q[0];
sx q[0];
rz(-0.38328538) q[0];
rz(0.17642994) q[2];
sx q[2];
rz(-0.33761218) q[2];
sx q[2];
rz(-0.38844973) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3568748) q[1];
sx q[1];
rz(-2.3349635) q[1];
sx q[1];
rz(-0.23826092) q[1];
rz(-2.1257524) q[3];
sx q[3];
rz(-1.1557126) q[3];
sx q[3];
rz(-0.89158981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3670369) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(-3.0920933) q[2];
rz(0.1285304) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(-0.11894225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13609919) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(-2.8438925) q[0];
rz(2.659335) q[1];
sx q[1];
rz(-2.3869956) q[1];
sx q[1];
rz(0.94435) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0194861) q[0];
sx q[0];
rz(-1.6158316) q[0];
sx q[0];
rz(-1.7800063) q[0];
rz(0.40446754) q[2];
sx q[2];
rz(-0.81768113) q[2];
sx q[2];
rz(0.58194619) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7625092) q[1];
sx q[1];
rz(-0.06113872) q[1];
sx q[1];
rz(1.6054543) q[1];
rz(-pi) q[2];
rz(-0.51575757) q[3];
sx q[3];
rz(-2.6532647) q[3];
sx q[3];
rz(0.26667903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.258761) q[2];
sx q[2];
rz(-1.0927039) q[2];
sx q[2];
rz(-3.0333701) q[2];
rz(3.1392858) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7047983) q[0];
sx q[0];
rz(-2.695485) q[0];
sx q[0];
rz(0.58445245) q[0];
rz(0.8862409) q[1];
sx q[1];
rz(-0.61683547) q[1];
sx q[1];
rz(3.086673) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4024825) q[0];
sx q[0];
rz(-1.4882898) q[0];
sx q[0];
rz(-1.8399747) q[0];
rz(-pi) q[1];
rz(-0.018718406) q[2];
sx q[2];
rz(-1.8939549) q[2];
sx q[2];
rz(-1.2013555) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1153032) q[1];
sx q[1];
rz(-2.4541828) q[1];
sx q[1];
rz(-0.61647146) q[1];
x q[2];
rz(-0.46338007) q[3];
sx q[3];
rz(-1.0372835) q[3];
sx q[3];
rz(0.86138553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3871258) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(-1.0167271) q[2];
rz(-0.54404849) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(1.3325161) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.465437) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(0.28453919) q[0];
rz(-2.1971205) q[1];
sx q[1];
rz(-1.1962793) q[1];
sx q[1];
rz(-0.91032666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0383354) q[0];
sx q[0];
rz(-3.0568125) q[0];
sx q[0];
rz(2.2708562) q[0];
rz(2.2482713) q[2];
sx q[2];
rz(-2.6258694) q[2];
sx q[2];
rz(-0.90781462) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6701723) q[1];
sx q[1];
rz(-0.52392611) q[1];
sx q[1];
rz(0.36797221) q[1];
x q[2];
rz(-2.4207553) q[3];
sx q[3];
rz(-1.1158873) q[3];
sx q[3];
rz(0.37973675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7515144) q[2];
sx q[2];
rz(-0.063515924) q[2];
sx q[2];
rz(0.92203036) q[2];
rz(0.56728029) q[3];
sx q[3];
rz(-1.7276948) q[3];
sx q[3];
rz(-2.1218307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89408016) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(0.12938736) q[0];
rz(-2.5091876) q[1];
sx q[1];
rz(-1.0267195) q[1];
sx q[1];
rz(-0.30050373) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8757272) q[0];
sx q[0];
rz(-2.2458796) q[0];
sx q[0];
rz(2.6595594) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95327611) q[2];
sx q[2];
rz(-0.91070181) q[2];
sx q[2];
rz(-2.2613139) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6943372) q[1];
sx q[1];
rz(-1.2303423) q[1];
sx q[1];
rz(-1.3854331) q[1];
rz(-pi) q[2];
rz(-2.825533) q[3];
sx q[3];
rz(-0.83871597) q[3];
sx q[3];
rz(2.4162606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.58632103) q[2];
sx q[2];
rz(-0.99493146) q[2];
sx q[2];
rz(-2.3596181) q[2];
rz(2.590495) q[3];
sx q[3];
rz(-1.7588153) q[3];
sx q[3];
rz(0.15792318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5732116) q[0];
sx q[0];
rz(-1.9848354) q[0];
sx q[0];
rz(0.12776275) q[0];
rz(-0.54221517) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(-2.382747) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1369143) q[0];
sx q[0];
rz(-1.5877962) q[0];
sx q[0];
rz(-1.1466115) q[0];
rz(-pi) q[1];
rz(1.3779638) q[2];
sx q[2];
rz(-0.97059965) q[2];
sx q[2];
rz(-0.45229518) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54770494) q[1];
sx q[1];
rz(-0.87915671) q[1];
sx q[1];
rz(-1.3681075) q[1];
x q[2];
rz(-1.5893448) q[3];
sx q[3];
rz(-2.4759001) q[3];
sx q[3];
rz(1.9513643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1252497) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(-0.49003595) q[2];
rz(-1.4222493) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(1.0629883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35995099) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(3.066257) q[0];
rz(0.8967337) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(2.5316701) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41994914) q[0];
sx q[0];
rz(-2.3099265) q[0];
sx q[0];
rz(1.9589817) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5209849) q[2];
sx q[2];
rz(-0.45244103) q[2];
sx q[2];
rz(-1.1021745) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6590609) q[1];
sx q[1];
rz(-0.7809124) q[1];
sx q[1];
rz(-2.798583) q[1];
x q[2];
rz(-2.6033953) q[3];
sx q[3];
rz(-2.3231069) q[3];
sx q[3];
rz(-1.324211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.909409) q[2];
sx q[2];
rz(-2.3164618) q[2];
sx q[2];
rz(0.71371901) q[2];
rz(2.7632726) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(2.2617214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80355766) q[0];
sx q[0];
rz(-1.9914347) q[0];
sx q[0];
rz(1.5557355) q[0];
rz(0.65080416) q[1];
sx q[1];
rz(-1.6497859) q[1];
sx q[1];
rz(-0.12129687) q[1];
rz(-1.7238293) q[2];
sx q[2];
rz(-2.9433708) q[2];
sx q[2];
rz(2.3797258) q[2];
rz(-0.25272947) q[3];
sx q[3];
rz(-2.0287632) q[3];
sx q[3];
rz(-0.91026929) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];