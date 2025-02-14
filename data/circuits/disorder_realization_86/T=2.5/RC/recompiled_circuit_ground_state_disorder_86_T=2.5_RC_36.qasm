OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0107083) q[0];
sx q[0];
rz(-2.4680128) q[0];
sx q[0];
rz(-2.3186865) q[0];
rz(1.4505439) q[1];
sx q[1];
rz(-0.6757285) q[1];
sx q[1];
rz(0.20767173) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0820706) q[0];
sx q[0];
rz(-2.2769109) q[0];
sx q[0];
rz(0.45990277) q[0];
x q[1];
rz(2.093621) q[2];
sx q[2];
rz(-0.62467161) q[2];
sx q[2];
rz(-2.5324948) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.365578) q[1];
sx q[1];
rz(-1.509359) q[1];
sx q[1];
rz(3.0173124) q[1];
rz(-pi) q[2];
x q[2];
rz(1.200233) q[3];
sx q[3];
rz(-1.6374805) q[3];
sx q[3];
rz(-0.62180623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9660008) q[2];
sx q[2];
rz(-1.5755743) q[2];
sx q[2];
rz(1.9123745) q[2];
rz(-1.2980596) q[3];
sx q[3];
rz(-1.3763873) q[3];
sx q[3];
rz(-1.1840597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.9640279) q[0];
sx q[0];
rz(-2.8819045) q[0];
sx q[0];
rz(2.2143256) q[0];
rz(-0.19993965) q[1];
sx q[1];
rz(-1.6688469) q[1];
sx q[1];
rz(0.98958579) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7829076) q[0];
sx q[0];
rz(-1.9316422) q[0];
sx q[0];
rz(0.11519687) q[0];
rz(-pi) q[1];
rz(0.43439718) q[2];
sx q[2];
rz(-1.6100307) q[2];
sx q[2];
rz(2.5920282) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0285988) q[1];
sx q[1];
rz(-1.494164) q[1];
sx q[1];
rz(-3.0930685) q[1];
rz(-pi) q[2];
rz(2.5603676) q[3];
sx q[3];
rz(-0.94454256) q[3];
sx q[3];
rz(-3.0581829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2135311) q[2];
sx q[2];
rz(-1.4175043) q[2];
sx q[2];
rz(-0.35378635) q[2];
rz(2.1900322) q[3];
sx q[3];
rz(-2.3944201) q[3];
sx q[3];
rz(-2.814754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3554409) q[0];
sx q[0];
rz(-0.94796258) q[0];
sx q[0];
rz(1.0107262) q[0];
rz(3.082869) q[1];
sx q[1];
rz(-1.6207691) q[1];
sx q[1];
rz(-0.40930632) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3407281) q[0];
sx q[0];
rz(-1.1363863) q[0];
sx q[0];
rz(0.77462642) q[0];
x q[1];
rz(1.9366197) q[2];
sx q[2];
rz(-0.43667781) q[2];
sx q[2];
rz(-0.15966378) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9524556) q[1];
sx q[1];
rz(-0.5824832) q[1];
sx q[1];
rz(2.7791695) q[1];
rz(-pi) q[2];
x q[2];
rz(1.549996) q[3];
sx q[3];
rz(-2.0796418) q[3];
sx q[3];
rz(1.5326981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.056444082) q[2];
sx q[2];
rz(-1.2480382) q[2];
sx q[2];
rz(-0.15963456) q[2];
rz(-1.8917482) q[3];
sx q[3];
rz(-1.4188473) q[3];
sx q[3];
rz(-1.0861402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98642629) q[0];
sx q[0];
rz(-2.7041589) q[0];
sx q[0];
rz(2.0470108) q[0];
rz(1.1711586) q[1];
sx q[1];
rz(-1.1181701) q[1];
sx q[1];
rz(-1.0135244) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6863368) q[0];
sx q[0];
rz(-1.6091954) q[0];
sx q[0];
rz(-3.0976035) q[0];
x q[1];
rz(2.065218) q[2];
sx q[2];
rz(-1.567588) q[2];
sx q[2];
rz(0.50960827) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2029548) q[1];
sx q[1];
rz(-1.5302916) q[1];
sx q[1];
rz(0.26085965) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1356527) q[3];
sx q[3];
rz(-1.6683104) q[3];
sx q[3];
rz(1.9103736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1363498) q[2];
sx q[2];
rz(-1.3815657) q[2];
sx q[2];
rz(1.8278149) q[2];
rz(0.93787307) q[3];
sx q[3];
rz(-1.8504986) q[3];
sx q[3];
rz(-3.042799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2336642) q[0];
sx q[0];
rz(-1.6133244) q[0];
sx q[0];
rz(-2.0966356) q[0];
rz(0.16432556) q[1];
sx q[1];
rz(-1.798809) q[1];
sx q[1];
rz(-2.9716861) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63498492) q[0];
sx q[0];
rz(-1.3847376) q[0];
sx q[0];
rz(-3.0691514) q[0];
rz(-0.33268945) q[2];
sx q[2];
rz(-1.2342808) q[2];
sx q[2];
rz(-0.83007407) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.53105866) q[1];
sx q[1];
rz(-1.8137168) q[1];
sx q[1];
rz(-2.2353735) q[1];
x q[2];
rz(-0.14753647) q[3];
sx q[3];
rz(-0.73158598) q[3];
sx q[3];
rz(-0.1583038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.56875151) q[2];
sx q[2];
rz(-0.70256394) q[2];
sx q[2];
rz(2.11002) q[2];
rz(2.0555563) q[3];
sx q[3];
rz(-0.7998172) q[3];
sx q[3];
rz(-1.4097376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80424911) q[0];
sx q[0];
rz(-2.0954837) q[0];
sx q[0];
rz(-1.401249) q[0];
rz(-0.42964545) q[1];
sx q[1];
rz(-2.255217) q[1];
sx q[1];
rz(1.3053798) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0383549) q[0];
sx q[0];
rz(-1.2631386) q[0];
sx q[0];
rz(-2.3820113) q[0];
rz(-pi) q[1];
x q[1];
rz(1.13187) q[2];
sx q[2];
rz(-1.5749075) q[2];
sx q[2];
rz(2.4356349) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4652047) q[1];
sx q[1];
rz(-2.891245) q[1];
sx q[1];
rz(-3.0167104) q[1];
x q[2];
rz(2.1277818) q[3];
sx q[3];
rz(-0.28388043) q[3];
sx q[3];
rz(-0.59040961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1244916) q[2];
sx q[2];
rz(-1.2064826) q[2];
sx q[2];
rz(-0.98480946) q[2];
rz(1.5752327) q[3];
sx q[3];
rz(-1.5896348) q[3];
sx q[3];
rz(0.28520939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.259909) q[0];
sx q[0];
rz(-0.30170983) q[0];
sx q[0];
rz(1.786422) q[0];
rz(-0.024070865) q[1];
sx q[1];
rz(-1.5211952) q[1];
sx q[1];
rz(-2.7640061) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.628869) q[0];
sx q[0];
rz(-2.2170904) q[0];
sx q[0];
rz(1.3747526) q[0];
rz(-pi) q[1];
rz(0.17654769) q[2];
sx q[2];
rz(-1.7804885) q[2];
sx q[2];
rz(3.049946) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8319893) q[1];
sx q[1];
rz(-1.1909232) q[1];
sx q[1];
rz(1.1516476) q[1];
x q[2];
rz(0.92711512) q[3];
sx q[3];
rz(-2.6171631) q[3];
sx q[3];
rz(-0.35628179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54921237) q[2];
sx q[2];
rz(-0.33752957) q[2];
sx q[2];
rz(2.737992) q[2];
rz(-2.9523383) q[3];
sx q[3];
rz(-2.0892102) q[3];
sx q[3];
rz(-0.45690593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52983487) q[0];
sx q[0];
rz(-2.5921322) q[0];
sx q[0];
rz(-2.8305565) q[0];
rz(0.95651904) q[1];
sx q[1];
rz(-1.4776968) q[1];
sx q[1];
rz(-0.32435736) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7081147) q[0];
sx q[0];
rz(-1.38033) q[0];
sx q[0];
rz(0.51149997) q[0];
rz(3.0359984) q[2];
sx q[2];
rz(-2.6289231) q[2];
sx q[2];
rz(1.6216506) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.94453401) q[1];
sx q[1];
rz(-2.1558216) q[1];
sx q[1];
rz(-2.6385175) q[1];
rz(2.764176) q[3];
sx q[3];
rz(-2.2126865) q[3];
sx q[3];
rz(2.9270594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4057464) q[2];
sx q[2];
rz(-0.79068557) q[2];
sx q[2];
rz(-2.9131367) q[2];
rz(-0.53260803) q[3];
sx q[3];
rz(-2.3753128) q[3];
sx q[3];
rz(0.84958357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5936977) q[0];
sx q[0];
rz(-0.51012796) q[0];
sx q[0];
rz(1.8713895) q[0];
rz(-2.5529329) q[1];
sx q[1];
rz(-1.9344067) q[1];
sx q[1];
rz(-0.38280815) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24341881) q[0];
sx q[0];
rz(-0.23056689) q[0];
sx q[0];
rz(-1.4617993) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8201939) q[2];
sx q[2];
rz(-1.1487085) q[2];
sx q[2];
rz(-2.7673134) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8296903) q[1];
sx q[1];
rz(-1.6964579) q[1];
sx q[1];
rz(-0.71255334) q[1];
rz(1.0939264) q[3];
sx q[3];
rz(-2.5989957) q[3];
sx q[3];
rz(2.2063856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1560912) q[2];
sx q[2];
rz(-1.5727377) q[2];
sx q[2];
rz(2.7413979) q[2];
rz(2.8943446) q[3];
sx q[3];
rz(-1.6797545) q[3];
sx q[3];
rz(-0.39321536) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3006712) q[0];
sx q[0];
rz(-1.8074169) q[0];
sx q[0];
rz(-2.1260496) q[0];
rz(0.66889846) q[1];
sx q[1];
rz(-1.4543507) q[1];
sx q[1];
rz(-1.5060172) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43011623) q[0];
sx q[0];
rz(-1.073494) q[0];
sx q[0];
rz(0.74691746) q[0];
x q[1];
rz(2.593562) q[2];
sx q[2];
rz(-1.9491674) q[2];
sx q[2];
rz(-3.1032094) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4989711) q[1];
sx q[1];
rz(-1.6097415) q[1];
sx q[1];
rz(-1.0907111) q[1];
rz(2.5745939) q[3];
sx q[3];
rz(-2.1666489) q[3];
sx q[3];
rz(1.8351549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16613913) q[2];
sx q[2];
rz(-1.0756805) q[2];
sx q[2];
rz(-1.5258741) q[2];
rz(-0.29159355) q[3];
sx q[3];
rz(-0.44277954) q[3];
sx q[3];
rz(2.7899138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6780554) q[0];
sx q[0];
rz(-0.96374496) q[0];
sx q[0];
rz(-2.5441334) q[0];
rz(-2.8529104) q[1];
sx q[1];
rz(-0.79816993) q[1];
sx q[1];
rz(0.023963902) q[1];
rz(0.63715061) q[2];
sx q[2];
rz(-2.646614) q[2];
sx q[2];
rz(2.6029233) q[2];
rz(-1.8221832) q[3];
sx q[3];
rz(-2.0181927) q[3];
sx q[3];
rz(-0.041139091) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
