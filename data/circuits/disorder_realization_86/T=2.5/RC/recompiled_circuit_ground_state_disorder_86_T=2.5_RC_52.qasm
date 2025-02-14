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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0595221) q[0];
sx q[0];
rz(-0.86468177) q[0];
sx q[0];
rz(2.6816899) q[0];
rz(-2.1291586) q[2];
sx q[2];
rz(-1.867138) q[2];
sx q[2];
rz(-1.3989965) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.77601469) q[1];
sx q[1];
rz(-1.509359) q[1];
sx q[1];
rz(-3.0173124) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7531597) q[3];
sx q[3];
rz(-2.7653501) q[3];
sx q[3];
rz(1.1188504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1755918) q[2];
sx q[2];
rz(-1.5660183) q[2];
sx q[2];
rz(1.2292181) q[2];
rz(1.843533) q[3];
sx q[3];
rz(-1.7652054) q[3];
sx q[3];
rz(1.1840597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17756473) q[0];
sx q[0];
rz(-0.25968817) q[0];
sx q[0];
rz(-0.92726707) q[0];
rz(-0.19993965) q[1];
sx q[1];
rz(-1.4727458) q[1];
sx q[1];
rz(2.1520069) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9703116) q[0];
sx q[0];
rz(-1.6785445) q[0];
sx q[0];
rz(1.9338444) q[0];
rz(-1.5275498) q[2];
sx q[2];
rz(-2.0048365) q[2];
sx q[2];
rz(-1.0394281) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4633219) q[1];
sx q[1];
rz(-3.0509147) q[1];
sx q[1];
rz(2.1342127) q[1];
rz(0.85736147) q[3];
sx q[3];
rz(-2.0318084) q[3];
sx q[3];
rz(2.0217413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2135311) q[2];
sx q[2];
rz(-1.4175043) q[2];
sx q[2];
rz(0.35378635) q[2];
rz(2.1900322) q[3];
sx q[3];
rz(-0.74717251) q[3];
sx q[3];
rz(2.814754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78615171) q[0];
sx q[0];
rz(-0.94796258) q[0];
sx q[0];
rz(-2.1308664) q[0];
rz(0.058723681) q[1];
sx q[1];
rz(-1.6207691) q[1];
sx q[1];
rz(-2.7322863) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6362013) q[0];
sx q[0];
rz(-0.86544468) q[0];
sx q[0];
rz(-0.58569293) q[0];
rz(1.9818166) q[2];
sx q[2];
rz(-1.7226698) q[2];
sx q[2];
rz(-1.077026) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9524556) q[1];
sx q[1];
rz(-0.5824832) q[1];
sx q[1];
rz(0.36242318) q[1];
rz(-pi) q[2];
rz(-1.549996) q[3];
sx q[3];
rz(-1.0619508) q[3];
sx q[3];
rz(1.5326981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.056444082) q[2];
sx q[2];
rz(-1.2480382) q[2];
sx q[2];
rz(2.9819581) q[2];
rz(-1.8917482) q[3];
sx q[3];
rz(-1.7227453) q[3];
sx q[3];
rz(1.0861402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.98642629) q[0];
sx q[0];
rz(-0.43743375) q[0];
sx q[0];
rz(-2.0470108) q[0];
rz(1.9704341) q[1];
sx q[1];
rz(-1.1181701) q[1];
sx q[1];
rz(1.0135244) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60166925) q[0];
sx q[0];
rz(-3.0832096) q[0];
sx q[0];
rz(-2.4235382) q[0];
x q[1];
rz(1.0763747) q[2];
sx q[2];
rz(-1.567588) q[2];
sx q[2];
rz(2.6319844) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5202433) q[1];
sx q[1];
rz(-1.3101556) q[1];
sx q[1];
rz(1.5288749) q[1];
x q[2];
rz(-1.3427686) q[3];
sx q[3];
rz(-2.6963391) q[3];
sx q[3];
rz(2.5955615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0052428) q[2];
sx q[2];
rz(-1.760027) q[2];
sx q[2];
rz(1.3137777) q[2];
rz(2.2037196) q[3];
sx q[3];
rz(-1.2910941) q[3];
sx q[3];
rz(-3.042799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90792847) q[0];
sx q[0];
rz(-1.6133244) q[0];
sx q[0];
rz(-1.044957) q[0];
rz(-2.9772671) q[1];
sx q[1];
rz(-1.798809) q[1];
sx q[1];
rz(-2.9716861) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2192046) q[0];
sx q[0];
rz(-1.6419852) q[0];
sx q[0];
rz(1.3842596) q[0];
rz(-pi) q[1];
rz(-1.2163148) q[2];
sx q[2];
rz(-1.8841685) q[2];
sx q[2];
rz(-0.8543259) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8039717) q[1];
sx q[1];
rz(-0.70120431) q[1];
sx q[1];
rz(-1.9528725) q[1];
rz(-pi) q[2];
rz(0.72615926) q[3];
sx q[3];
rz(-1.4724331) q[3];
sx q[3];
rz(-1.5226328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5728411) q[2];
sx q[2];
rz(-0.70256394) q[2];
sx q[2];
rz(-2.11002) q[2];
rz(1.0860363) q[3];
sx q[3];
rz(-2.3417754) q[3];
sx q[3];
rz(-1.4097376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80424911) q[0];
sx q[0];
rz(-2.0954837) q[0];
sx q[0];
rz(-1.401249) q[0];
rz(2.7119472) q[1];
sx q[1];
rz(-2.255217) q[1];
sx q[1];
rz(1.3053798) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0383549) q[0];
sx q[0];
rz(-1.2631386) q[0];
sx q[0];
rz(0.75958138) q[0];
x q[1];
rz(-3.1370509) q[2];
sx q[2];
rz(-1.1318739) q[2];
sx q[2];
rz(2.2786841) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.77338058) q[1];
sx q[1];
rz(-1.6016593) q[1];
sx q[1];
rz(-0.24847757) q[1];
rz(-pi) q[2];
rz(-1.8135728) q[3];
sx q[3];
rz(-1.7194028) q[3];
sx q[3];
rz(-2.7000303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1244916) q[2];
sx q[2];
rz(-1.93511) q[2];
sx q[2];
rz(0.98480946) q[2];
rz(1.5752327) q[3];
sx q[3];
rz(-1.5896348) q[3];
sx q[3];
rz(0.28520939) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8816836) q[0];
sx q[0];
rz(-2.8398828) q[0];
sx q[0];
rz(-1.786422) q[0];
rz(0.024070865) q[1];
sx q[1];
rz(-1.5211952) q[1];
sx q[1];
rz(2.7640061) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9474038) q[0];
sx q[0];
rz(-0.67127514) q[0];
sx q[0];
rz(-0.25269521) q[0];
rz(2.965045) q[2];
sx q[2];
rz(-1.7804885) q[2];
sx q[2];
rz(-3.049946) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.43289513) q[1];
sx q[1];
rz(-0.55799498) q[1];
sx q[1];
rz(-2.3466228) q[1];
rz(-pi) q[2];
rz(-2.8074516) q[3];
sx q[3];
rz(-1.982882) q[3];
sx q[3];
rz(-1.0704294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5923803) q[2];
sx q[2];
rz(-2.8040631) q[2];
sx q[2];
rz(0.40360061) q[2];
rz(-2.9523383) q[3];
sx q[3];
rz(-1.0523825) q[3];
sx q[3];
rz(-2.6846867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6117578) q[0];
sx q[0];
rz(-2.5921322) q[0];
sx q[0];
rz(-2.8305565) q[0];
rz(-0.95651904) q[1];
sx q[1];
rz(-1.4776968) q[1];
sx q[1];
rz(0.32435736) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8984025) q[0];
sx q[0];
rz(-2.0721738) q[0];
sx q[0];
rz(1.3531951) q[0];
rz(-pi) q[1];
rz(-3.0359984) q[2];
sx q[2];
rz(-2.6289231) q[2];
sx q[2];
rz(-1.6216506) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1970586) q[1];
sx q[1];
rz(-0.98577104) q[1];
sx q[1];
rz(2.6385175) q[1];
x q[2];
rz(-0.37741669) q[3];
sx q[3];
rz(-0.92890611) q[3];
sx q[3];
rz(-2.9270594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4057464) q[2];
sx q[2];
rz(-0.79068557) q[2];
sx q[2];
rz(2.9131367) q[2];
rz(0.53260803) q[3];
sx q[3];
rz(-2.3753128) q[3];
sx q[3];
rz(2.2920091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
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
rz(-1.207186) q[1];
sx q[1];
rz(-2.7587845) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24341881) q[0];
sx q[0];
rz(-2.9110258) q[0];
sx q[0];
rz(-1.4617993) q[0];
rz(-pi) q[1];
rz(0.43395502) q[2];
sx q[2];
rz(-1.343691) q[2];
sx q[2];
rz(1.3004829) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0267566) q[1];
sx q[1];
rz(-0.72163218) q[1];
sx q[1];
rz(-0.19087692) q[1];
rz(-pi) q[2];
rz(2.8715897) q[3];
sx q[3];
rz(-2.0473891) q[3];
sx q[3];
rz(1.4780413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1560912) q[2];
sx q[2];
rz(-1.5727377) q[2];
sx q[2];
rz(-0.4001948) q[2];
rz(-2.8943446) q[3];
sx q[3];
rz(-1.6797545) q[3];
sx q[3];
rz(-2.7483773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3006712) q[0];
sx q[0];
rz(-1.3341757) q[0];
sx q[0];
rz(1.0155431) q[0];
rz(2.4726942) q[1];
sx q[1];
rz(-1.6872419) q[1];
sx q[1];
rz(1.6355754) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5251314) q[0];
sx q[0];
rz(-2.2716953) q[0];
sx q[0];
rz(-0.67411311) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0066543) q[2];
sx q[2];
rz(-1.0653989) q[2];
sx q[2];
rz(1.3874029) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.64262154) q[1];
sx q[1];
rz(-1.6097415) q[1];
sx q[1];
rz(-1.0907111) q[1];
rz(-pi) q[2];
rz(-0.90088441) q[3];
sx q[3];
rz(-2.3436119) q[3];
sx q[3];
rz(0.98679286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16613913) q[2];
sx q[2];
rz(-2.0659122) q[2];
sx q[2];
rz(-1.6157185) q[2];
rz(2.8499991) q[3];
sx q[3];
rz(-2.6988131) q[3];
sx q[3];
rz(0.35167882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.6780554) q[0];
sx q[0];
rz(-0.96374496) q[0];
sx q[0];
rz(-2.5441334) q[0];
rz(-0.28868227) q[1];
sx q[1];
rz(-2.3434227) q[1];
sx q[1];
rz(-3.1176288) q[1];
rz(-1.8815251) q[2];
sx q[2];
rz(-1.962553) q[2];
sx q[2];
rz(0.16061781) q[2];
rz(2.6816159) q[3];
sx q[3];
rz(-1.3446076) q[3];
sx q[3];
rz(1.4190058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
