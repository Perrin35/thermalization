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
rz(-2.6429472) q[0];
sx q[0];
rz(-0.54583025) q[0];
sx q[0];
rz(0.11830615) q[0];
rz(0.91953295) q[1];
sx q[1];
rz(-1.8383205) q[1];
sx q[1];
rz(-1.2143171) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8144381) q[0];
sx q[0];
rz(-2.0025578) q[0];
sx q[0];
rz(-2.9080954) q[0];
rz(-pi) q[1];
rz(2.1636398) q[2];
sx q[2];
rz(-1.0314157) q[2];
sx q[2];
rz(-3.116089) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2976297) q[1];
sx q[1];
rz(-1.7073892) q[1];
sx q[1];
rz(0.94874391) q[1];
rz(1.494058) q[3];
sx q[3];
rz(-0.27394781) q[3];
sx q[3];
rz(2.5016967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.168657) q[2];
sx q[2];
rz(-1.1309705) q[2];
sx q[2];
rz(-1.0966148) q[2];
rz(-2.9003669) q[3];
sx q[3];
rz(-1.7719519) q[3];
sx q[3];
rz(0.84478861) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0519401) q[0];
sx q[0];
rz(-1.0900499) q[0];
sx q[0];
rz(-2.7570778) q[0];
rz(0.3081201) q[1];
sx q[1];
rz(-2.6607951) q[1];
sx q[1];
rz(-2.5164129) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2848672) q[0];
sx q[0];
rz(-1.856904) q[0];
sx q[0];
rz(1.5377783) q[0];
rz(-pi) q[1];
rz(1.6982295) q[2];
sx q[2];
rz(-2.2837167) q[2];
sx q[2];
rz(2.8486203) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5601666) q[1];
sx q[1];
rz(-2.4871965) q[1];
sx q[1];
rz(0.1398211) q[1];
x q[2];
rz(-0.92722882) q[3];
sx q[3];
rz(-0.72690287) q[3];
sx q[3];
rz(2.79799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.008931) q[2];
sx q[2];
rz(-2.2098358) q[2];
sx q[2];
rz(-2.8399732) q[2];
rz(-0.53145069) q[3];
sx q[3];
rz(-2.2558432) q[3];
sx q[3];
rz(-2.0187995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5136435) q[0];
sx q[0];
rz(-1.0114089) q[0];
sx q[0];
rz(-1.1847786) q[0];
rz(-0.14492497) q[1];
sx q[1];
rz(-1.2624319) q[1];
sx q[1];
rz(-1.5473993) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74565174) q[0];
sx q[0];
rz(-2.3303623) q[0];
sx q[0];
rz(2.439789) q[0];
rz(-pi) q[1];
rz(2.9582439) q[2];
sx q[2];
rz(-1.1969222) q[2];
sx q[2];
rz(0.093106769) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.84128296) q[1];
sx q[1];
rz(-2.1108859) q[1];
sx q[1];
rz(-0.45068355) q[1];
x q[2];
rz(-1.8933442) q[3];
sx q[3];
rz(-1.1841906) q[3];
sx q[3];
rz(-0.69535386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3323815) q[2];
sx q[2];
rz(-2.1766162) q[2];
sx q[2];
rz(1.8154209) q[2];
rz(0.91340804) q[3];
sx q[3];
rz(-1.923442) q[3];
sx q[3];
rz(-0.53328812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42295414) q[0];
sx q[0];
rz(-0.47684968) q[0];
sx q[0];
rz(1.3612716) q[0];
rz(2.4286229) q[1];
sx q[1];
rz(-2.0013335) q[1];
sx q[1];
rz(-2.4580809) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8585629) q[0];
sx q[0];
rz(-0.57335317) q[0];
sx q[0];
rz(-1.552215) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2095378) q[2];
sx q[2];
rz(-1.7895797) q[2];
sx q[2];
rz(1.0971341) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5826368) q[1];
sx q[1];
rz(-0.92933944) q[1];
sx q[1];
rz(0.90822405) q[1];
rz(-2.8972647) q[3];
sx q[3];
rz(-1.667898) q[3];
sx q[3];
rz(-1.8116784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5061364) q[2];
sx q[2];
rz(-2.7244302) q[2];
sx q[2];
rz(0.48640856) q[2];
rz(-0.5168612) q[3];
sx q[3];
rz(-1.1812527) q[3];
sx q[3];
rz(1.0335056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085623398) q[0];
sx q[0];
rz(-1.509868) q[0];
sx q[0];
rz(1.7686718) q[0];
rz(-1.7347451) q[1];
sx q[1];
rz(-2.5739248) q[1];
sx q[1];
rz(-2.6645606) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6179028) q[0];
sx q[0];
rz(-1.7672024) q[0];
sx q[0];
rz(1.1454885) q[0];
rz(-pi) q[1];
rz(1.050694) q[2];
sx q[2];
rz(-1.2251266) q[2];
sx q[2];
rz(-0.006397853) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7030695) q[1];
sx q[1];
rz(-2.065899) q[1];
sx q[1];
rz(-0.34214051) q[1];
rz(-pi) q[2];
rz(2.2955016) q[3];
sx q[3];
rz(-2.2109004) q[3];
sx q[3];
rz(-2.7654331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0341952) q[2];
sx q[2];
rz(-2.154978) q[2];
sx q[2];
rz(-2.6743215) q[2];
rz(-0.16921903) q[3];
sx q[3];
rz(-1.7305814) q[3];
sx q[3];
rz(-0.30538487) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6071103) q[0];
sx q[0];
rz(-2.2619673) q[0];
sx q[0];
rz(2.95209) q[0];
rz(-2.8684008) q[1];
sx q[1];
rz(-2.0676282) q[1];
sx q[1];
rz(0.80064076) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9439745) q[0];
sx q[0];
rz(-1.5276143) q[0];
sx q[0];
rz(1.3373242) q[0];
rz(-pi) q[1];
rz(1.6612946) q[2];
sx q[2];
rz(-1.3234089) q[2];
sx q[2];
rz(0.43898459) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2792017) q[1];
sx q[1];
rz(-1.6995879) q[1];
sx q[1];
rz(1.3736082) q[1];
x q[2];
rz(-1.0687105) q[3];
sx q[3];
rz(-1.0866345) q[3];
sx q[3];
rz(-1.7418777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.93929401) q[2];
sx q[2];
rz(-0.17387986) q[2];
sx q[2];
rz(0.536971) q[2];
rz(1.2823074) q[3];
sx q[3];
rz(-1.8351646) q[3];
sx q[3];
rz(-1.4793226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.031161664) q[0];
sx q[0];
rz(-2.3261676) q[0];
sx q[0];
rz(1.331331) q[0];
rz(-1.4415461) q[1];
sx q[1];
rz(-1.4794289) q[1];
sx q[1];
rz(1.9190681) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0095897) q[0];
sx q[0];
rz(-1.9870269) q[0];
sx q[0];
rz(2.8819109) q[0];
rz(0.55180727) q[2];
sx q[2];
rz(-1.6344317) q[2];
sx q[2];
rz(2.3752874) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.26204208) q[1];
sx q[1];
rz(-2.4051504) q[1];
sx q[1];
rz(1.8950736) q[1];
x q[2];
rz(-3.0223373) q[3];
sx q[3];
rz(-1.1290765) q[3];
sx q[3];
rz(0.91631232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9407201) q[2];
sx q[2];
rz(-2.6437289) q[2];
sx q[2];
rz(-2.8087924) q[2];
rz(0.29916549) q[3];
sx q[3];
rz(-1.9009512) q[3];
sx q[3];
rz(3.1201194) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12896319) q[0];
sx q[0];
rz(-2.3733932) q[0];
sx q[0];
rz(0.27498883) q[0];
rz(-3.0601652) q[1];
sx q[1];
rz(-2.3402201) q[1];
sx q[1];
rz(1.8191232) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24182651) q[0];
sx q[0];
rz(-1.6303008) q[0];
sx q[0];
rz(-1.8023947) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0103578) q[2];
sx q[2];
rz(-0.73652525) q[2];
sx q[2];
rz(-0.16939746) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.15582196) q[1];
sx q[1];
rz(-0.81822526) q[1];
sx q[1];
rz(-2.6236412) q[1];
rz(-2.4422936) q[3];
sx q[3];
rz(-0.85996503) q[3];
sx q[3];
rz(-3.1131203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.72121173) q[2];
sx q[2];
rz(-3.0215441) q[2];
sx q[2];
rz(-1.7397286) q[2];
rz(1.8574572) q[3];
sx q[3];
rz(-1.1893585) q[3];
sx q[3];
rz(-3.106626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5818802) q[0];
sx q[0];
rz(-1.5927097) q[0];
sx q[0];
rz(-2.6522719) q[0];
rz(-3.1210461) q[1];
sx q[1];
rz(-1.2362044) q[1];
sx q[1];
rz(2.4403341) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1398292) q[0];
sx q[0];
rz(-2.1368623) q[0];
sx q[0];
rz(-0.28621121) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5145153) q[2];
sx q[2];
rz(-0.70437925) q[2];
sx q[2];
rz(3.0075108) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3534782) q[1];
sx q[1];
rz(-0.81190049) q[1];
sx q[1];
rz(-2.8142156) q[1];
rz(-0.55072983) q[3];
sx q[3];
rz(-2.3756873) q[3];
sx q[3];
rz(-0.85635161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4150841) q[2];
sx q[2];
rz(-2.5231611) q[2];
sx q[2];
rz(-2.6611967) q[2];
rz(1.4091617) q[3];
sx q[3];
rz(-1.7879854) q[3];
sx q[3];
rz(-2.8044146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82520634) q[0];
sx q[0];
rz(-1.9012863) q[0];
sx q[0];
rz(0.41148841) q[0];
rz(-1.8972137) q[1];
sx q[1];
rz(-2.0020516) q[1];
sx q[1];
rz(-2.8172475) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8702983) q[0];
sx q[0];
rz(-2.76519) q[0];
sx q[0];
rz(1.4270622) q[0];
rz(-1.2191992) q[2];
sx q[2];
rz(-1.2484387) q[2];
sx q[2];
rz(-0.11478648) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.38144073) q[1];
sx q[1];
rz(-0.57059971) q[1];
sx q[1];
rz(-1.4165611) q[1];
rz(3.0203811) q[3];
sx q[3];
rz(-1.5838266) q[3];
sx q[3];
rz(-2.7708294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1856498) q[2];
sx q[2];
rz(-1.1142542) q[2];
sx q[2];
rz(-3.0246217) q[2];
rz(2.1864435) q[3];
sx q[3];
rz(-2.9626466) q[3];
sx q[3];
rz(2.5999033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87880001) q[0];
sx q[0];
rz(-0.87651064) q[0];
sx q[0];
rz(0.87690092) q[0];
rz(-1.1204489) q[1];
sx q[1];
rz(-2.3255377) q[1];
sx q[1];
rz(-1.9768523) q[1];
rz(-1.2016313) q[2];
sx q[2];
rz(-2.6045447) q[2];
sx q[2];
rz(2.4640026) q[2];
rz(3.1163327) q[3];
sx q[3];
rz(-1.3002445) q[3];
sx q[3];
rz(0.04095622) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
