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
rz(-2.6540304) q[0];
sx q[0];
rz(-2.7609479) q[0];
sx q[0];
rz(-0.08925499) q[0];
rz(-0.46269497) q[1];
sx q[1];
rz(-2.8609639) q[1];
sx q[1];
rz(1.7791003) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.270705) q[0];
sx q[0];
rz(-1.1977559) q[0];
sx q[0];
rz(-1.5482658) q[0];
rz(0.2398165) q[2];
sx q[2];
rz(-0.15286286) q[2];
sx q[2];
rz(2.3279362) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.465724) q[1];
sx q[1];
rz(-1.3480519) q[1];
sx q[1];
rz(-2.5506819) q[1];
rz(-pi) q[2];
rz(1.6600577) q[3];
sx q[3];
rz(-1.1162471) q[3];
sx q[3];
rz(-2.4069569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6224299) q[2];
sx q[2];
rz(-0.76202718) q[2];
sx q[2];
rz(-2.2117174) q[2];
rz(-0.30396384) q[3];
sx q[3];
rz(-1.8094939) q[3];
sx q[3];
rz(-0.16873321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25794432) q[0];
sx q[0];
rz(-1.5771663) q[0];
sx q[0];
rz(1.8899348) q[0];
rz(2.9257863) q[1];
sx q[1];
rz(-2.1025175) q[1];
sx q[1];
rz(0.13914093) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.62876) q[0];
sx q[0];
rz(-0.91673512) q[0];
sx q[0];
rz(-2.9399228) q[0];
rz(-0.47454796) q[2];
sx q[2];
rz(-2.5926771) q[2];
sx q[2];
rz(0.35473706) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.58932585) q[1];
sx q[1];
rz(-1.8084452) q[1];
sx q[1];
rz(2.4219805) q[1];
x q[2];
rz(2.3288547) q[3];
sx q[3];
rz(-2.2366282) q[3];
sx q[3];
rz(-2.8259371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4591878) q[2];
sx q[2];
rz(-2.6725957) q[2];
sx q[2];
rz(1.0914618) q[2];
rz(1.5952716) q[3];
sx q[3];
rz(-1.073444) q[3];
sx q[3];
rz(-0.44927621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.062155398) q[0];
sx q[0];
rz(-1.2261483) q[0];
sx q[0];
rz(-0.011818258) q[0];
rz(-0.91105175) q[1];
sx q[1];
rz(-1.751519) q[1];
sx q[1];
rz(-1.8131088) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8810711) q[0];
sx q[0];
rz(-1.0539454) q[0];
sx q[0];
rz(-1.757574) q[0];
x q[1];
rz(-2.6027868) q[2];
sx q[2];
rz(-2.9201815) q[2];
sx q[2];
rz(-0.087269727) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7367626) q[1];
sx q[1];
rz(-1.1859845) q[1];
sx q[1];
rz(0.75508187) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6823009) q[3];
sx q[3];
rz(-1.498025) q[3];
sx q[3];
rz(-2.7804216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.38982424) q[2];
sx q[2];
rz(-1.6168892) q[2];
sx q[2];
rz(-0.41979182) q[2];
rz(-1.7828364) q[3];
sx q[3];
rz(-2.8163781) q[3];
sx q[3];
rz(-1.9512008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.457394) q[0];
sx q[0];
rz(-1.1379108) q[0];
sx q[0];
rz(2.6275291) q[0];
rz(-0.38814107) q[1];
sx q[1];
rz(-2.5297647) q[1];
sx q[1];
rz(-1.911389) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82321754) q[0];
sx q[0];
rz(-1.3127232) q[0];
sx q[0];
rz(2.8283872) q[0];
rz(-2.7204334) q[2];
sx q[2];
rz(-2.371726) q[2];
sx q[2];
rz(1.0854967) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19691218) q[1];
sx q[1];
rz(-2.3397603) q[1];
sx q[1];
rz(-2.6242328) q[1];
rz(-1.7994653) q[3];
sx q[3];
rz(-0.91311753) q[3];
sx q[3];
rz(-0.87008892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6396883) q[2];
sx q[2];
rz(-2.1268667) q[2];
sx q[2];
rz(-2.225303) q[2];
rz(0.4942975) q[3];
sx q[3];
rz(-1.9435792) q[3];
sx q[3];
rz(0.77427197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2275527) q[0];
sx q[0];
rz(-0.95714772) q[0];
sx q[0];
rz(-0.48602948) q[0];
rz(2.5388429) q[1];
sx q[1];
rz(-1.9662247) q[1];
sx q[1];
rz(-2.026162) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1971216) q[0];
sx q[0];
rz(-2.174518) q[0];
sx q[0];
rz(-2.05462) q[0];
x q[1];
rz(-2.3472842) q[2];
sx q[2];
rz(-1.8657078) q[2];
sx q[2];
rz(-0.25943929) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.095905) q[1];
sx q[1];
rz(-2.3153911) q[1];
sx q[1];
rz(-2.5627441) q[1];
rz(-2.8776725) q[3];
sx q[3];
rz(-1.1047398) q[3];
sx q[3];
rz(1.6845244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1439765) q[2];
sx q[2];
rz(-2.82085) q[2];
sx q[2];
rz(2.0799267) q[2];
rz(1.5205421) q[3];
sx q[3];
rz(-1.1989667) q[3];
sx q[3];
rz(-0.89402136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068037085) q[0];
sx q[0];
rz(-3.1017922) q[0];
sx q[0];
rz(-1.918248) q[0];
rz(-2.6262737) q[1];
sx q[1];
rz(-0.66910187) q[1];
sx q[1];
rz(1.7759148) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4752303) q[0];
sx q[0];
rz(-1.7710553) q[0];
sx q[0];
rz(-3.0259575) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0480644) q[2];
sx q[2];
rz(-1.3892738) q[2];
sx q[2];
rz(0.43147555) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.20582304) q[1];
sx q[1];
rz(-2.6403815) q[1];
sx q[1];
rz(-2.8233714) q[1];
rz(-pi) q[2];
rz(-3.0240718) q[3];
sx q[3];
rz(-0.29942808) q[3];
sx q[3];
rz(2.260331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0627275) q[2];
sx q[2];
rz(-1.3389503) q[2];
sx q[2];
rz(-1.3272746) q[2];
rz(1.4158538) q[3];
sx q[3];
rz(-0.073181987) q[3];
sx q[3];
rz(2.2590526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13407229) q[0];
sx q[0];
rz(-2.6044758) q[0];
sx q[0];
rz(0.19675955) q[0];
rz(-2.9255731) q[1];
sx q[1];
rz(-1.0088423) q[1];
sx q[1];
rz(-2.3541727) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8674054) q[0];
sx q[0];
rz(-1.3816773) q[0];
sx q[0];
rz(1.8195137) q[0];
x q[1];
rz(0.49968736) q[2];
sx q[2];
rz(-2.5737107) q[2];
sx q[2];
rz(-2.5497921) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1087139) q[1];
sx q[1];
rz(-1.6163663) q[1];
sx q[1];
rz(2.5497863) q[1];
rz(-pi) q[2];
rz(1.2489572) q[3];
sx q[3];
rz(-1.3909381) q[3];
sx q[3];
rz(-1.091979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.54062033) q[2];
sx q[2];
rz(-1.629402) q[2];
sx q[2];
rz(0.38360325) q[2];
rz(-2.2331734) q[3];
sx q[3];
rz(-2.3288265) q[3];
sx q[3];
rz(0.15373357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28021321) q[0];
sx q[0];
rz(-2.6160243) q[0];
sx q[0];
rz(-2.5049765) q[0];
rz(0.55066291) q[1];
sx q[1];
rz(-1.6267585) q[1];
sx q[1];
rz(-2.6689463) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.11926) q[0];
sx q[0];
rz(-0.68281931) q[0];
sx q[0];
rz(-0.14108087) q[0];
rz(-pi) q[1];
rz(2.6246895) q[2];
sx q[2];
rz(-3.0344525) q[2];
sx q[2];
rz(-0.60984367) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7618666) q[1];
sx q[1];
rz(-1.7032924) q[1];
sx q[1];
rz(-1.6220868) q[1];
rz(2.74284) q[3];
sx q[3];
rz(-0.67091432) q[3];
sx q[3];
rz(0.89445597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4281281) q[2];
sx q[2];
rz(-1.8992004) q[2];
sx q[2];
rz(-0.098527519) q[2];
rz(3.1108917) q[3];
sx q[3];
rz(-1.5701598) q[3];
sx q[3];
rz(-2.9876685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5310265) q[0];
sx q[0];
rz(-0.95777804) q[0];
sx q[0];
rz(2.4089693) q[0];
rz(-0.0064119617) q[1];
sx q[1];
rz(-1.0259722) q[1];
sx q[1];
rz(1.7195255) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.812987) q[0];
sx q[0];
rz(-1.7610794) q[0];
sx q[0];
rz(-1.8237795) q[0];
rz(-pi) q[1];
rz(2.6426475) q[2];
sx q[2];
rz(-1.8932668) q[2];
sx q[2];
rz(2.6889963) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3141503) q[1];
sx q[1];
rz(-1.536751) q[1];
sx q[1];
rz(2.9877547) q[1];
rz(-pi) q[2];
rz(1.645817) q[3];
sx q[3];
rz(-1.0524233) q[3];
sx q[3];
rz(-0.013766001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0942568) q[2];
sx q[2];
rz(-1.7862659) q[2];
sx q[2];
rz(2.6195841) q[2];
rz(3.1212854) q[3];
sx q[3];
rz(-2.4419407) q[3];
sx q[3];
rz(2.844753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8727528) q[0];
sx q[0];
rz(-1.5248542) q[0];
sx q[0];
rz(1.1312477) q[0];
rz(-0.77992431) q[1];
sx q[1];
rz(-1.2029388) q[1];
sx q[1];
rz(0.47526971) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49233046) q[0];
sx q[0];
rz(-1.2075675) q[0];
sx q[0];
rz(-0.80357768) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9067801) q[2];
sx q[2];
rz(-0.91660344) q[2];
sx q[2];
rz(1.393911) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.77716) q[1];
sx q[1];
rz(-1.8231943) q[1];
sx q[1];
rz(-2.9168058) q[1];
rz(-pi) q[2];
rz(2.9679139) q[3];
sx q[3];
rz(-0.70391897) q[3];
sx q[3];
rz(-2.6736163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37107006) q[2];
sx q[2];
rz(-1.5757685) q[2];
sx q[2];
rz(-1.2062262) q[2];
rz(1.3953588) q[3];
sx q[3];
rz(-0.54308707) q[3];
sx q[3];
rz(-0.079308184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4974925) q[0];
sx q[0];
rz(-0.19484367) q[0];
sx q[0];
rz(2.3051443) q[0];
rz(-0.99776987) q[1];
sx q[1];
rz(-2.0496968) q[1];
sx q[1];
rz(-3.0539378) q[1];
rz(1.3523819) q[2];
sx q[2];
rz(-1.6383258) q[2];
sx q[2];
rz(0.93371423) q[2];
rz(-1.6818123) q[3];
sx q[3];
rz(-0.89499499) q[3];
sx q[3];
rz(0.53193308) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
