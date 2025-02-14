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
rz(1.3156112) q[0];
sx q[0];
rz(5.4551107) q[0];
sx q[0];
rz(8.5760737) q[0];
rz(1.4915713) q[1];
sx q[1];
rz(-0.68869156) q[1];
sx q[1];
rz(0.42660776) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2036982) q[0];
sx q[0];
rz(-0.12591322) q[0];
sx q[0];
rz(2.4975169) q[0];
rz(-3.0399917) q[2];
sx q[2];
rz(-1.6825321) q[2];
sx q[2];
rz(2.8174809) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6411031) q[1];
sx q[1];
rz(-1.4979384) q[1];
sx q[1];
rz(0.44568731) q[1];
rz(0.95152609) q[3];
sx q[3];
rz(-2.2689156) q[3];
sx q[3];
rz(-1.2941456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9244869) q[2];
sx q[2];
rz(-1.3873528) q[2];
sx q[2];
rz(-3.1021049) q[2];
rz(-1.0314137) q[3];
sx q[3];
rz(-2.1125427) q[3];
sx q[3];
rz(1.2379117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51204387) q[0];
sx q[0];
rz(-1.6144253) q[0];
sx q[0];
rz(1.7426096) q[0];
rz(1.6276739) q[1];
sx q[1];
rz(-2.2986423) q[1];
sx q[1];
rz(-1.7248076) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48081452) q[0];
sx q[0];
rz(-2.3153458) q[0];
sx q[0];
rz(-1.2673402) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8292839) q[2];
sx q[2];
rz(-2.333039) q[2];
sx q[2];
rz(1.1201348) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21525225) q[1];
sx q[1];
rz(-1.7071586) q[1];
sx q[1];
rz(-2.7424314) q[1];
x q[2];
rz(-1.6187158) q[3];
sx q[3];
rz(-1.1635492) q[3];
sx q[3];
rz(-1.6919235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.17025718) q[2];
sx q[2];
rz(-1.1032666) q[2];
sx q[2];
rz(-2.4309168) q[2];
rz(-3.1273048) q[3];
sx q[3];
rz(-0.24737869) q[3];
sx q[3];
rz(1.8054731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22064848) q[0];
sx q[0];
rz(-0.9592239) q[0];
sx q[0];
rz(0.4162108) q[0];
rz(-1.2843708) q[1];
sx q[1];
rz(-1.4764079) q[1];
sx q[1];
rz(-2.823337) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98781768) q[0];
sx q[0];
rz(-1.9526228) q[0];
sx q[0];
rz(2.8433958) q[0];
rz(2.2936965) q[2];
sx q[2];
rz(-1.0506786) q[2];
sx q[2];
rz(-2.0682356) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-15/(11*pi)) q[1];
sx q[1];
rz(-2.4726082) q[1];
sx q[1];
rz(-2.9218052) q[1];
rz(0.24782039) q[3];
sx q[3];
rz(-1.2711002) q[3];
sx q[3];
rz(1.1775985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2661813) q[2];
sx q[2];
rz(-0.77989945) q[2];
sx q[2];
rz(-1.8511339) q[2];
rz(-0.053704638) q[3];
sx q[3];
rz(-1.8219681) q[3];
sx q[3];
rz(1.3857589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5829492) q[0];
sx q[0];
rz(-2.4531893) q[0];
sx q[0];
rz(2.131856) q[0];
rz(2.2263777) q[1];
sx q[1];
rz(-1.6123632) q[1];
sx q[1];
rz(-0.42542747) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64176004) q[0];
sx q[0];
rz(-1.4041889) q[0];
sx q[0];
rz(0.84979041) q[0];
rz(-0.82773005) q[2];
sx q[2];
rz(-1.5166188) q[2];
sx q[2];
rz(2.3119213) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.47290651) q[1];
sx q[1];
rz(-0.91835178) q[1];
sx q[1];
rz(-2.621306) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6050379) q[3];
sx q[3];
rz(-2.9786907) q[3];
sx q[3];
rz(2.3690641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2414744) q[2];
sx q[2];
rz(-1.5510617) q[2];
sx q[2];
rz(3.0305064) q[2];
rz(0.19242081) q[3];
sx q[3];
rz(-2.7399053) q[3];
sx q[3];
rz(-2.4470611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7531994) q[0];
sx q[0];
rz(-2.41112) q[0];
sx q[0];
rz(-0.86517349) q[0];
rz(2.6490037) q[1];
sx q[1];
rz(-1.0961696) q[1];
sx q[1];
rz(-1.385484) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2619721) q[0];
sx q[0];
rz(-0.79085717) q[0];
sx q[0];
rz(1.7636931) q[0];
rz(-pi) q[1];
rz(-2.7215459) q[2];
sx q[2];
rz(-0.85424747) q[2];
sx q[2];
rz(1.6796215) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.046987586) q[1];
sx q[1];
rz(-1.054227) q[1];
sx q[1];
rz(2.0265685) q[1];
x q[2];
rz(2.4597857) q[3];
sx q[3];
rz(-2.1958133) q[3];
sx q[3];
rz(2.8375219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8869141) q[2];
sx q[2];
rz(-2.6723537) q[2];
sx q[2];
rz(0.706642) q[2];
rz(0.32736579) q[3];
sx q[3];
rz(-1.2923765) q[3];
sx q[3];
rz(2.4842026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8364328) q[0];
sx q[0];
rz(-1.3494455) q[0];
sx q[0];
rz(-2.7850372) q[0];
rz(0.56218475) q[1];
sx q[1];
rz(-2.0088582) q[1];
sx q[1];
rz(-0.98639375) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0092457744) q[0];
sx q[0];
rz(-2.2725687) q[0];
sx q[0];
rz(2.470507) q[0];
x q[1];
rz(2.7964727) q[2];
sx q[2];
rz(-1.9202931) q[2];
sx q[2];
rz(2.0728759) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.84556678) q[1];
sx q[1];
rz(-2.3289776) q[1];
sx q[1];
rz(-1.6065099) q[1];
rz(-pi) q[2];
rz(-1.8142734) q[3];
sx q[3];
rz(-1.7469353) q[3];
sx q[3];
rz(1.0276664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76546136) q[2];
sx q[2];
rz(-2.4262846) q[2];
sx q[2];
rz(-0.7005271) q[2];
rz(-0.96945196) q[3];
sx q[3];
rz(-1.0337831) q[3];
sx q[3];
rz(-2.2112924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1932061) q[0];
sx q[0];
rz(-0.32873118) q[0];
sx q[0];
rz(0.1513924) q[0];
rz(1.6311749) q[1];
sx q[1];
rz(-2.0163586) q[1];
sx q[1];
rz(-1.6815394) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3356381) q[0];
sx q[0];
rz(-2.3344669) q[0];
sx q[0];
rz(-0.7025848) q[0];
rz(-pi) q[1];
rz(1.0437772) q[2];
sx q[2];
rz(-1.8879381) q[2];
sx q[2];
rz(0.70995599) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9402079) q[1];
sx q[1];
rz(-1.0321069) q[1];
sx q[1];
rz(2.687665) q[1];
rz(-pi) q[2];
rz(-1.7655108) q[3];
sx q[3];
rz(-1.4275744) q[3];
sx q[3];
rz(-2.9778874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2363362) q[2];
sx q[2];
rz(-2.4589804) q[2];
sx q[2];
rz(0.35161099) q[2];
rz(2.6998399) q[3];
sx q[3];
rz(-0.92919246) q[3];
sx q[3];
rz(-0.15171224) q[3];
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
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041895954) q[0];
sx q[0];
rz(-1.8459039) q[0];
sx q[0];
rz(1.9805441) q[0];
rz(0.64758045) q[1];
sx q[1];
rz(-1.7627962) q[1];
sx q[1];
rz(2.3191998) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1402991) q[0];
sx q[0];
rz(-1.1736794) q[0];
sx q[0];
rz(2.1236093) q[0];
x q[1];
rz(-2.0089125) q[2];
sx q[2];
rz(-2.2575976) q[2];
sx q[2];
rz(0.39114726) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5288189) q[1];
sx q[1];
rz(-2.3102323) q[1];
sx q[1];
rz(-2.7465613) q[1];
rz(-pi) q[2];
rz(2.3052182) q[3];
sx q[3];
rz(-0.55787239) q[3];
sx q[3];
rz(-1.316586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9245727) q[2];
sx q[2];
rz(-1.9968888) q[2];
sx q[2];
rz(1.3059957) q[2];
rz(0.71803391) q[3];
sx q[3];
rz(-0.82457232) q[3];
sx q[3];
rz(0.048862783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7056535) q[0];
sx q[0];
rz(-2.0581364) q[0];
sx q[0];
rz(-3.1245681) q[0];
rz(1.1964993) q[1];
sx q[1];
rz(-0.39533177) q[1];
sx q[1];
rz(2.8065525) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25680731) q[0];
sx q[0];
rz(-1.5809158) q[0];
sx q[0];
rz(-1.0027893) q[0];
rz(-pi) q[1];
x q[1];
rz(0.64249383) q[2];
sx q[2];
rz(-1.9064925) q[2];
sx q[2];
rz(3.1106126) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2372053) q[1];
sx q[1];
rz(-0.6977302) q[1];
sx q[1];
rz(-1.8664788) q[1];
rz(-1.0862971) q[3];
sx q[3];
rz(-2.006975) q[3];
sx q[3];
rz(2.0791813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.42387858) q[2];
sx q[2];
rz(-0.7462036) q[2];
sx q[2];
rz(2.8625873) q[2];
rz(-1.3689857) q[3];
sx q[3];
rz(-1.5149346) q[3];
sx q[3];
rz(-2.2122808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6245215) q[0];
sx q[0];
rz(-2.9947424) q[0];
sx q[0];
rz(0.76706925) q[0];
rz(0.37297878) q[1];
sx q[1];
rz(-1.6423128) q[1];
sx q[1];
rz(2.9708718) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45324651) q[0];
sx q[0];
rz(-0.73536721) q[0];
sx q[0];
rz(-3.0985996) q[0];
x q[1];
rz(-0.41497322) q[2];
sx q[2];
rz(-1.4230721) q[2];
sx q[2];
rz(1.8208675) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.77826277) q[1];
sx q[1];
rz(-2.7365541) q[1];
sx q[1];
rz(2.2447849) q[1];
rz(-pi) q[2];
rz(-2.6737718) q[3];
sx q[3];
rz(-1.7138283) q[3];
sx q[3];
rz(0.19090548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3501222) q[2];
sx q[2];
rz(-2.5779724) q[2];
sx q[2];
rz(-3.138809) q[2];
rz(2.2400098) q[3];
sx q[3];
rz(-1.0222579) q[3];
sx q[3];
rz(2.3367052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0881385) q[0];
sx q[0];
rz(-2.8202941) q[0];
sx q[0];
rz(1.1711076) q[0];
rz(0.83012719) q[1];
sx q[1];
rz(-1.8142038) q[1];
sx q[1];
rz(1.289191) q[1];
rz(0.96021658) q[2];
sx q[2];
rz(-1.1398906) q[2];
sx q[2];
rz(0.7791145) q[2];
rz(-2.4975111) q[3];
sx q[3];
rz(-0.91337503) q[3];
sx q[3];
rz(2.5144284) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
