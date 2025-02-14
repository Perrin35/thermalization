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
rz(0.25207818) q[0];
sx q[0];
rz(-2.5660388) q[0];
sx q[0];
rz(-0.12230305) q[0];
rz(1.1072493) q[1];
sx q[1];
rz(-0.9613494) q[1];
sx q[1];
rz(-0.038318757) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81403804) q[0];
sx q[0];
rz(-2.4310242) q[0];
sx q[0];
rz(-2.1585805) q[0];
rz(0.11876583) q[2];
sx q[2];
rz(-1.6890172) q[2];
sx q[2];
rz(1.880065) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3342383) q[1];
sx q[1];
rz(-1.0206725) q[1];
sx q[1];
rz(-0.38205876) q[1];
rz(-pi) q[2];
rz(1.5502717) q[3];
sx q[3];
rz(-1.8181268) q[3];
sx q[3];
rz(-2.1065245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.64500874) q[2];
sx q[2];
rz(-1.4270447) q[2];
sx q[2];
rz(1.4487779) q[2];
rz(2.951176) q[3];
sx q[3];
rz(-2.0864291) q[3];
sx q[3];
rz(-2.9922488) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9377624) q[0];
sx q[0];
rz(-1.8730524) q[0];
sx q[0];
rz(-0.25800905) q[0];
rz(1.5915271) q[1];
sx q[1];
rz(-1.240088) q[1];
sx q[1];
rz(-2.9749427) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.193522) q[0];
sx q[0];
rz(-1.4753454) q[0];
sx q[0];
rz(0.58346099) q[0];
rz(-pi) q[1];
rz(-2.5160997) q[2];
sx q[2];
rz(-2.6776367) q[2];
sx q[2];
rz(-2.6044012) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2981603) q[1];
sx q[1];
rz(-0.86151228) q[1];
sx q[1];
rz(1.405836) q[1];
x q[2];
rz(2.5949778) q[3];
sx q[3];
rz(-1.0738465) q[3];
sx q[3];
rz(0.95595804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0698645) q[2];
sx q[2];
rz(-2.0191777) q[2];
sx q[2];
rz(1.2321164) q[2];
rz(2.2169436) q[3];
sx q[3];
rz(-0.93622127) q[3];
sx q[3];
rz(-1.1054976) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.449618) q[0];
sx q[0];
rz(-1.469935) q[0];
sx q[0];
rz(0.0019419226) q[0];
rz(0.034189668) q[1];
sx q[1];
rz(-1.9296153) q[1];
sx q[1];
rz(1.5984104) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7352075) q[0];
sx q[0];
rz(-1.3133366) q[0];
sx q[0];
rz(-2.6763112) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9667186) q[2];
sx q[2];
rz(-0.87475077) q[2];
sx q[2];
rz(-0.093122236) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0405827) q[1];
sx q[1];
rz(-1.2173182) q[1];
sx q[1];
rz(1.6124875) q[1];
rz(-1.8887599) q[3];
sx q[3];
rz(-0.67287669) q[3];
sx q[3];
rz(-1.8102136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2446642) q[2];
sx q[2];
rz(-2.7074773) q[2];
sx q[2];
rz(-1.0373235) q[2];
rz(-2.0364929) q[3];
sx q[3];
rz(-1.6048071) q[3];
sx q[3];
rz(1.1387811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26745519) q[0];
sx q[0];
rz(-0.48478165) q[0];
sx q[0];
rz(-1.4111891) q[0];
rz(-0.63938582) q[1];
sx q[1];
rz(-1.6849898) q[1];
sx q[1];
rz(3.0944518) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1829202) q[0];
sx q[0];
rz(-2.483568) q[0];
sx q[0];
rz(0.47112314) q[0];
x q[1];
rz(-2.7732255) q[2];
sx q[2];
rz(-2.6410111) q[2];
sx q[2];
rz(2.6427302) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6005391) q[1];
sx q[1];
rz(-2.4324634) q[1];
sx q[1];
rz(-2.1313271) q[1];
x q[2];
rz(-0.29508884) q[3];
sx q[3];
rz(-0.86936823) q[3];
sx q[3];
rz(1.6158582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.822927) q[2];
sx q[2];
rz(-1.1456127) q[2];
sx q[2];
rz(-1.9226496) q[2];
rz(2.8259891) q[3];
sx q[3];
rz(-0.054840755) q[3];
sx q[3];
rz(0.1489197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3217992) q[0];
sx q[0];
rz(-2.5892374) q[0];
sx q[0];
rz(-2.7929982) q[0];
rz(-1.8888585) q[1];
sx q[1];
rz(-1.1628954) q[1];
sx q[1];
rz(-1.8399651) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063066517) q[0];
sx q[0];
rz(-1.6498897) q[0];
sx q[0];
rz(-2.3344759) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61321494) q[2];
sx q[2];
rz(-0.65468951) q[2];
sx q[2];
rz(1.6729205) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.0685454) q[1];
sx q[1];
rz(-1.3497735) q[1];
sx q[1];
rz(0.18494341) q[1];
rz(-pi) q[2];
rz(-1.0318448) q[3];
sx q[3];
rz(-2.5205118) q[3];
sx q[3];
rz(0.10824848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9584413) q[2];
sx q[2];
rz(-2.8330467) q[2];
sx q[2];
rz(-1.4754254) q[2];
rz(2.4557377) q[3];
sx q[3];
rz(-2.1619022) q[3];
sx q[3];
rz(0.53387749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0297861) q[0];
sx q[0];
rz(-2.989558) q[0];
sx q[0];
rz(-1.6492122) q[0];
rz(-0.97914186) q[1];
sx q[1];
rz(-1.6056332) q[1];
sx q[1];
rz(-1.917256) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79681841) q[0];
sx q[0];
rz(-0.20032702) q[0];
sx q[0];
rz(-0.15423675) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54927214) q[2];
sx q[2];
rz(-1.3830997) q[2];
sx q[2];
rz(2.930738) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2883207) q[1];
sx q[1];
rz(-1.558312) q[1];
sx q[1];
rz(-0.95697831) q[1];
rz(-pi) q[2];
rz(-0.34747261) q[3];
sx q[3];
rz(-2.1712448) q[3];
sx q[3];
rz(0.97755177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7225723) q[2];
sx q[2];
rz(-1.4973065) q[2];
sx q[2];
rz(-2.2255619) q[2];
rz(-1.3845059) q[3];
sx q[3];
rz(-1.8839096) q[3];
sx q[3];
rz(-0.70703435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1322121) q[0];
sx q[0];
rz(-3.0896602) q[0];
sx q[0];
rz(-0.95426553) q[0];
rz(2.7047899) q[1];
sx q[1];
rz(-1.5780459) q[1];
sx q[1];
rz(-0.11016914) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018579114) q[0];
sx q[0];
rz(-1.8948104) q[0];
sx q[0];
rz(-1.3186245) q[0];
x q[1];
rz(2.4442441) q[2];
sx q[2];
rz(-1.1432054) q[2];
sx q[2];
rz(-2.9128437) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.71437826) q[1];
sx q[1];
rz(-2.7958779) q[1];
sx q[1];
rz(-1.5835254) q[1];
rz(-pi) q[2];
rz(-2.0721335) q[3];
sx q[3];
rz(-1.7946464) q[3];
sx q[3];
rz(-1.7373794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43180141) q[2];
sx q[2];
rz(-3.1352037) q[2];
sx q[2];
rz(-2.1978281) q[2];
rz(-0.56378311) q[3];
sx q[3];
rz(-1.0958593) q[3];
sx q[3];
rz(-1.110466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7875882) q[0];
sx q[0];
rz(-0.28689757) q[0];
sx q[0];
rz(-2.7960844) q[0];
rz(1.0954789) q[1];
sx q[1];
rz(-1.6637207) q[1];
sx q[1];
rz(1.5798205) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6876707) q[0];
sx q[0];
rz(-0.78886388) q[0];
sx q[0];
rz(-3.0973017) q[0];
x q[1];
rz(3.079633) q[2];
sx q[2];
rz(-2.2270348) q[2];
sx q[2];
rz(1.7554754) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.53198481) q[1];
sx q[1];
rz(-2.4188305) q[1];
sx q[1];
rz(-1.0015798) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6437279) q[3];
sx q[3];
rz(-1.5988873) q[3];
sx q[3];
rz(-1.3974578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.76274189) q[2];
sx q[2];
rz(-2.8801535) q[2];
sx q[2];
rz(-1.4307107) q[2];
rz(0.53449574) q[3];
sx q[3];
rz(-1.4662687) q[3];
sx q[3];
rz(2.1889595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6119824) q[0];
sx q[0];
rz(-3.072325) q[0];
sx q[0];
rz(-1.8310504) q[0];
rz(2.2546841) q[1];
sx q[1];
rz(-0.84019089) q[1];
sx q[1];
rz(1.2695674) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5909605) q[0];
sx q[0];
rz(-2.7661588) q[0];
sx q[0];
rz(-1.2143137) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1940895) q[2];
sx q[2];
rz(-2.7393638) q[2];
sx q[2];
rz(2.5144983) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.090987) q[1];
sx q[1];
rz(-0.36105737) q[1];
sx q[1];
rz(2.2772842) q[1];
x q[2];
rz(2.9501602) q[3];
sx q[3];
rz(-2.9342954) q[3];
sx q[3];
rz(-1.8199004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5294007) q[2];
sx q[2];
rz(-1.8058913) q[2];
sx q[2];
rz(0.23294918) q[2];
rz(-1.8565146) q[3];
sx q[3];
rz(-0.67684567) q[3];
sx q[3];
rz(-2.7555833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9832298) q[0];
sx q[0];
rz(-0.60929275) q[0];
sx q[0];
rz(-3.0443211) q[0];
rz(0.80398792) q[1];
sx q[1];
rz(-2.4318047) q[1];
sx q[1];
rz(2.9248617) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50057208) q[0];
sx q[0];
rz(-1.4002698) q[0];
sx q[0];
rz(-3.0855623) q[0];
rz(-0.14171504) q[2];
sx q[2];
rz(-1.5248305) q[2];
sx q[2];
rz(-1.7739997) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4137553) q[1];
sx q[1];
rz(-2.1247132) q[1];
sx q[1];
rz(2.3604849) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9166462) q[3];
sx q[3];
rz(-2.2035363) q[3];
sx q[3];
rz(2.2010397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29297605) q[2];
sx q[2];
rz(-0.28102195) q[2];
sx q[2];
rz(2.6463553) q[2];
rz(-1.5826591) q[3];
sx q[3];
rz(-2.0571183) q[3];
sx q[3];
rz(2.9787279) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0611298) q[0];
sx q[0];
rz(-1.5123788) q[0];
sx q[0];
rz(2.6176591) q[0];
rz(1.0864661) q[1];
sx q[1];
rz(-0.51849425) q[1];
sx q[1];
rz(0.35324221) q[1];
rz(-1.5389961) q[2];
sx q[2];
rz(-1.563579) q[2];
sx q[2];
rz(-2.210571) q[2];
rz(-2.2829655) q[3];
sx q[3];
rz(-1.6391907) q[3];
sx q[3];
rz(-1.0427604) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
