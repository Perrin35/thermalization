OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3333617) q[0];
sx q[0];
rz(-0.97167492) q[0];
sx q[0];
rz(1.6605641) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(-0.08637698) q[1];
sx q[1];
rz(0.02286214) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1296596) q[0];
sx q[0];
rz(-1.4225905) q[0];
sx q[0];
rz(-3.0236142) q[0];
x q[1];
rz(2.0055662) q[2];
sx q[2];
rz(-2.044894) q[2];
sx q[2];
rz(0.20027645) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.37576807) q[1];
sx q[1];
rz(-1.7205015) q[1];
sx q[1];
rz(-1.9181812) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95211897) q[3];
sx q[3];
rz(-1.7784356) q[3];
sx q[3];
rz(1.4116956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7227398) q[2];
sx q[2];
rz(-2.3143694) q[2];
sx q[2];
rz(-2.375405) q[2];
rz(2.9700759) q[3];
sx q[3];
rz(-2.4092716) q[3];
sx q[3];
rz(2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8121174) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(0.086659327) q[0];
rz(-2.2791729) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(-3.0564953) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25227308) q[0];
sx q[0];
rz(-0.62250455) q[0];
sx q[0];
rz(3.1134393) q[0];
x q[1];
rz(2.1727824) q[2];
sx q[2];
rz(-1.3170871) q[2];
sx q[2];
rz(1.837455) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.47463402) q[1];
sx q[1];
rz(-2.9999795) q[1];
sx q[1];
rz(0.9160362) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1477633) q[3];
sx q[3];
rz(-0.54518632) q[3];
sx q[3];
rz(2.546437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.26248419) q[2];
sx q[2];
rz(-1.5905249) q[2];
sx q[2];
rz(-1.6764199) q[2];
rz(-1.1761459) q[3];
sx q[3];
rz(-2.2361103) q[3];
sx q[3];
rz(-0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48433205) q[0];
sx q[0];
rz(-0.13467877) q[0];
sx q[0];
rz(-2.235967) q[0];
rz(-2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(1.7571626) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15554409) q[0];
sx q[0];
rz(-1.6533924) q[0];
sx q[0];
rz(2.6655469) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92817523) q[2];
sx q[2];
rz(-2.536142) q[2];
sx q[2];
rz(1.5290934) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3570909) q[1];
sx q[1];
rz(-0.75575268) q[1];
sx q[1];
rz(-3.109039) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82377394) q[3];
sx q[3];
rz(-0.71411055) q[3];
sx q[3];
rz(-0.87574524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8973792) q[2];
sx q[2];
rz(-0.34003568) q[2];
sx q[2];
rz(1.8133694) q[2];
rz(1.7437079) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(0.11169294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(3.0760647) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(-2.4705825) q[0];
rz(-1.3440075) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(2.9342594) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5642693) q[0];
sx q[0];
rz(-0.042357001) q[0];
sx q[0];
rz(-1.3576515) q[0];
x q[1];
rz(-1.7447128) q[2];
sx q[2];
rz(-1.5408278) q[2];
sx q[2];
rz(0.44142351) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.15370788) q[1];
sx q[1];
rz(-2.8578836) q[1];
sx q[1];
rz(0.99516408) q[1];
rz(-0.51165032) q[3];
sx q[3];
rz(-1.5582065) q[3];
sx q[3];
rz(-2.1714641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.31071445) q[2];
sx q[2];
rz(-1.0881311) q[2];
sx q[2];
rz(0.80835289) q[2];
rz(2.3338142) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(-0.19367735) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59072524) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(-2.2461058) q[0];
rz(0.26516178) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(2.6079544) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14347178) q[0];
sx q[0];
rz(-0.40092418) q[0];
sx q[0];
rz(-0.801416) q[0];
rz(-2.9075378) q[2];
sx q[2];
rz(-2.3599527) q[2];
sx q[2];
rz(-2.3603338) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4662019) q[1];
sx q[1];
rz(-1.7168105) q[1];
sx q[1];
rz(-0.22202613) q[1];
x q[2];
rz(0.94138937) q[3];
sx q[3];
rz(-1.7377059) q[3];
sx q[3];
rz(-3.133579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2395893) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(-0.56838244) q[2];
rz(1.9249632) q[3];
sx q[3];
rz(-2.670848) q[3];
sx q[3];
rz(-2.7980428) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50826532) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(-1.9653962) q[0];
rz(-1.5856702) q[1];
sx q[1];
rz(-0.95247477) q[1];
sx q[1];
rz(1.0046545) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83197901) q[0];
sx q[0];
rz(-1.2011315) q[0];
sx q[0];
rz(-1.9395104) q[0];
rz(-1.651628) q[2];
sx q[2];
rz(-1.1751375) q[2];
sx q[2];
rz(0.8880907) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5914454) q[1];
sx q[1];
rz(-1.0284371) q[1];
sx q[1];
rz(-1.28246) q[1];
rz(-1.5980914) q[3];
sx q[3];
rz(-0.41397646) q[3];
sx q[3];
rz(2.9881145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4040318) q[2];
sx q[2];
rz(-0.99016756) q[2];
sx q[2];
rz(-0.44815865) q[2];
rz(-0.29799497) q[3];
sx q[3];
rz(-2.4500676) q[3];
sx q[3];
rz(2.7929849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3649243) q[0];
sx q[0];
rz(-1.8402599) q[0];
sx q[0];
rz(-0.65814322) q[0];
rz(1.2843885) q[1];
sx q[1];
rz(-2.6874266) q[1];
sx q[1];
rz(-2.4694494) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2015398) q[0];
sx q[0];
rz(-1.7367898) q[0];
sx q[0];
rz(0.0075017651) q[0];
x q[1];
rz(-2.0124112) q[2];
sx q[2];
rz(-1.1814983) q[2];
sx q[2];
rz(-1.0600952) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0587412) q[1];
sx q[1];
rz(-0.62742678) q[1];
sx q[1];
rz(-2.4861366) q[1];
rz(-1.7147195) q[3];
sx q[3];
rz(-2.2669499) q[3];
sx q[3];
rz(-1.6203984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.60823524) q[2];
sx q[2];
rz(-2.50864) q[2];
sx q[2];
rz(-2.7427924) q[2];
rz(2.3170025) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(-0.59857541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44052112) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(-2.9833802) q[0];
rz(-1.0868866) q[1];
sx q[1];
rz(-0.67651665) q[1];
sx q[1];
rz(-0.80668443) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2465625) q[0];
sx q[0];
rz(-2.6916305) q[0];
sx q[0];
rz(1.297784) q[0];
rz(0.20571795) q[2];
sx q[2];
rz(-1.6763902) q[2];
sx q[2];
rz(-0.92668698) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.403277) q[1];
sx q[1];
rz(-2.3936831) q[1];
sx q[1];
rz(2.7110389) q[1];
rz(-pi) q[2];
rz(-1.673647) q[3];
sx q[3];
rz(-0.30234435) q[3];
sx q[3];
rz(-2.9000988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.57758254) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(-2.1419443) q[2];
rz(-3.0380761) q[3];
sx q[3];
rz(-0.13705702) q[3];
sx q[3];
rz(2.0172393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0209811) q[0];
sx q[0];
rz(-0.7779026) q[0];
sx q[0];
rz(-2.4840684) q[0];
rz(-2.855037) q[1];
sx q[1];
rz(-2.2189238) q[1];
sx q[1];
rz(0.34067571) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0176312) q[0];
sx q[0];
rz(-0.88273662) q[0];
sx q[0];
rz(-2.7244363) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1358143) q[2];
sx q[2];
rz(-1.3823595) q[2];
sx q[2];
rz(0.68067683) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9790736) q[1];
sx q[1];
rz(-0.58774978) q[1];
sx q[1];
rz(1.2390922) q[1];
x q[2];
rz(-2.9433603) q[3];
sx q[3];
rz(-1.7361589) q[3];
sx q[3];
rz(-2.4076622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.75491536) q[2];
sx q[2];
rz(-1.200054) q[2];
sx q[2];
rz(-0.49794751) q[2];
rz(-2.8399816) q[3];
sx q[3];
rz(-0.40062723) q[3];
sx q[3];
rz(2.6224459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72934812) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(0.27858946) q[0];
rz(2.5623698) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(-3.0648807) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3007606) q[0];
sx q[0];
rz(-1.4474086) q[0];
sx q[0];
rz(1.8873909) q[0];
x q[1];
rz(2.9701091) q[2];
sx q[2];
rz(-2.137261) q[2];
sx q[2];
rz(-0.4085853) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6381166) q[1];
sx q[1];
rz(-1.5917115) q[1];
sx q[1];
rz(-0.36550826) q[1];
rz(-pi) q[2];
x q[2];
rz(1.65927) q[3];
sx q[3];
rz(-1.7469329) q[3];
sx q[3];
rz(-0.74151553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3672093) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(0.45551604) q[2];
rz(0.43141836) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(0.063522696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9615622) q[0];
sx q[0];
rz(-1.4008235) q[0];
sx q[0];
rz(0.93749198) q[0];
rz(0.58615276) q[1];
sx q[1];
rz(-1.502232) q[1];
sx q[1];
rz(-1.4486817) q[1];
rz(2.8557599) q[2];
sx q[2];
rz(-2.0252068) q[2];
sx q[2];
rz(0.037638738) q[2];
rz(1.205411) q[3];
sx q[3];
rz(-1.972651) q[3];
sx q[3];
rz(-2.7046711) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];