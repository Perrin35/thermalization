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
rz(0.60740745) q[0];
sx q[0];
rz(2.6976801) q[0];
sx q[0];
rz(9.1023268) q[0];
rz(-1.9880265) q[1];
sx q[1];
rz(-1.1918951) q[1];
sx q[1];
rz(0.21790394) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8810824) q[0];
sx q[0];
rz(-2.5985607) q[0];
sx q[0];
rz(0.61123993) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9548151) q[2];
sx q[2];
rz(-1.412027) q[2];
sx q[2];
rz(-1.1940317) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.17927781) q[1];
sx q[1];
rz(-1.4554108) q[1];
sx q[1];
rz(-2.8820267) q[1];
rz(2.3090906) q[3];
sx q[3];
rz(-0.41838405) q[3];
sx q[3];
rz(2.7184021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0107062) q[2];
sx q[2];
rz(-0.43121269) q[2];
sx q[2];
rz(-0.19403379) q[2];
rz(-2.7583097) q[3];
sx q[3];
rz(-2.2563939) q[3];
sx q[3];
rz(-1.6398199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5782769) q[0];
sx q[0];
rz(-1.0551772) q[0];
sx q[0];
rz(-1.1042327) q[0];
rz(0.68179321) q[1];
sx q[1];
rz(-1.3673404) q[1];
sx q[1];
rz(-1.7368447) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0007219) q[0];
sx q[0];
rz(-2.8551176) q[0];
sx q[0];
rz(2.6764144) q[0];
rz(-2.4656336) q[2];
sx q[2];
rz(-1.1164719) q[2];
sx q[2];
rz(0.76269645) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4841744) q[1];
sx q[1];
rz(-0.60125105) q[1];
sx q[1];
rz(0.57415331) q[1];
x q[2];
rz(2.8867763) q[3];
sx q[3];
rz(-1.4597893) q[3];
sx q[3];
rz(1.0274344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76689395) q[2];
sx q[2];
rz(-2.1465492) q[2];
sx q[2];
rz(-1.2332756) q[2];
rz(0.79902664) q[3];
sx q[3];
rz(-0.34501758) q[3];
sx q[3];
rz(-3.0123805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98053539) q[0];
sx q[0];
rz(-1.072847) q[0];
sx q[0];
rz(-0.55915731) q[0];
rz(2.0237538) q[1];
sx q[1];
rz(-1.585377) q[1];
sx q[1];
rz(-3.0303755) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2905544) q[0];
sx q[0];
rz(-1.1946329) q[0];
sx q[0];
rz(1.4567503) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3005914) q[2];
sx q[2];
rz(-0.96365721) q[2];
sx q[2];
rz(2.6532113) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4065985) q[1];
sx q[1];
rz(-1.7593699) q[1];
sx q[1];
rz(0.85822924) q[1];
rz(-pi) q[2];
rz(1.4239156) q[3];
sx q[3];
rz(-1.960117) q[3];
sx q[3];
rz(-1.9922272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8289651) q[2];
sx q[2];
rz(-1.4451005) q[2];
sx q[2];
rz(-1.1980537) q[2];
rz(1.7184006) q[3];
sx q[3];
rz(-1.6953902) q[3];
sx q[3];
rz(2.6660582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51197416) q[0];
sx q[0];
rz(-0.29933512) q[0];
sx q[0];
rz(-2.3053115) q[0];
rz(1.5760999) q[1];
sx q[1];
rz(-2.0471768) q[1];
sx q[1];
rz(2.2907168) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99573409) q[0];
sx q[0];
rz(-0.99829095) q[0];
sx q[0];
rz(0.03972223) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81340547) q[2];
sx q[2];
rz(-0.69077531) q[2];
sx q[2];
rz(-1.1479614) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5961413) q[1];
sx q[1];
rz(-0.79789466) q[1];
sx q[1];
rz(-1.376838) q[1];
x q[2];
rz(2.0381884) q[3];
sx q[3];
rz(-0.43741832) q[3];
sx q[3];
rz(-1.9269772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0389082) q[2];
sx q[2];
rz(-0.85997283) q[2];
sx q[2];
rz(2.6343708) q[2];
rz(0.66633362) q[3];
sx q[3];
rz(-0.75080502) q[3];
sx q[3];
rz(-2.2255118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0136593) q[0];
sx q[0];
rz(-1.4480696) q[0];
sx q[0];
rz(1.5284982) q[0];
rz(-1.5616034) q[1];
sx q[1];
rz(-2.0785619) q[1];
sx q[1];
rz(2.9772421) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30475475) q[0];
sx q[0];
rz(-1.4753818) q[0];
sx q[0];
rz(-1.6518023) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5098624) q[2];
sx q[2];
rz(-2.1513878) q[2];
sx q[2];
rz(2.3931062) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0009672) q[1];
sx q[1];
rz(-2.3552309) q[1];
sx q[1];
rz(-1.0746075) q[1];
x q[2];
rz(0.23779121) q[3];
sx q[3];
rz(-1.8565408) q[3];
sx q[3];
rz(-1.745949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5932172) q[2];
sx q[2];
rz(-2.2943353) q[2];
sx q[2];
rz(0.66620052) q[2];
rz(3.0001452) q[3];
sx q[3];
rz(-1.5513159) q[3];
sx q[3];
rz(-0.9945873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87404609) q[0];
sx q[0];
rz(-2.5135437) q[0];
sx q[0];
rz(-0.48126599) q[0];
rz(-2.4080243) q[1];
sx q[1];
rz(-1.8736519) q[1];
sx q[1];
rz(-0.092863277) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32321445) q[0];
sx q[0];
rz(-2.0131454) q[0];
sx q[0];
rz(2.0967006) q[0];
rz(-pi) q[1];
rz(-1.5728371) q[2];
sx q[2];
rz(-1.5321595) q[2];
sx q[2];
rz(-2.3734951) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0561547) q[1];
sx q[1];
rz(-0.78821086) q[1];
sx q[1];
rz(-2.5484392) q[1];
x q[2];
rz(1.2081356) q[3];
sx q[3];
rz(-1.3247962) q[3];
sx q[3];
rz(-0.73482692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7218472) q[2];
sx q[2];
rz(-1.6215308) q[2];
sx q[2];
rz(1.9160371) q[2];
rz(-0.078977481) q[3];
sx q[3];
rz(-2.0073399) q[3];
sx q[3];
rz(-1.2823766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6135898) q[0];
sx q[0];
rz(-2.4230175) q[0];
sx q[0];
rz(-0.80793107) q[0];
rz(-1.5241874) q[1];
sx q[1];
rz(-0.15847358) q[1];
sx q[1];
rz(-2.4375516) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2794802) q[0];
sx q[0];
rz(-0.021678019) q[0];
sx q[0];
rz(1.9901748) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1104326) q[2];
sx q[2];
rz(-0.0093447091) q[2];
sx q[2];
rz(1.921794) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0615394) q[1];
sx q[1];
rz(-2.1916917) q[1];
sx q[1];
rz(3.1136334) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3392884) q[3];
sx q[3];
rz(-0.62776126) q[3];
sx q[3];
rz(-0.14428142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0713221) q[2];
sx q[2];
rz(-0.68838781) q[2];
sx q[2];
rz(-2.2570611) q[2];
rz(0.5136579) q[3];
sx q[3];
rz(-1.9653886) q[3];
sx q[3];
rz(0.48532143) q[3];
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
rz(-3.1270444) q[0];
sx q[0];
rz(-2.3764648) q[0];
sx q[0];
rz(0.86084086) q[0];
rz(0.035482081) q[1];
sx q[1];
rz(-1.2018964) q[1];
sx q[1];
rz(-1.6532345) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5527715) q[0];
sx q[0];
rz(-0.5824005) q[0];
sx q[0];
rz(-0.78314535) q[0];
rz(-pi) q[1];
rz(2.9126337) q[2];
sx q[2];
rz(-1.3334074) q[2];
sx q[2];
rz(1.6473532) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1775629) q[1];
sx q[1];
rz(-0.42516758) q[1];
sx q[1];
rz(-2.0079101) q[1];
rz(2.8163696) q[3];
sx q[3];
rz(-2.8517493) q[3];
sx q[3];
rz(-1.0018283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.615768) q[2];
sx q[2];
rz(-2.0066228) q[2];
sx q[2];
rz(1.5110678) q[2];
rz(2.3903971) q[3];
sx q[3];
rz(-0.54371756) q[3];
sx q[3];
rz(2.3403919) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5649331) q[0];
sx q[0];
rz(-1.3929921) q[0];
sx q[0];
rz(-0.15889731) q[0];
rz(-1.502602) q[1];
sx q[1];
rz(-1.331012) q[1];
sx q[1];
rz(-1.9749036) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5432271) q[0];
sx q[0];
rz(-2.1413714) q[0];
sx q[0];
rz(1.8680509) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8966271) q[2];
sx q[2];
rz(-1.3052757) q[2];
sx q[2];
rz(-0.42761432) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.37196463) q[1];
sx q[1];
rz(-1.5329016) q[1];
sx q[1];
rz(-0.016788646) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1205584) q[3];
sx q[3];
rz(-1.5411177) q[3];
sx q[3];
rz(-0.51845394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4022973) q[2];
sx q[2];
rz(-0.6604971) q[2];
sx q[2];
rz(0.26563773) q[2];
rz(-1.9499251) q[3];
sx q[3];
rz(-1.8040413) q[3];
sx q[3];
rz(-0.55408293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3091076) q[0];
sx q[0];
rz(-0.6159679) q[0];
sx q[0];
rz(-0.68514222) q[0];
rz(2.5849672) q[1];
sx q[1];
rz(-1.6770505) q[1];
sx q[1];
rz(-2.305078) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5852528) q[0];
sx q[0];
rz(-1.9808021) q[0];
sx q[0];
rz(-2.7707731) q[0];
rz(-pi) q[1];
rz(-0.30843251) q[2];
sx q[2];
rz(-2.0698202) q[2];
sx q[2];
rz(-0.67828808) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.592258) q[1];
sx q[1];
rz(-0.95227949) q[1];
sx q[1];
rz(1.0299512) q[1];
rz(2.9850858) q[3];
sx q[3];
rz(-1.8409074) q[3];
sx q[3];
rz(2.6685985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4461925) q[2];
sx q[2];
rz(-1.7170898) q[2];
sx q[2];
rz(-2.2851473) q[2];
rz(0.19927464) q[3];
sx q[3];
rz(-0.98237413) q[3];
sx q[3];
rz(-0.91163951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0675426) q[0];
sx q[0];
rz(-2.1265125) q[0];
sx q[0];
rz(1.4763005) q[0];
rz(2.603727) q[1];
sx q[1];
rz(-1.8712416) q[1];
sx q[1];
rz(1.3457294) q[1];
rz(-2.4171038) q[2];
sx q[2];
rz(-2.2064887) q[2];
sx q[2];
rz(0.51842368) q[2];
rz(0.84705254) q[3];
sx q[3];
rz(-1.3624227) q[3];
sx q[3];
rz(-0.25400454) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
