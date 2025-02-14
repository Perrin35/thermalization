OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2151467) q[0];
sx q[0];
rz(0.67987052) q[0];
sx q[0];
rz(9.3240919) q[0];
rz(-2.6868532) q[1];
sx q[1];
rz(-2.4603031) q[1];
sx q[1];
rz(1.0256306) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5012623) q[0];
sx q[0];
rz(-0.20168951) q[0];
sx q[0];
rz(1.8224688) q[0];
rz(1.6221912) q[2];
sx q[2];
rz(-1.8784461) q[2];
sx q[2];
rz(-2.4884698) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.611022) q[1];
sx q[1];
rz(-2.321302) q[1];
sx q[1];
rz(1.6802701) q[1];
x q[2];
rz(2.4769931) q[3];
sx q[3];
rz(-2.0076111) q[3];
sx q[3];
rz(-2.7561881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.92496282) q[2];
sx q[2];
rz(-1.7731885) q[2];
sx q[2];
rz(1.6538357) q[2];
rz(1.9803068) q[3];
sx q[3];
rz(-1.0043283) q[3];
sx q[3];
rz(-2.0465093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0733114) q[0];
sx q[0];
rz(-2.735205) q[0];
sx q[0];
rz(-0.43346369) q[0];
rz(2.0794226) q[1];
sx q[1];
rz(-1.5824317) q[1];
sx q[1];
rz(-0.72057048) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2395011) q[0];
sx q[0];
rz(-1.6974309) q[0];
sx q[0];
rz(-3.104371) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8938204) q[2];
sx q[2];
rz(-1.0508014) q[2];
sx q[2];
rz(-1.8682075) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.51746619) q[1];
sx q[1];
rz(-1.7192914) q[1];
sx q[1];
rz(-0.35617574) q[1];
rz(-pi) q[2];
rz(-1.0834915) q[3];
sx q[3];
rz(-1.3089027) q[3];
sx q[3];
rz(-1.1074806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.34708193) q[2];
sx q[2];
rz(-1.9158659) q[2];
sx q[2];
rz(0.98428717) q[2];
rz(-2.2847564) q[3];
sx q[3];
rz(-0.58647668) q[3];
sx q[3];
rz(1.2543359) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1089351) q[0];
sx q[0];
rz(-0.28626838) q[0];
sx q[0];
rz(-0.74288595) q[0];
rz(-0.0097097857) q[1];
sx q[1];
rz(-1.5747036) q[1];
sx q[1];
rz(0.28365338) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8946202) q[0];
sx q[0];
rz(-0.63480154) q[0];
sx q[0];
rz(-0.60946861) q[0];
rz(-pi) q[1];
rz(2.9417324) q[2];
sx q[2];
rz(-1.4132573) q[2];
sx q[2];
rz(-3.1197366) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3481118) q[1];
sx q[1];
rz(-1.2657796) q[1];
sx q[1];
rz(-0.17615894) q[1];
rz(-pi) q[2];
x q[2];
rz(1.795122) q[3];
sx q[3];
rz(-3.0375518) q[3];
sx q[3];
rz(1.6999754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11188406) q[2];
sx q[2];
rz(-1.6969029) q[2];
sx q[2];
rz(2.4857944) q[2];
rz(-0.2119952) q[3];
sx q[3];
rz(-2.5097804) q[3];
sx q[3];
rz(-1.9688212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(3.1126605) q[0];
sx q[0];
rz(-0.49429587) q[0];
sx q[0];
rz(1.5546881) q[0];
rz(0.2298062) q[1];
sx q[1];
rz(-0.72139144) q[1];
sx q[1];
rz(-1.8386286) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4939319) q[0];
sx q[0];
rz(-1.126673) q[0];
sx q[0];
rz(2.537389) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26891687) q[2];
sx q[2];
rz(-1.8559716) q[2];
sx q[2];
rz(-2.0320924) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.11889549) q[1];
sx q[1];
rz(-2.6175584) q[1];
sx q[1];
rz(-3.0985718) q[1];
x q[2];
rz(2.7888621) q[3];
sx q[3];
rz(-1.0364854) q[3];
sx q[3];
rz(-1.9580458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.135123) q[2];
sx q[2];
rz(-2.3389356) q[2];
sx q[2];
rz(-0.50372493) q[2];
rz(1.6040365) q[3];
sx q[3];
rz(-1.9683014) q[3];
sx q[3];
rz(1.5776207) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2876494) q[0];
sx q[0];
rz(-0.33677736) q[0];
sx q[0];
rz(-2.7119998) q[0];
rz(0.26643878) q[1];
sx q[1];
rz(-2.3315505) q[1];
sx q[1];
rz(-2.0727167) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80238587) q[0];
sx q[0];
rz(-2.2369719) q[0];
sx q[0];
rz(2.1016913) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2524302) q[2];
sx q[2];
rz(-0.7944383) q[2];
sx q[2];
rz(-0.81046644) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7366745) q[1];
sx q[1];
rz(-2.6231986) q[1];
sx q[1];
rz(-2.6035518) q[1];
x q[2];
rz(2.3924218) q[3];
sx q[3];
rz(-2.7862644) q[3];
sx q[3];
rz(0.65566777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8916919) q[2];
sx q[2];
rz(-0.72600681) q[2];
sx q[2];
rz(-0.1740087) q[2];
rz(1.2627259) q[3];
sx q[3];
rz(-1.6014674) q[3];
sx q[3];
rz(1.5866097) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0083017666) q[0];
sx q[0];
rz(-1.0187848) q[0];
sx q[0];
rz(-2.2513576) q[0];
rz(1.7091735) q[1];
sx q[1];
rz(-2.1232257) q[1];
sx q[1];
rz(0.11996809) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4910866) q[0];
sx q[0];
rz(-1.0179369) q[0];
sx q[0];
rz(-1.4006459) q[0];
rz(-pi) q[1];
rz(0.89996296) q[2];
sx q[2];
rz(-2.6113178) q[2];
sx q[2];
rz(-1.2880304) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.97291291) q[1];
sx q[1];
rz(-2.4434103) q[1];
sx q[1];
rz(2.7538128) q[1];
rz(-pi) q[2];
rz(2.6011068) q[3];
sx q[3];
rz(-1.7367762) q[3];
sx q[3];
rz(-2.9946208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1098108) q[2];
sx q[2];
rz(-1.4795156) q[2];
sx q[2];
rz(-0.39153448) q[2];
rz(-1.3854965) q[3];
sx q[3];
rz(-2.9031495) q[3];
sx q[3];
rz(2.5029206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2077797) q[0];
sx q[0];
rz(-0.55140984) q[0];
sx q[0];
rz(-1.3054003) q[0];
rz(0.50453672) q[1];
sx q[1];
rz(-1.7326109) q[1];
sx q[1];
rz(-0.22444589) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11672606) q[0];
sx q[0];
rz(-3.118692) q[0];
sx q[0];
rz(-2.3842627) q[0];
rz(-pi) q[1];
rz(2.7573757) q[2];
sx q[2];
rz(-0.52500341) q[2];
sx q[2];
rz(0.25431654) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.45221165) q[1];
sx q[1];
rz(-1.2465917) q[1];
sx q[1];
rz(0.28181847) q[1];
rz(-pi) q[2];
rz(2.8185063) q[3];
sx q[3];
rz(-1.9226908) q[3];
sx q[3];
rz(3.0754967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.18181431) q[2];
sx q[2];
rz(-1.10428) q[2];
sx q[2];
rz(-2.1941198) q[2];
rz(-3.0090561) q[3];
sx q[3];
rz(-0.65968502) q[3];
sx q[3];
rz(2.1555677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.2333616) q[0];
sx q[0];
rz(-2.7955604) q[0];
sx q[0];
rz(0.37286266) q[0];
rz(-2.8064959) q[1];
sx q[1];
rz(-1.2058039) q[1];
sx q[1];
rz(1.0710517) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27386943) q[0];
sx q[0];
rz(-3.0715471) q[0];
sx q[0];
rz(-2.9788736) q[0];
rz(-0.55166371) q[2];
sx q[2];
rz(-1.3886189) q[2];
sx q[2];
rz(0.38924402) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0275299) q[1];
sx q[1];
rz(-2.8004871) q[1];
sx q[1];
rz(2.3807664) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.849412) q[3];
sx q[3];
rz(-2.8077586) q[3];
sx q[3];
rz(-0.95164559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4272473) q[2];
sx q[2];
rz(-1.5421966) q[2];
sx q[2];
rz(2.074312) q[2];
rz(2.769477) q[3];
sx q[3];
rz(-0.99388638) q[3];
sx q[3];
rz(2.9414162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2315955) q[0];
sx q[0];
rz(-0.59760004) q[0];
sx q[0];
rz(-1.7344612) q[0];
rz(2.097997) q[1];
sx q[1];
rz(-0.93799543) q[1];
sx q[1];
rz(2.6503906) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1192899) q[0];
sx q[0];
rz(-0.91141846) q[0];
sx q[0];
rz(-0.36880606) q[0];
rz(-pi) q[1];
rz(0.78989002) q[2];
sx q[2];
rz(-0.23755506) q[2];
sx q[2];
rz(-1.4804763) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8206827) q[1];
sx q[1];
rz(-1.9100128) q[1];
sx q[1];
rz(-1.4312107) q[1];
rz(-pi) q[2];
rz(2.7647385) q[3];
sx q[3];
rz(-1.7389193) q[3];
sx q[3];
rz(0.73790077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4444943) q[2];
sx q[2];
rz(-1.5105057) q[2];
sx q[2];
rz(-2.0972882) q[2];
rz(-0.77110428) q[3];
sx q[3];
rz(-1.5325357) q[3];
sx q[3];
rz(0.3573629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20030178) q[0];
sx q[0];
rz(-1.2933949) q[0];
sx q[0];
rz(-1.6218761) q[0];
rz(-1.3007851) q[1];
sx q[1];
rz(-1.8404574) q[1];
sx q[1];
rz(-2.4470952) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3169169) q[0];
sx q[0];
rz(-1.3429214) q[0];
sx q[0];
rz(3.1275463) q[0];
rz(-pi) q[1];
rz(-1.7037665) q[2];
sx q[2];
rz(-2.0914234) q[2];
sx q[2];
rz(-0.43308738) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3292103) q[1];
sx q[1];
rz(-0.33420104) q[1];
sx q[1];
rz(1.5370395) q[1];
x q[2];
rz(2.3894074) q[3];
sx q[3];
rz(-1.859741) q[3];
sx q[3];
rz(-2.1352701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0684315) q[2];
sx q[2];
rz(-1.7794926) q[2];
sx q[2];
rz(-0.51363242) q[2];
rz(-0.3178151) q[3];
sx q[3];
rz(-1.603926) q[3];
sx q[3];
rz(-2.1224799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0014872) q[0];
sx q[0];
rz(-1.6248063) q[0];
sx q[0];
rz(-0.90909062) q[0];
rz(1.9152676) q[1];
sx q[1];
rz(-1.6188123) q[1];
sx q[1];
rz(0.34368044) q[1];
rz(-0.69777674) q[2];
sx q[2];
rz(-1.991376) q[2];
sx q[2];
rz(-2.1817653) q[2];
rz(-1.0448643) q[3];
sx q[3];
rz(-1.5689701) q[3];
sx q[3];
rz(0.29472385) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
