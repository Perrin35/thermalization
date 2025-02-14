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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2446063) q[0];
sx q[0];
rz(-1.7660487) q[0];
sx q[0];
rz(0.050873916) q[0];
rz(-pi) q[1];
rz(1.6221912) q[2];
sx q[2];
rz(-1.2631466) q[2];
sx q[2];
rz(2.4884698) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.611022) q[1];
sx q[1];
rz(-2.321302) q[1];
sx q[1];
rz(1.4613226) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4769931) q[3];
sx q[3];
rz(-2.0076111) q[3];
sx q[3];
rz(2.7561881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.92496282) q[2];
sx q[2];
rz(-1.7731885) q[2];
sx q[2];
rz(-1.487757) q[2];
rz(-1.9803068) q[3];
sx q[3];
rz(-1.0043283) q[3];
sx q[3];
rz(-1.0950834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0682812) q[0];
sx q[0];
rz(-0.40638766) q[0];
sx q[0];
rz(-0.43346369) q[0];
rz(-2.0794226) q[1];
sx q[1];
rz(-1.5824317) q[1];
sx q[1];
rz(0.72057048) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90209157) q[0];
sx q[0];
rz(-1.4441617) q[0];
sx q[0];
rz(0.037221639) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24777221) q[2];
sx q[2];
rz(-1.0508014) q[2];
sx q[2];
rz(1.2733851) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6241265) q[1];
sx q[1];
rz(-1.7192914) q[1];
sx q[1];
rz(-2.7854169) q[1];
rz(-pi) q[2];
rz(-1.0834915) q[3];
sx q[3];
rz(-1.8326899) q[3];
sx q[3];
rz(-2.0341121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7945107) q[2];
sx q[2];
rz(-1.2257267) q[2];
sx q[2];
rz(-2.1573055) q[2];
rz(-0.85683626) q[3];
sx q[3];
rz(-0.58647668) q[3];
sx q[3];
rz(-1.2543359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03265753) q[0];
sx q[0];
rz(-0.28626838) q[0];
sx q[0];
rz(-0.74288595) q[0];
rz(-3.1318829) q[1];
sx q[1];
rz(-1.5747036) q[1];
sx q[1];
rz(2.8579393) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8359287) q[0];
sx q[0];
rz(-1.2244512) q[0];
sx q[0];
rz(0.54327528) q[0];
rz(-pi) q[1];
rz(-2.9417324) q[2];
sx q[2];
rz(-1.4132573) q[2];
sx q[2];
rz(3.1197366) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3481118) q[1];
sx q[1];
rz(-1.2657796) q[1];
sx q[1];
rz(-0.17615894) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.795122) q[3];
sx q[3];
rz(-3.0375518) q[3];
sx q[3];
rz(1.4416172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0297086) q[2];
sx q[2];
rz(-1.6969029) q[2];
sx q[2];
rz(2.4857944) q[2];
rz(0.2119952) q[3];
sx q[3];
rz(-2.5097804) q[3];
sx q[3];
rz(-1.1727715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1126605) q[0];
sx q[0];
rz(-0.49429587) q[0];
sx q[0];
rz(1.5546881) q[0];
rz(2.9117865) q[1];
sx q[1];
rz(-2.4202012) q[1];
sx q[1];
rz(1.302964) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36520805) q[0];
sx q[0];
rz(-1.0321277) q[0];
sx q[0];
rz(-1.046565) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26891687) q[2];
sx q[2];
rz(-1.8559716) q[2];
sx q[2];
rz(-2.0320924) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.069217056) q[1];
sx q[1];
rz(-1.0472968) q[1];
sx q[1];
rz(1.5956466) q[1];
rz(-pi) q[2];
rz(-2.0992364) q[3];
sx q[3];
rz(-2.5109249) q[3];
sx q[3];
rz(1.809465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.0064696781) q[2];
sx q[2];
rz(-0.80265704) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2876494) q[0];
sx q[0];
rz(-2.8048153) q[0];
sx q[0];
rz(2.7119998) q[0];
rz(-2.8751539) q[1];
sx q[1];
rz(-0.8100422) q[1];
sx q[1];
rz(-1.0688759) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7212413) q[0];
sx q[0];
rz(-1.9801894) q[0];
sx q[0];
rz(2.4024525) q[0];
rz(-pi) q[1];
rz(-1.8891625) q[2];
sx q[2];
rz(-0.7944383) q[2];
sx q[2];
rz(2.3311262) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4539448) q[1];
sx q[1];
rz(-1.8275211) q[1];
sx q[1];
rz(-2.686108) q[1];
rz(2.3924218) q[3];
sx q[3];
rz(-0.35532829) q[3];
sx q[3];
rz(2.4859249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8916919) q[2];
sx q[2];
rz(-0.72600681) q[2];
sx q[2];
rz(-0.1740087) q[2];
rz(1.8788667) q[3];
sx q[3];
rz(-1.5401253) q[3];
sx q[3];
rz(1.5866097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1332909) q[0];
sx q[0];
rz(-1.0187848) q[0];
sx q[0];
rz(-0.89023501) q[0];
rz(-1.7091735) q[1];
sx q[1];
rz(-2.1232257) q[1];
sx q[1];
rz(-0.11996809) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9667119) q[0];
sx q[0];
rz(-0.57583664) q[0];
sx q[0];
rz(2.8737646) q[0];
rz(0.89996296) q[2];
sx q[2];
rz(-2.6113178) q[2];
sx q[2];
rz(1.8535623) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6585582) q[1];
sx q[1];
rz(-0.93340988) q[1];
sx q[1];
rz(-1.263522) q[1];
x q[2];
rz(-2.8268485) q[3];
sx q[3];
rz(-2.5786244) q[3];
sx q[3];
rz(-1.986435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.031781901) q[2];
sx q[2];
rz(-1.6620771) q[2];
sx q[2];
rz(0.39153448) q[2];
rz(-1.7560962) q[3];
sx q[3];
rz(-2.9031495) q[3];
sx q[3];
rz(-2.5029206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93381298) q[0];
sx q[0];
rz(-2.5901828) q[0];
sx q[0];
rz(-1.3054003) q[0];
rz(0.50453672) q[1];
sx q[1];
rz(-1.4089818) q[1];
sx q[1];
rz(-2.9171468) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69687122) q[0];
sx q[0];
rz(-1.5550647) q[0];
sx q[0];
rz(-3.12495) q[0];
x q[1];
rz(2.6487892) q[2];
sx q[2];
rz(-1.3818008) q[2];
sx q[2];
rz(2.1616621) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.45221165) q[1];
sx q[1];
rz(-1.895001) q[1];
sx q[1];
rz(-0.28181847) q[1];
rz(-pi) q[2];
rz(-0.32308635) q[3];
sx q[3];
rz(-1.2189019) q[3];
sx q[3];
rz(0.066095933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18181431) q[2];
sx q[2];
rz(-1.10428) q[2];
sx q[2];
rz(-0.9474729) q[2];
rz(3.0090561) q[3];
sx q[3];
rz(-2.4819076) q[3];
sx q[3];
rz(-0.98602492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9082311) q[0];
sx q[0];
rz(-2.7955604) q[0];
sx q[0];
rz(-2.76873) q[0];
rz(0.33509675) q[1];
sx q[1];
rz(-1.9357888) q[1];
sx q[1];
rz(2.0705409) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.459254) q[0];
sx q[0];
rz(-1.5821348) q[0];
sx q[0];
rz(3.0724694) q[0];
rz(2.5899289) q[2];
sx q[2];
rz(-1.3886189) q[2];
sx q[2];
rz(-2.7523486) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1879857) q[1];
sx q[1];
rz(-1.803557) q[1];
sx q[1];
rz(0.25165073) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.095094918) q[3];
sx q[3];
rz(-1.8912867) q[3];
sx q[3];
rz(1.2456428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4272473) q[2];
sx q[2];
rz(-1.5421966) q[2];
sx q[2];
rz(1.0672807) q[2];
rz(-2.769477) q[3];
sx q[3];
rz(-0.99388638) q[3];
sx q[3];
rz(-2.9414162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9099971) q[0];
sx q[0];
rz(-2.5439926) q[0];
sx q[0];
rz(-1.4071314) q[0];
rz(2.097997) q[1];
sx q[1];
rz(-2.2035972) q[1];
sx q[1];
rz(-2.6503906) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.585116) q[0];
sx q[0];
rz(-2.399754) q[0];
sx q[0];
rz(1.1354574) q[0];
rz(0.16881659) q[2];
sx q[2];
rz(-1.4028609) q[2];
sx q[2];
rz(2.4562648) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.32090998) q[1];
sx q[1];
rz(-1.9100128) q[1];
sx q[1];
rz(-1.7103819) q[1];
rz(-pi) q[2];
rz(2.7094531) q[3];
sx q[3];
rz(-2.7305805) q[3];
sx q[3];
rz(2.7087045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4444943) q[2];
sx q[2];
rz(-1.5105057) q[2];
sx q[2];
rz(1.0443045) q[2];
rz(2.3704884) q[3];
sx q[3];
rz(-1.6090569) q[3];
sx q[3];
rz(-0.3573629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20030178) q[0];
sx q[0];
rz(-1.8481978) q[0];
sx q[0];
rz(1.6218761) q[0];
rz(1.3007851) q[1];
sx q[1];
rz(-1.3011353) q[1];
sx q[1];
rz(-2.4470952) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25481564) q[0];
sx q[0];
rz(-0.22829994) q[0];
sx q[0];
rz(1.5103025) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7037665) q[2];
sx q[2];
rz(-1.0501692) q[2];
sx q[2];
rz(0.43308738) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.84811454) q[1];
sx q[1];
rz(-1.9047996) q[1];
sx q[1];
rz(3.1298742) q[1];
rz(-2.7312134) q[3];
sx q[3];
rz(-2.346092) q[3];
sx q[3];
rz(0.26887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0684315) q[2];
sx q[2];
rz(-1.7794926) q[2];
sx q[2];
rz(0.51363242) q[2];
rz(0.3178151) q[3];
sx q[3];
rz(-1.5376667) q[3];
sx q[3];
rz(1.0191127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0014872) q[0];
sx q[0];
rz(-1.6248063) q[0];
sx q[0];
rz(-0.90909062) q[0];
rz(-1.2263251) q[1];
sx q[1];
rz(-1.6188123) q[1];
sx q[1];
rz(0.34368044) q[1];
rz(2.4438159) q[2];
sx q[2];
rz(-1.991376) q[2];
sx q[2];
rz(-2.1817653) q[2];
rz(-2.0967284) q[3];
sx q[3];
rz(-1.5726225) q[3];
sx q[3];
rz(-2.8468688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
