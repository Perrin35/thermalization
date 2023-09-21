OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.37378398) q[0];
sx q[0];
rz(-2.7019579) q[0];
sx q[0];
rz(0.08131942) q[0];
rz(0.65027872) q[1];
sx q[1];
rz(-1.283409) q[1];
sx q[1];
rz(-2.3587956) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62277943) q[0];
sx q[0];
rz(-0.25758994) q[0];
sx q[0];
rz(1.8151087) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15687234) q[2];
sx q[2];
rz(-1.1638767) q[2];
sx q[2];
rz(3.0770609) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1374955) q[1];
sx q[1];
rz(-1.4427408) q[1];
sx q[1];
rz(1.1346243) q[1];
rz(-2.8257915) q[3];
sx q[3];
rz(-1.6557906) q[3];
sx q[3];
rz(1.0335361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89447442) q[2];
sx q[2];
rz(-2.1366182) q[2];
sx q[2];
rz(-3.0257814) q[2];
rz(1.5995021) q[3];
sx q[3];
rz(-0.096297979) q[3];
sx q[3];
rz(1.0533062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88749921) q[0];
sx q[0];
rz(-0.54953456) q[0];
sx q[0];
rz(-0.19533531) q[0];
rz(-0.37503606) q[1];
sx q[1];
rz(-1.6655567) q[1];
sx q[1];
rz(-2.9017752) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1447906) q[0];
sx q[0];
rz(-1.4142493) q[0];
sx q[0];
rz(-0.24897225) q[0];
rz(-0.96599483) q[2];
sx q[2];
rz(-2.7808393) q[2];
sx q[2];
rz(-2.989907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5259243) q[1];
sx q[1];
rz(-2.3500372) q[1];
sx q[1];
rz(0.073992373) q[1];
rz(-1.5870985) q[3];
sx q[3];
rz(-1.5566412) q[3];
sx q[3];
rz(-2.6027038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5043162) q[2];
sx q[2];
rz(-0.80233032) q[2];
sx q[2];
rz(1.8117388) q[2];
rz(-1.7999533) q[3];
sx q[3];
rz(-1.4989217) q[3];
sx q[3];
rz(1.3247066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2770237) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(-2.4734316) q[0];
rz(1.6502624) q[1];
sx q[1];
rz(-0.69258339) q[1];
sx q[1];
rz(1.0659165) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6324368) q[0];
sx q[0];
rz(-1.8189438) q[0];
sx q[0];
rz(-3.0993673) q[0];
x q[1];
rz(-1.714528) q[2];
sx q[2];
rz(-2.1215237) q[2];
sx q[2];
rz(1.3155754) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6560116) q[1];
sx q[1];
rz(-2.1335568) q[1];
sx q[1];
rz(-2.807711) q[1];
rz(-3.021574) q[3];
sx q[3];
rz(-2.2285322) q[3];
sx q[3];
rz(-1.1413871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9540017) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(-3.1090453) q[2];
rz(0.35999808) q[3];
sx q[3];
rz(-2.0149752) q[3];
sx q[3];
rz(-0.79157296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86768326) q[0];
sx q[0];
rz(-1.5861347) q[0];
sx q[0];
rz(-2.4348863) q[0];
rz(1.2061521) q[1];
sx q[1];
rz(-0.34148347) q[1];
sx q[1];
rz(1.6548086) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1409722) q[0];
sx q[0];
rz(-1.2396493) q[0];
sx q[0];
rz(-1.2739869) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47866486) q[2];
sx q[2];
rz(-2.475111) q[2];
sx q[2];
rz(2.037231) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0123803) q[1];
sx q[1];
rz(-2.817569) q[1];
sx q[1];
rz(-2.654241) q[1];
rz(-pi) q[2];
rz(-2.3151822) q[3];
sx q[3];
rz(-1.5658169) q[3];
sx q[3];
rz(1.6299562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5450181) q[2];
sx q[2];
rz(-0.89349616) q[2];
sx q[2];
rz(-1.0085227) q[2];
rz(-2.0452943) q[3];
sx q[3];
rz(-1.9129646) q[3];
sx q[3];
rz(3.0116459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0356045) q[0];
sx q[0];
rz(-0.85177079) q[0];
sx q[0];
rz(1.6812356) q[0];
rz(1.5885072) q[1];
sx q[1];
rz(-1.9057143) q[1];
sx q[1];
rz(-3.1255186) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92656266) q[0];
sx q[0];
rz(-0.83850551) q[0];
sx q[0];
rz(2.3123884) q[0];
rz(-0.0083382567) q[2];
sx q[2];
rz(-1.1401046) q[2];
sx q[2];
rz(-0.24239937) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2542418) q[1];
sx q[1];
rz(-2.4130531) q[1];
sx q[1];
rz(0.90708797) q[1];
rz(0.36230476) q[3];
sx q[3];
rz(-1.575003) q[3];
sx q[3];
rz(-1.6844974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1092704) q[2];
sx q[2];
rz(-1.0291928) q[2];
sx q[2];
rz(-0.62189046) q[2];
rz(-2.0444929) q[3];
sx q[3];
rz(-2.3648839) q[3];
sx q[3];
rz(-2.9212852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9428403) q[0];
sx q[0];
rz(-0.0033012882) q[0];
sx q[0];
rz(-0.90676701) q[0];
rz(-2.3268907) q[1];
sx q[1];
rz(-0.68836132) q[1];
sx q[1];
rz(-1.9168568) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46366102) q[0];
sx q[0];
rz(-1.3555129) q[0];
sx q[0];
rz(-2.2023375) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0133063) q[2];
sx q[2];
rz(-1.3820717) q[2];
sx q[2];
rz(2.92885) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.743222) q[1];
sx q[1];
rz(-2.1757158) q[1];
sx q[1];
rz(0.066992316) q[1];
rz(1.6812427) q[3];
sx q[3];
rz(-1.1718281) q[3];
sx q[3];
rz(-0.16050592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1427052) q[2];
sx q[2];
rz(-2.3110516) q[2];
sx q[2];
rz(-2.2018946) q[2];
rz(0.21329221) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(-1.3903769) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9796824) q[0];
sx q[0];
rz(-0.96452159) q[0];
sx q[0];
rz(0.58037037) q[0];
rz(-1.0549818) q[1];
sx q[1];
rz(-1.4529198) q[1];
sx q[1];
rz(-2.4408128) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0133007) q[0];
sx q[0];
rz(-1.4585146) q[0];
sx q[0];
rz(-2.0237472) q[0];
rz(-pi) q[1];
rz(-1.4540265) q[2];
sx q[2];
rz(-1.8883369) q[2];
sx q[2];
rz(1.5955398) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.62774819) q[1];
sx q[1];
rz(-2.5690837) q[1];
sx q[1];
rz(-1.4036914) q[1];
rz(2.3015162) q[3];
sx q[3];
rz(-1.5080161) q[3];
sx q[3];
rz(-1.9907794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7769527) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(3.11943) q[2];
rz(-0.74470216) q[3];
sx q[3];
rz(-2.0245168) q[3];
sx q[3];
rz(-2.7409592) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3547524) q[0];
sx q[0];
rz(-2.1180034) q[0];
sx q[0];
rz(1.4165075) q[0];
rz(1.3757061) q[1];
sx q[1];
rz(-1.7388652) q[1];
sx q[1];
rz(0.8917121) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69670024) q[0];
sx q[0];
rz(-0.80388821) q[0];
sx q[0];
rz(-0.21960396) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5614369) q[2];
sx q[2];
rz(-1.9556502) q[2];
sx q[2];
rz(-0.58154026) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2878694) q[1];
sx q[1];
rz(-2.3239377) q[1];
sx q[1];
rz(-0.57662782) q[1];
x q[2];
rz(0.91120054) q[3];
sx q[3];
rz(-1.824114) q[3];
sx q[3];
rz(2.5384056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4884168) q[2];
sx q[2];
rz(-1.9501053) q[2];
sx q[2];
rz(0.78424224) q[2];
rz(-2.6358321) q[3];
sx q[3];
rz(-0.85251802) q[3];
sx q[3];
rz(0.23323664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6655675) q[0];
sx q[0];
rz(-1.5371756) q[0];
sx q[0];
rz(0.72203565) q[0];
rz(-0.33323914) q[1];
sx q[1];
rz(-1.9457341) q[1];
sx q[1];
rz(-1.3649712) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46699539) q[0];
sx q[0];
rz(-0.82775263) q[0];
sx q[0];
rz(-1.6270431) q[0];
rz(-pi) q[1];
rz(-3.0579371) q[2];
sx q[2];
rz(-1.9803279) q[2];
sx q[2];
rz(1.4510029) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2024467) q[1];
sx q[1];
rz(-2.2447963) q[1];
sx q[1];
rz(2.0161611) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2737208) q[3];
sx q[3];
rz(-1.7617412) q[3];
sx q[3];
rz(1.8553268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0329131) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(1.2314679) q[2];
rz(3.1075297) q[3];
sx q[3];
rz(-1.27682) q[3];
sx q[3];
rz(-0.6347707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.05474) q[0];
sx q[0];
rz(-0.56607902) q[0];
sx q[0];
rz(1.4779133) q[0];
rz(1.0832896) q[1];
sx q[1];
rz(-1.7419086) q[1];
sx q[1];
rz(0.96819425) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9111239) q[0];
sx q[0];
rz(-1.9096806) q[0];
sx q[0];
rz(1.1662657) q[0];
rz(-pi) q[1];
rz(0.88629006) q[2];
sx q[2];
rz(-2.5286525) q[2];
sx q[2];
rz(0.71000368) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72659513) q[1];
sx q[1];
rz(-1.8907049) q[1];
sx q[1];
rz(1.725561) q[1];
rz(-pi) q[2];
rz(3.1413583) q[3];
sx q[3];
rz(-1.8929385) q[3];
sx q[3];
rz(-2.5664267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0633462) q[2];
sx q[2];
rz(-2.4137256) q[2];
sx q[2];
rz(-3.1402804) q[2];
rz(-2.0007658) q[3];
sx q[3];
rz(-1.8246548) q[3];
sx q[3];
rz(1.3425945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7447727) q[0];
sx q[0];
rz(-2.0347432) q[0];
sx q[0];
rz(1.9532935) q[0];
rz(2.7753579) q[1];
sx q[1];
rz(-1.2013422) q[1];
sx q[1];
rz(1.3399301) q[1];
rz(1.2583854) q[2];
sx q[2];
rz(-2.2956143) q[2];
sx q[2];
rz(2.14239) q[2];
rz(-1.0352186) q[3];
sx q[3];
rz(-0.62789161) q[3];
sx q[3];
rz(-0.37554489) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];