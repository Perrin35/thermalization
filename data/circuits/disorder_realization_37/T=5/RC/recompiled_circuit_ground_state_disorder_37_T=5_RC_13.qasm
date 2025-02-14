OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6310298) q[0];
sx q[0];
rz(3.6462311) q[0];
sx q[0];
rz(12.991821) q[0];
rz(-2.5692441) q[1];
sx q[1];
rz(-1.0971789) q[1];
sx q[1];
rz(1.9414577) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5050738) q[0];
sx q[0];
rz(-2.8905792) q[0];
sx q[0];
rz(0.039029718) q[0];
x q[1];
rz(-0.27165551) q[2];
sx q[2];
rz(-1.9052398) q[2];
sx q[2];
rz(-1.3525427) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.81600707) q[1];
sx q[1];
rz(-0.43445265) q[1];
sx q[1];
rz(-0.52838051) q[1];
rz(-pi) q[2];
rz(-2.4704504) q[3];
sx q[3];
rz(-1.6832388) q[3];
sx q[3];
rz(-2.6813316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9878865) q[2];
sx q[2];
rz(-1.1517297) q[2];
sx q[2];
rz(0.64725867) q[2];
rz(2.9636532) q[3];
sx q[3];
rz(-1.8837594) q[3];
sx q[3];
rz(-1.9378763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4030289) q[0];
sx q[0];
rz(-0.084349923) q[0];
sx q[0];
rz(2.6179598) q[0];
rz(-1.2298443) q[1];
sx q[1];
rz(-2.3975027) q[1];
sx q[1];
rz(2.8935166) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1988293) q[0];
sx q[0];
rz(-2.307245) q[0];
sx q[0];
rz(1.7146066) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43681983) q[2];
sx q[2];
rz(-1.9736276) q[2];
sx q[2];
rz(-2.6879329) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4011127) q[1];
sx q[1];
rz(-1.4251627) q[1];
sx q[1];
rz(0.039197368) q[1];
rz(-2.5026191) q[3];
sx q[3];
rz(-0.77497122) q[3];
sx q[3];
rz(2.4955622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6123885) q[2];
sx q[2];
rz(-2.2233621) q[2];
sx q[2];
rz(-2.6017792) q[2];
rz(-0.16048935) q[3];
sx q[3];
rz(-1.548998) q[3];
sx q[3];
rz(-0.81015712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4983343) q[0];
sx q[0];
rz(-0.67688268) q[0];
sx q[0];
rz(0.13370378) q[0];
rz(1.5191822) q[1];
sx q[1];
rz(-1.8459873) q[1];
sx q[1];
rz(-2.349283) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80817079) q[0];
sx q[0];
rz(-1.9229888) q[0];
sx q[0];
rz(2.9461224) q[0];
x q[1];
rz(1.7025204) q[2];
sx q[2];
rz(-1.4274397) q[2];
sx q[2];
rz(-1.6241339) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.49323248) q[1];
sx q[1];
rz(-0.91714749) q[1];
sx q[1];
rz(-0.75809376) q[1];
rz(-pi) q[2];
rz(-0.30012975) q[3];
sx q[3];
rz(-1.5709166) q[3];
sx q[3];
rz(2.4693268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2866659) q[2];
sx q[2];
rz(-2.4714405) q[2];
sx q[2];
rz(-0.83038846) q[2];
rz(1.7928127) q[3];
sx q[3];
rz(-2.1068137) q[3];
sx q[3];
rz(2.5428298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26375672) q[0];
sx q[0];
rz(-2.1649375) q[0];
sx q[0];
rz(2.7184955) q[0];
rz(1.9844203) q[1];
sx q[1];
rz(-1.4856228) q[1];
sx q[1];
rz(2.5194936) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3836648) q[0];
sx q[0];
rz(-2.1970941) q[0];
sx q[0];
rz(2.6561198) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1272264) q[2];
sx q[2];
rz(-2.6831919) q[2];
sx q[2];
rz(2.1200402) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2492952) q[1];
sx q[1];
rz(-1.9782606) q[1];
sx q[1];
rz(-3.1235126) q[1];
x q[2];
rz(3.1142523) q[3];
sx q[3];
rz(-2.2961444) q[3];
sx q[3];
rz(3.0375061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.28079924) q[2];
sx q[2];
rz(-1.7223027) q[2];
sx q[2];
rz(-0.61817509) q[2];
rz(-1.482796) q[3];
sx q[3];
rz(-1.7890472) q[3];
sx q[3];
rz(-1.6331204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16267714) q[0];
sx q[0];
rz(-1.3293043) q[0];
sx q[0];
rz(2.3090114) q[0];
rz(0.32514462) q[1];
sx q[1];
rz(-2.145576) q[1];
sx q[1];
rz(0.064373374) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8125809) q[0];
sx q[0];
rz(-1.4509177) q[0];
sx q[0];
rz(-2.8146312) q[0];
rz(0.60026786) q[2];
sx q[2];
rz(-1.1986599) q[2];
sx q[2];
rz(2.0780217) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.87582237) q[1];
sx q[1];
rz(-2.0811305) q[1];
sx q[1];
rz(-0.93632621) q[1];
x q[2];
rz(-2.0470555) q[3];
sx q[3];
rz(-2.4570997) q[3];
sx q[3];
rz(1.5409926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9084106) q[2];
sx q[2];
rz(-0.57047129) q[2];
sx q[2];
rz(1.849966) q[2];
rz(-0.49606797) q[3];
sx q[3];
rz(-0.78947624) q[3];
sx q[3];
rz(1.9991649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18601501) q[0];
sx q[0];
rz(-0.29009524) q[0];
sx q[0];
rz(0.41906038) q[0];
rz(2.8298607) q[1];
sx q[1];
rz(-1.5429976) q[1];
sx q[1];
rz(-2.9980803) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79663178) q[0];
sx q[0];
rz(-1.2141742) q[0];
sx q[0];
rz(0.48194569) q[0];
rz(-2.8673929) q[2];
sx q[2];
rz(-1.7975217) q[2];
sx q[2];
rz(-2.271914) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1309424) q[1];
sx q[1];
rz(-0.82320628) q[1];
sx q[1];
rz(0.44013104) q[1];
x q[2];
rz(-2.6426598) q[3];
sx q[3];
rz(-0.29424516) q[3];
sx q[3];
rz(0.1732451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.36025563) q[2];
sx q[2];
rz(-2.2402253) q[2];
sx q[2];
rz(2.0085013) q[2];
rz(-2.0356778) q[3];
sx q[3];
rz(-1.7224885) q[3];
sx q[3];
rz(2.6667986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.356242) q[0];
sx q[0];
rz(-0.79371912) q[0];
sx q[0];
rz(1.445795) q[0];
rz(0.71714199) q[1];
sx q[1];
rz(-1.8627867) q[1];
sx q[1];
rz(-0.76146567) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1061072) q[0];
sx q[0];
rz(-0.64752827) q[0];
sx q[0];
rz(-2.8962563) q[0];
rz(-pi) q[1];
rz(2.285073) q[2];
sx q[2];
rz(-1.1978784) q[2];
sx q[2];
rz(-0.97717092) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.95132212) q[1];
sx q[1];
rz(-0.99748625) q[1];
sx q[1];
rz(-3.0828397) q[1];
x q[2];
rz(-1.3827356) q[3];
sx q[3];
rz(-1.4131964) q[3];
sx q[3];
rz(3.0732875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7722499) q[2];
sx q[2];
rz(-0.49392924) q[2];
sx q[2];
rz(2.210145) q[2];
rz(-2.5233968) q[3];
sx q[3];
rz(-1.5074916) q[3];
sx q[3];
rz(-0.96562323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6922927) q[0];
sx q[0];
rz(-1.7789142) q[0];
sx q[0];
rz(-1.8071254) q[0];
rz(-0.5131228) q[1];
sx q[1];
rz(-0.98315364) q[1];
sx q[1];
rz(-1.9122874) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2199101) q[0];
sx q[0];
rz(-1.1510885) q[0];
sx q[0];
rz(2.3786663) q[0];
rz(0.092809794) q[2];
sx q[2];
rz(-1.3121737) q[2];
sx q[2];
rz(-1.5802204) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9521515) q[1];
sx q[1];
rz(-1.5501889) q[1];
sx q[1];
rz(0.15592305) q[1];
rz(-pi) q[2];
x q[2];
rz(3.119) q[3];
sx q[3];
rz(-0.66231662) q[3];
sx q[3];
rz(-0.8099156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8339771) q[2];
sx q[2];
rz(-0.5270842) q[2];
sx q[2];
rz(2.9311467) q[2];
rz(0.065841913) q[3];
sx q[3];
rz(-0.47271553) q[3];
sx q[3];
rz(-0.79894799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.598572) q[0];
sx q[0];
rz(-0.6544756) q[0];
sx q[0];
rz(1.8684813) q[0];
rz(-1.2741362) q[1];
sx q[1];
rz(-0.68398634) q[1];
sx q[1];
rz(-0.88409105) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7403075) q[0];
sx q[0];
rz(-1.570556) q[0];
sx q[0];
rz(2.1876641) q[0];
rz(-pi) q[1];
rz(-0.52427788) q[2];
sx q[2];
rz(-1.345621) q[2];
sx q[2];
rz(-0.88722992) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6581304) q[1];
sx q[1];
rz(-2.6331055) q[1];
sx q[1];
rz(2.1668651) q[1];
rz(-pi) q[2];
rz(0.69334778) q[3];
sx q[3];
rz(-2.7487488) q[3];
sx q[3];
rz(1.0885914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1511496) q[2];
sx q[2];
rz(-2.1572025) q[2];
sx q[2];
rz(-1.3746064) q[2];
rz(-0.19045842) q[3];
sx q[3];
rz(-2.3951267) q[3];
sx q[3];
rz(1.1838574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9780438) q[0];
sx q[0];
rz(-1.6885641) q[0];
sx q[0];
rz(1.8769886) q[0];
rz(3.0679852) q[1];
sx q[1];
rz(-1.0818447) q[1];
sx q[1];
rz(0.61990613) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093599565) q[0];
sx q[0];
rz(-1.7353829) q[0];
sx q[0];
rz(-0.82407804) q[0];
rz(-pi) q[1];
rz(2.4700353) q[2];
sx q[2];
rz(-2.4118556) q[2];
sx q[2];
rz(0.31873733) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.82295361) q[1];
sx q[1];
rz(-1.3711437) q[1];
sx q[1];
rz(2.9266788) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4510433) q[3];
sx q[3];
rz(-2.6313734) q[3];
sx q[3];
rz(0.0039781489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8668883) q[2];
sx q[2];
rz(-2.4394749) q[2];
sx q[2];
rz(0.14599027) q[2];
rz(2.919096) q[3];
sx q[3];
rz(-0.22495088) q[3];
sx q[3];
rz(2.4116624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9324026) q[0];
sx q[0];
rz(-1.1840191) q[0];
sx q[0];
rz(-2.9970072) q[0];
rz(-1.2296386) q[1];
sx q[1];
rz(-1.3508136) q[1];
sx q[1];
rz(-1.2523686) q[1];
rz(0.31730439) q[2];
sx q[2];
rz(-2.5226421) q[2];
sx q[2];
rz(0.68963827) q[2];
rz(-3.1093507) q[3];
sx q[3];
rz(-1.7809636) q[3];
sx q[3];
rz(-1.8973593) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
