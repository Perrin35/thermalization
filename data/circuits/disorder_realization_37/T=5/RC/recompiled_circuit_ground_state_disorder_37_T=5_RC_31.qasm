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
rz(0.57234859) q[1];
sx q[1];
rz(4.2387716) q[1];
sx q[1];
rz(7.4833202) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6365189) q[0];
sx q[0];
rz(-2.8905792) q[0];
sx q[0];
rz(-0.039029718) q[0];
x q[1];
rz(-2.8699371) q[2];
sx q[2];
rz(-1.9052398) q[2];
sx q[2];
rz(1.3525427) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3878183) q[1];
sx q[1];
rz(-1.1987616) q[1];
sx q[1];
rz(1.8005936) q[1];
rz(-pi) q[2];
rz(-1.7140017) q[3];
sx q[3];
rz(-2.2369336) q[3];
sx q[3];
rz(-1.0216658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1537062) q[2];
sx q[2];
rz(-1.1517297) q[2];
sx q[2];
rz(2.494334) q[2];
rz(-0.17793947) q[3];
sx q[3];
rz(-1.2578332) q[3];
sx q[3];
rz(-1.2037163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7385638) q[0];
sx q[0];
rz(-3.0572427) q[0];
sx q[0];
rz(-0.52363288) q[0];
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
rz(-2.6726674) q[0];
sx q[0];
rz(-1.4644196) q[0];
sx q[0];
rz(0.74161462) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7047728) q[2];
sx q[2];
rz(-1.9736276) q[2];
sx q[2];
rz(2.6879329) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6650508) q[1];
sx q[1];
rz(-0.15078031) q[1];
sx q[1];
rz(-1.3097179) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0421776) q[3];
sx q[3];
rz(-2.1671767) q[3];
sx q[3];
rz(-1.6906052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5292042) q[2];
sx q[2];
rz(-2.2233621) q[2];
sx q[2];
rz(2.6017792) q[2];
rz(0.16048935) q[3];
sx q[3];
rz(-1.548998) q[3];
sx q[3];
rz(-2.3314355) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4983343) q[0];
sx q[0];
rz(-0.67688268) q[0];
sx q[0];
rz(-0.13370378) q[0];
rz(1.5191822) q[1];
sx q[1];
rz(-1.8459873) q[1];
sx q[1];
rz(0.79230961) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2871029) q[0];
sx q[0];
rz(-0.40081319) q[0];
sx q[0];
rz(2.0569747) q[0];
rz(-2.4033264) q[2];
sx q[2];
rz(-2.9472137) q[2];
sx q[2];
rz(0.76972085) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.49323248) q[1];
sx q[1];
rz(-2.2244452) q[1];
sx q[1];
rz(2.3834989) q[1];
x q[2];
rz(-1.5706704) q[3];
sx q[3];
rz(-1.8709261) q[3];
sx q[3];
rz(-2.2430994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85492674) q[2];
sx q[2];
rz(-2.4714405) q[2];
sx q[2];
rz(-2.3112042) q[2];
rz(1.7928127) q[3];
sx q[3];
rz(-1.034779) q[3];
sx q[3];
rz(-2.5428298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8778359) q[0];
sx q[0];
rz(-2.1649375) q[0];
sx q[0];
rz(2.7184955) q[0];
rz(-1.9844203) q[1];
sx q[1];
rz(-1.6559699) q[1];
sx q[1];
rz(-0.62209904) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7579278) q[0];
sx q[0];
rz(-0.94449857) q[0];
sx q[0];
rz(0.48547283) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2086925) q[2];
sx q[2];
rz(-1.1596173) q[2];
sx q[2];
rz(1.5087939) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2492952) q[1];
sx q[1];
rz(-1.163332) q[1];
sx q[1];
rz(-3.1235126) q[1];
rz(0.84526269) q[3];
sx q[3];
rz(-1.5912531) q[3];
sx q[3];
rz(1.6930228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8607934) q[2];
sx q[2];
rz(-1.7223027) q[2];
sx q[2];
rz(-2.5234176) q[2];
rz(1.482796) q[3];
sx q[3];
rz(-1.3525454) q[3];
sx q[3];
rz(-1.6331204) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16267714) q[0];
sx q[0];
rz(-1.3293043) q[0];
sx q[0];
rz(-0.83258122) q[0];
rz(-2.816448) q[1];
sx q[1];
rz(-2.145576) q[1];
sx q[1];
rz(-3.0772193) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3290117) q[0];
sx q[0];
rz(-1.4509177) q[0];
sx q[0];
rz(-2.8146312) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5369297) q[2];
sx q[2];
rz(-0.69398601) q[2];
sx q[2];
rz(-0.99550216) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2657703) q[1];
sx q[1];
rz(-2.0811305) q[1];
sx q[1];
rz(2.2052664) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0470555) q[3];
sx q[3];
rz(-0.68449293) q[3];
sx q[3];
rz(-1.5409926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9084106) q[2];
sx q[2];
rz(-0.57047129) q[2];
sx q[2];
rz(-1.2916267) q[2];
rz(0.49606797) q[3];
sx q[3];
rz(-2.3521164) q[3];
sx q[3];
rz(-1.1424278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9555776) q[0];
sx q[0];
rz(-0.29009524) q[0];
sx q[0];
rz(-2.7225323) q[0];
rz(0.31173197) q[1];
sx q[1];
rz(-1.5985951) q[1];
sx q[1];
rz(0.14351235) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3449609) q[0];
sx q[0];
rz(-1.2141742) q[0];
sx q[0];
rz(-2.659647) q[0];
rz(-2.8673929) q[2];
sx q[2];
rz(-1.7975217) q[2];
sx q[2];
rz(-2.271914) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12998768) q[1];
sx q[1];
rz(-1.2530348) q[1];
sx q[1];
rz(-2.3684273) q[1];
x q[2];
rz(-0.49893283) q[3];
sx q[3];
rz(-0.29424516) q[3];
sx q[3];
rz(-0.1732451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.781337) q[2];
sx q[2];
rz(-0.90136734) q[2];
sx q[2];
rz(2.0085013) q[2];
rz(2.0356778) q[3];
sx q[3];
rz(-1.4191041) q[3];
sx q[3];
rz(-0.47479409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7853506) q[0];
sx q[0];
rz(-2.3478735) q[0];
sx q[0];
rz(1.6957977) q[0];
rz(2.4244507) q[1];
sx q[1];
rz(-1.8627867) q[1];
sx q[1];
rz(0.76146567) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1061072) q[0];
sx q[0];
rz(-2.4940644) q[0];
sx q[0];
rz(2.8962563) q[0];
rz(-2.285073) q[2];
sx q[2];
rz(-1.1978784) q[2];
sx q[2];
rz(-2.1644217) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1902705) q[1];
sx q[1];
rz(-0.99748625) q[1];
sx q[1];
rz(0.058752937) q[1];
x q[2];
rz(-0.16038068) q[3];
sx q[3];
rz(-1.7564991) q[3];
sx q[3];
rz(-1.4726313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.36934272) q[2];
sx q[2];
rz(-0.49392924) q[2];
sx q[2];
rz(-0.93144766) q[2];
rz(2.5233968) q[3];
sx q[3];
rz(-1.5074916) q[3];
sx q[3];
rz(0.96562323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44929993) q[0];
sx q[0];
rz(-1.7789142) q[0];
sx q[0];
rz(-1.8071254) q[0];
rz(2.6284699) q[1];
sx q[1];
rz(-0.98315364) q[1];
sx q[1];
rz(1.2293053) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052027651) q[0];
sx q[0];
rz(-2.2917245) q[0];
sx q[0];
rz(-2.5682279) q[0];
x q[1];
rz(1.9077577) q[2];
sx q[2];
rz(-2.8671727) q[2];
sx q[2];
rz(1.2123337) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.38459435) q[1];
sx q[1];
rz(-1.726686) q[1];
sx q[1];
rz(-1.5916567) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66219285) q[3];
sx q[3];
rz(-1.5569038) q[3];
sx q[3];
rz(2.3628949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3076155) q[2];
sx q[2];
rz(-0.5270842) q[2];
sx q[2];
rz(-2.9311467) q[2];
rz(3.0757507) q[3];
sx q[3];
rz(-0.47271553) q[3];
sx q[3];
rz(-2.3426447) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.598572) q[0];
sx q[0];
rz(-0.6544756) q[0];
sx q[0];
rz(1.8684813) q[0];
rz(1.8674564) q[1];
sx q[1];
rz(-0.68398634) q[1];
sx q[1];
rz(2.2575016) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40128517) q[0];
sx q[0];
rz(-1.5710366) q[0];
sx q[0];
rz(0.9539286) q[0];
x q[1];
rz(-1.312125) q[2];
sx q[2];
rz(-1.0610559) q[2];
sx q[2];
rz(-2.3296251) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4834622) q[1];
sx q[1];
rz(-2.6331055) q[1];
sx q[1];
rz(0.97472753) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8330634) q[3];
sx q[3];
rz(-1.3236227) q[3];
sx q[3];
rz(2.9690774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.990443) q[2];
sx q[2];
rz(-2.1572025) q[2];
sx q[2];
rz(-1.3746064) q[2];
rz(0.19045842) q[3];
sx q[3];
rz(-2.3951267) q[3];
sx q[3];
rz(-1.1838574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9780438) q[0];
sx q[0];
rz(-1.4530285) q[0];
sx q[0];
rz(1.2646041) q[0];
rz(3.0679852) q[1];
sx q[1];
rz(-2.059748) q[1];
sx q[1];
rz(-0.61990613) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3267075) q[0];
sx q[0];
rz(-2.305079) q[0];
sx q[0];
rz(-0.22255465) q[0];
x q[1];
rz(-0.61087278) q[2];
sx q[2];
rz(-1.9985285) q[2];
sx q[2];
rz(-1.7868702) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.79110791) q[1];
sx q[1];
rz(-1.3602166) q[1];
sx q[1];
rz(-1.3665707) q[1];
rz(1.3315174) q[3];
sx q[3];
rz(-1.115723) q[3];
sx q[3];
rz(2.6309155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2747043) q[2];
sx q[2];
rz(-2.4394749) q[2];
sx q[2];
rz(-2.9956024) q[2];
rz(-2.919096) q[3];
sx q[3];
rz(-0.22495088) q[3];
sx q[3];
rz(0.72993025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9324026) q[0];
sx q[0];
rz(-1.1840191) q[0];
sx q[0];
rz(-2.9970072) q[0];
rz(-1.9119541) q[1];
sx q[1];
rz(-1.790779) q[1];
sx q[1];
rz(1.889224) q[1];
rz(-0.59496224) q[2];
sx q[2];
rz(-1.3887726) q[2];
sx q[2];
rz(-0.61979823) q[2];
rz(1.4208117) q[3];
sx q[3];
rz(-2.9290028) q[3];
sx q[3];
rz(1.3976189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
