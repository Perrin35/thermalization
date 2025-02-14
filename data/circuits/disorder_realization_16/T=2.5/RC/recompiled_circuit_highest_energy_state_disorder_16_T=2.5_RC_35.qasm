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
rz(0.056869153) q[0];
sx q[0];
rz(-0.19357227) q[0];
sx q[0];
rz(-0.40785664) q[0];
rz(-0.017539311) q[1];
sx q[1];
rz(-1.8899625) q[1];
sx q[1];
rz(-1.6012021) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9303794) q[0];
sx q[0];
rz(-1.5222852) q[0];
sx q[0];
rz(2.0082974) q[0];
rz(-pi) q[1];
rz(-2.8694334) q[2];
sx q[2];
rz(-0.46900422) q[2];
sx q[2];
rz(-2.9763165) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0527264) q[1];
sx q[1];
rz(-1.7824934) q[1];
sx q[1];
rz(-2.746341) q[1];
rz(-2.3204221) q[3];
sx q[3];
rz(-1.0945012) q[3];
sx q[3];
rz(-1.1007635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5876329) q[2];
sx q[2];
rz(-1.0593869) q[2];
sx q[2];
rz(-0.21767347) q[2];
rz(0.31146464) q[3];
sx q[3];
rz(-0.61190999) q[3];
sx q[3];
rz(-1.1658839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72164732) q[0];
sx q[0];
rz(-2.8657275) q[0];
sx q[0];
rz(-2.4496147) q[0];
rz(-0.75633374) q[1];
sx q[1];
rz(-2.6192009) q[1];
sx q[1];
rz(1.5501032) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1082008) q[0];
sx q[0];
rz(-2.4188359) q[0];
sx q[0];
rz(-2.9563144) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0372411) q[2];
sx q[2];
rz(-1.1747825) q[2];
sx q[2];
rz(1.9463152) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.34775388) q[1];
sx q[1];
rz(-1.4112541) q[1];
sx q[1];
rz(-1.4473587) q[1];
rz(-0.85568537) q[3];
sx q[3];
rz(-2.4484903) q[3];
sx q[3];
rz(0.86831304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99593607) q[2];
sx q[2];
rz(-2.9684976) q[2];
sx q[2];
rz(0.23925979) q[2];
rz(-1.4907106) q[3];
sx q[3];
rz(-1.2942634) q[3];
sx q[3];
rz(1.2283121) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8656798) q[0];
sx q[0];
rz(-2.2353807) q[0];
sx q[0];
rz(0.014634125) q[0];
rz(-2.2485661) q[1];
sx q[1];
rz(-2.1664186) q[1];
sx q[1];
rz(-0.21569529) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6071226) q[0];
sx q[0];
rz(-1.5754345) q[0];
sx q[0];
rz(-1.5195373) q[0];
rz(-pi) q[1];
rz(1.9142308) q[2];
sx q[2];
rz(-2.2015155) q[2];
sx q[2];
rz(1.3892) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.025192997) q[1];
sx q[1];
rz(-1.6982251) q[1];
sx q[1];
rz(-2.9078632) q[1];
rz(-pi) q[2];
rz(-2.3560432) q[3];
sx q[3];
rz(-2.2758898) q[3];
sx q[3];
rz(0.37032235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0594242) q[2];
sx q[2];
rz(-1.1801722) q[2];
sx q[2];
rz(-0.99349418) q[2];
rz(-0.70837402) q[3];
sx q[3];
rz(-2.0230484) q[3];
sx q[3];
rz(0.12349252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3941536) q[0];
sx q[0];
rz(-0.49462947) q[0];
sx q[0];
rz(-0.69394773) q[0];
rz(-2.8548062) q[1];
sx q[1];
rz(-1.5319805) q[1];
sx q[1];
rz(2.0538816) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3932753) q[0];
sx q[0];
rz(-1.4469742) q[0];
sx q[0];
rz(2.078915) q[0];
rz(-pi) q[1];
rz(-2.3871586) q[2];
sx q[2];
rz(-2.4409119) q[2];
sx q[2];
rz(2.5452819) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.81146152) q[1];
sx q[1];
rz(-0.68492678) q[1];
sx q[1];
rz(-2.0622938) q[1];
rz(-pi) q[2];
rz(-3.089043) q[3];
sx q[3];
rz(-0.71501117) q[3];
sx q[3];
rz(2.8831567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.22471681) q[2];
sx q[2];
rz(-1.9852873) q[2];
sx q[2];
rz(-1.3678331) q[2];
rz(1.9035089) q[3];
sx q[3];
rz(-1.5690469) q[3];
sx q[3];
rz(0.21627538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5991768) q[0];
sx q[0];
rz(-1.8734064) q[0];
sx q[0];
rz(-0.61781484) q[0];
rz(-1.1082331) q[1];
sx q[1];
rz(-2.8384659) q[1];
sx q[1];
rz(2.2584426) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1533494) q[0];
sx q[0];
rz(-0.81479615) q[0];
sx q[0];
rz(1.5809466) q[0];
x q[1];
rz(1.7332212) q[2];
sx q[2];
rz(-2.0018296) q[2];
sx q[2];
rz(-1.169432) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1146436) q[1];
sx q[1];
rz(-2.0609239) q[1];
sx q[1];
rz(3.0672706) q[1];
rz(-pi) q[2];
rz(-2.1260902) q[3];
sx q[3];
rz(-1.9262553) q[3];
sx q[3];
rz(2.773284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9250179) q[2];
sx q[2];
rz(-2.8897132) q[2];
sx q[2];
rz(1.3612755) q[2];
rz(1.8682293) q[3];
sx q[3];
rz(-1.7073771) q[3];
sx q[3];
rz(0.98289615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0971766) q[0];
sx q[0];
rz(-1.2460848) q[0];
sx q[0];
rz(0.62028766) q[0];
rz(-0.39707956) q[1];
sx q[1];
rz(-0.47080165) q[1];
sx q[1];
rz(1.9420067) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37017779) q[0];
sx q[0];
rz(-1.4645029) q[0];
sx q[0];
rz(-3.1205157) q[0];
x q[1];
rz(-1.7825837) q[2];
sx q[2];
rz(-1.7149936) q[2];
sx q[2];
rz(-1.5825001) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0565655) q[1];
sx q[1];
rz(-2.0904087) q[1];
sx q[1];
rz(1.3481793) q[1];
rz(-1.4523466) q[3];
sx q[3];
rz(-1.9129686) q[3];
sx q[3];
rz(2.7717048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.24641307) q[2];
sx q[2];
rz(-1.2440888) q[2];
sx q[2];
rz(2.4708774) q[2];
rz(0.29117584) q[3];
sx q[3];
rz(-0.61101919) q[3];
sx q[3];
rz(2.5197855) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7156242) q[0];
sx q[0];
rz(-1.427303) q[0];
sx q[0];
rz(-0.17844644) q[0];
rz(-2.2715691) q[1];
sx q[1];
rz(-0.88302892) q[1];
sx q[1];
rz(-2.3387486) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3960796) q[0];
sx q[0];
rz(-0.95142309) q[0];
sx q[0];
rz(1.5399163) q[0];
rz(-pi) q[1];
rz(-0.3608256) q[2];
sx q[2];
rz(-0.61491167) q[2];
sx q[2];
rz(0.71984839) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5685988) q[1];
sx q[1];
rz(-0.90403176) q[1];
sx q[1];
rz(1.1519037) q[1];
rz(-2.1986897) q[3];
sx q[3];
rz(-1.3602837) q[3];
sx q[3];
rz(0.26417662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7650083) q[2];
sx q[2];
rz(-1.2976126) q[2];
sx q[2];
rz(-2.2479748) q[2];
rz(-1.5417967) q[3];
sx q[3];
rz(-1.6518281) q[3];
sx q[3];
rz(2.5347575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97745085) q[0];
sx q[0];
rz(-2.7315388) q[0];
sx q[0];
rz(-2.9264911) q[0];
rz(-2.2970301) q[1];
sx q[1];
rz(-1.3203878) q[1];
sx q[1];
rz(0.79427687) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0117019) q[0];
sx q[0];
rz(-0.029057682) q[0];
sx q[0];
rz(-0.41883166) q[0];
rz(-1.2634981) q[2];
sx q[2];
rz(-1.1329044) q[2];
sx q[2];
rz(-1.4222933) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.25674838) q[1];
sx q[1];
rz(-1.0326385) q[1];
sx q[1];
rz(-3.0358008) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78171697) q[3];
sx q[3];
rz(-1.3355044) q[3];
sx q[3];
rz(0.58636753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8636785) q[2];
sx q[2];
rz(-1.9387551) q[2];
sx q[2];
rz(-1.4097144) q[2];
rz(2.388741) q[3];
sx q[3];
rz(-1.9949621) q[3];
sx q[3];
rz(1.0021771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79733217) q[0];
sx q[0];
rz(-2.7334038) q[0];
sx q[0];
rz(-2.1687188) q[0];
rz(1.6196039) q[1];
sx q[1];
rz(-1.3643967) q[1];
sx q[1];
rz(-0.63046986) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4578188) q[0];
sx q[0];
rz(-0.22916482) q[0];
sx q[0];
rz(0.650371) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7277355) q[2];
sx q[2];
rz(-1.5633768) q[2];
sx q[2];
rz(-0.8366303) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3296632) q[1];
sx q[1];
rz(-2.3721937) q[1];
sx q[1];
rz(-0.071746221) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34548605) q[3];
sx q[3];
rz(-0.18426963) q[3];
sx q[3];
rz(0.94217448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1260599) q[2];
sx q[2];
rz(-1.1522747) q[2];
sx q[2];
rz(1.3014334) q[2];
rz(-1.5850916) q[3];
sx q[3];
rz(-2.3157178) q[3];
sx q[3];
rz(-2.7263156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5179317) q[0];
sx q[0];
rz(-2.9908337) q[0];
sx q[0];
rz(-0.49322042) q[0];
rz(-1.8863691) q[1];
sx q[1];
rz(-1.8938277) q[1];
sx q[1];
rz(2.7769322) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4987558) q[0];
sx q[0];
rz(-0.36223534) q[0];
sx q[0];
rz(-1.1722159) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91835351) q[2];
sx q[2];
rz(-1.9404354) q[2];
sx q[2];
rz(0.015394848) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4019822) q[1];
sx q[1];
rz(-1.6270042) q[1];
sx q[1];
rz(-2.3422041) q[1];
x q[2];
rz(1.4480035) q[3];
sx q[3];
rz(-2.8229575) q[3];
sx q[3];
rz(-0.84240571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4517335) q[2];
sx q[2];
rz(-2.1351337) q[2];
sx q[2];
rz(-0.9922007) q[2];
rz(-0.36710468) q[3];
sx q[3];
rz(-1.7016943) q[3];
sx q[3];
rz(2.4632857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.622396) q[0];
sx q[0];
rz(-1.665103) q[0];
sx q[0];
rz(-1.7892224) q[0];
rz(0.92338152) q[1];
sx q[1];
rz(-0.8538178) q[1];
sx q[1];
rz(-1.1710844) q[1];
rz(0.14262234) q[2];
sx q[2];
rz(-1.3726948) q[2];
sx q[2];
rz(1.2683327) q[2];
rz(0.97040776) q[3];
sx q[3];
rz(-2.1526489) q[3];
sx q[3];
rz(-2.7027705) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
