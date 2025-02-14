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
rz(-1.8981847) q[0];
sx q[0];
rz(-1.564448) q[0];
sx q[0];
rz(2.2301883) q[0];
rz(-0.66943327) q[1];
sx q[1];
rz(-0.27639204) q[1];
sx q[1];
rz(-2.3120094) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2612635) q[0];
sx q[0];
rz(-1.9682878) q[0];
sx q[0];
rz(2.7842159) q[0];
rz(-pi) q[1];
rz(-0.35419365) q[2];
sx q[2];
rz(-2.8787432) q[2];
sx q[2];
rz(0.38989241) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1948283) q[1];
sx q[1];
rz(-0.84673893) q[1];
sx q[1];
rz(-1.9763234) q[1];
rz(-pi) q[2];
rz(2.7054993) q[3];
sx q[3];
rz(-1.4319518) q[3];
sx q[3];
rz(-2.258955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36898819) q[2];
sx q[2];
rz(-2.172894) q[2];
sx q[2];
rz(1.9770835) q[2];
rz(-2.4146967) q[3];
sx q[3];
rz(-1.3317069) q[3];
sx q[3];
rz(0.070778457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1233391) q[0];
sx q[0];
rz(-0.66272074) q[0];
sx q[0];
rz(-0.29139274) q[0];
rz(2.5413051) q[1];
sx q[1];
rz(-2.184506) q[1];
sx q[1];
rz(-0.16500638) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24848973) q[0];
sx q[0];
rz(-1.6418253) q[0];
sx q[0];
rz(-1.3079337) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8012456) q[2];
sx q[2];
rz(-2.682095) q[2];
sx q[2];
rz(0.27333096) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7126727) q[1];
sx q[1];
rz(-0.68376741) q[1];
sx q[1];
rz(3.0316975) q[1];
rz(1.6253774) q[3];
sx q[3];
rz(-1.3228058) q[3];
sx q[3];
rz(2.8021106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0773641) q[2];
sx q[2];
rz(-2.1891429) q[2];
sx q[2];
rz(-0.23951086) q[2];
rz(-2.8149878) q[3];
sx q[3];
rz(-2.9686847) q[3];
sx q[3];
rz(0.094836205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25642446) q[0];
sx q[0];
rz(-2.7056077) q[0];
sx q[0];
rz(-1.5727795) q[0];
rz(2.7478711) q[1];
sx q[1];
rz(-2.5865793) q[1];
sx q[1];
rz(2.6655925) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8825553) q[0];
sx q[0];
rz(-0.60094423) q[0];
sx q[0];
rz(1.1968234) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7287929) q[2];
sx q[2];
rz(-0.18922986) q[2];
sx q[2];
rz(-0.10052027) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6913298) q[1];
sx q[1];
rz(-2.1457377) q[1];
sx q[1];
rz(2.6460669) q[1];
rz(-pi) q[2];
x q[2];
rz(1.908444) q[3];
sx q[3];
rz(-1.8051355) q[3];
sx q[3];
rz(-0.657814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8126882) q[2];
sx q[2];
rz(-2.4220971) q[2];
sx q[2];
rz(2.1997814) q[2];
rz(1.4224667) q[3];
sx q[3];
rz(-2.7673281) q[3];
sx q[3];
rz(-2.4678738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0734237) q[0];
sx q[0];
rz(-0.24003679) q[0];
sx q[0];
rz(1.5330676) q[0];
rz(2.7948921) q[1];
sx q[1];
rz(-1.0418714) q[1];
sx q[1];
rz(2.9516721) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1153206) q[0];
sx q[0];
rz(-2.2824725) q[0];
sx q[0];
rz(1.9188966) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0937017) q[2];
sx q[2];
rz(-1.3253085) q[2];
sx q[2];
rz(-2.7141822) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9393438) q[1];
sx q[1];
rz(-0.96572438) q[1];
sx q[1];
rz(-3.0720965) q[1];
rz(0.40855405) q[3];
sx q[3];
rz(-1.9667278) q[3];
sx q[3];
rz(-2.8568639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7803663) q[2];
sx q[2];
rz(-2.2915338) q[2];
sx q[2];
rz(-0.38193199) q[2];
rz(-3.0319038) q[3];
sx q[3];
rz(-1.0332801) q[3];
sx q[3];
rz(3.0221353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1031951) q[0];
sx q[0];
rz(-1.1829475) q[0];
sx q[0];
rz(2.2947445) q[0];
rz(0.13941828) q[1];
sx q[1];
rz(-2.3250695) q[1];
sx q[1];
rz(-0.74904186) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12550929) q[0];
sx q[0];
rz(-2.2575543) q[0];
sx q[0];
rz(-1.7619049) q[0];
rz(-pi) q[1];
rz(1.7125569) q[2];
sx q[2];
rz(-1.0899915) q[2];
sx q[2];
rz(-2.9713809) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7498503) q[1];
sx q[1];
rz(-0.23255177) q[1];
sx q[1];
rz(-0.58156965) q[1];
x q[2];
rz(-1.7827634) q[3];
sx q[3];
rz(-2.6304768) q[3];
sx q[3];
rz(1.8960475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.9355115) q[2];
sx q[2];
rz(-1.8847909) q[2];
sx q[2];
rz(2.0743267) q[2];
rz(0.71077985) q[3];
sx q[3];
rz(-1.4830517) q[3];
sx q[3];
rz(1.7162286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7988605) q[0];
sx q[0];
rz(-1.6628168) q[0];
sx q[0];
rz(-2.1730098) q[0];
rz(0.69106483) q[1];
sx q[1];
rz(-1.7993118) q[1];
sx q[1];
rz(1.3074646) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1145612) q[0];
sx q[0];
rz(-1.8252715) q[0];
sx q[0];
rz(-2.4976394) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2910462) q[2];
sx q[2];
rz(-1.7564764) q[2];
sx q[2];
rz(-3.1094375) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1561574) q[1];
sx q[1];
rz(-1.623454) q[1];
sx q[1];
rz(-1.1833722) q[1];
rz(2.1627681) q[3];
sx q[3];
rz(-2.8205964) q[3];
sx q[3];
rz(0.61883607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6664194) q[2];
sx q[2];
rz(-2.9066777) q[2];
sx q[2];
rz(-1.2972181) q[2];
rz(1.3951067) q[3];
sx q[3];
rz(-1.3336983) q[3];
sx q[3];
rz(1.5313799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.34201) q[0];
sx q[0];
rz(-0.81556773) q[0];
sx q[0];
rz(-2.2014501) q[0];
rz(-0.6483342) q[1];
sx q[1];
rz(-1.2038566) q[1];
sx q[1];
rz(2.8727093) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20737442) q[0];
sx q[0];
rz(-1.3257799) q[0];
sx q[0];
rz(-2.5859358) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70272311) q[2];
sx q[2];
rz(-1.8112735) q[2];
sx q[2];
rz(-0.81288494) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9098879) q[1];
sx q[1];
rz(-0.99867601) q[1];
sx q[1];
rz(0.27426274) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65461378) q[3];
sx q[3];
rz(-2.6179492) q[3];
sx q[3];
rz(-0.80010993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1342643) q[2];
sx q[2];
rz(-1.8625883) q[2];
sx q[2];
rz(0.68592611) q[2];
rz(0.736233) q[3];
sx q[3];
rz(-1.0400306) q[3];
sx q[3];
rz(-2.9161016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4418884) q[0];
sx q[0];
rz(-1.9206973) q[0];
sx q[0];
rz(0.7487444) q[0];
rz(0.092983149) q[1];
sx q[1];
rz(-0.75983202) q[1];
sx q[1];
rz(-2.0704796) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7629246) q[0];
sx q[0];
rz(-1.0707741) q[0];
sx q[0];
rz(-2.0471694) q[0];
x q[1];
rz(2.5701809) q[2];
sx q[2];
rz(-0.73941747) q[2];
sx q[2];
rz(2.6020056) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.86102137) q[1];
sx q[1];
rz(-1.8771267) q[1];
sx q[1];
rz(0.88545274) q[1];
x q[2];
rz(-1.981826) q[3];
sx q[3];
rz(-1.618652) q[3];
sx q[3];
rz(1.4519297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6894655) q[2];
sx q[2];
rz(-0.63673821) q[2];
sx q[2];
rz(-2.9316736) q[2];
rz(-0.12957761) q[3];
sx q[3];
rz(-0.94730535) q[3];
sx q[3];
rz(1.5642222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0196446) q[0];
sx q[0];
rz(-0.26109281) q[0];
sx q[0];
rz(3.1280532) q[0];
rz(0.53567046) q[1];
sx q[1];
rz(-2.5745013) q[1];
sx q[1];
rz(-0.013669107) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5807728) q[0];
sx q[0];
rz(-0.81308621) q[0];
sx q[0];
rz(1.2570501) q[0];
rz(1.0837567) q[2];
sx q[2];
rz(-1.9197316) q[2];
sx q[2];
rz(-1.2020122) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0993578) q[1];
sx q[1];
rz(-0.7868979) q[1];
sx q[1];
rz(-0.46532793) q[1];
rz(-pi) q[2];
rz(-2.3914371) q[3];
sx q[3];
rz(-1.292308) q[3];
sx q[3];
rz(0.089547308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.51836625) q[2];
sx q[2];
rz(-1.7609111) q[2];
sx q[2];
rz(-1.3631442) q[2];
rz(2.3220883) q[3];
sx q[3];
rz(-2.6624811) q[3];
sx q[3];
rz(-0.68377408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0226456) q[0];
sx q[0];
rz(-0.85856694) q[0];
sx q[0];
rz(1.488142) q[0];
rz(2.7986774) q[1];
sx q[1];
rz(-1.5577134) q[1];
sx q[1];
rz(-2.3905579) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58758754) q[0];
sx q[0];
rz(-2.8556654) q[0];
sx q[0];
rz(1.343973) q[0];
rz(-pi) q[1];
rz(1.1840759) q[2];
sx q[2];
rz(-0.49418517) q[2];
sx q[2];
rz(-1.7750211) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0768834) q[1];
sx q[1];
rz(-0.61184498) q[1];
sx q[1];
rz(2.9348899) q[1];
x q[2];
rz(3.1084177) q[3];
sx q[3];
rz(-1.4806595) q[3];
sx q[3];
rz(1.7189004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.670383) q[2];
sx q[2];
rz(-2.1977916) q[2];
sx q[2];
rz(-0.38718265) q[2];
rz(-2.6662628) q[3];
sx q[3];
rz(-1.1673704) q[3];
sx q[3];
rz(-0.79132426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.249007) q[0];
sx q[0];
rz(-0.96207608) q[0];
sx q[0];
rz(0.67208653) q[0];
rz(0.36920209) q[1];
sx q[1];
rz(-0.75405706) q[1];
sx q[1];
rz(-1.3934607) q[1];
rz(-0.71699981) q[2];
sx q[2];
rz(-1.430027) q[2];
sx q[2];
rz(1.2363557) q[2];
rz(-1.8440856) q[3];
sx q[3];
rz(-1.6798301) q[3];
sx q[3];
rz(-1.9435431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
