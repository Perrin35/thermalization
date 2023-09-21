OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0857467) q[0];
sx q[0];
rz(-0.081781713) q[0];
sx q[0];
rz(-2.6401289) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(0.3224386) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3381391) q[0];
sx q[0];
rz(-0.33284602) q[0];
sx q[0];
rz(1.8344318) q[0];
x q[1];
rz(2.5311004) q[2];
sx q[2];
rz(-0.98449003) q[2];
sx q[2];
rz(2.5551978) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.3024219) q[1];
sx q[1];
rz(-1.9170554) q[1];
sx q[1];
rz(0.94567169) q[1];
rz(-pi) q[2];
rz(-0.97500719) q[3];
sx q[3];
rz(-2.4132204) q[3];
sx q[3];
rz(-1.7155852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6364608) q[2];
sx q[2];
rz(-0.59288609) q[2];
sx q[2];
rz(0.55603975) q[2];
rz(-2.3089144) q[3];
sx q[3];
rz(-1.6502389) q[3];
sx q[3];
rz(0.94579831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44822025) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(0.15727501) q[0];
rz(-0.26113025) q[1];
sx q[1];
rz(-1.3477247) q[1];
sx q[1];
rz(-0.10903407) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63576525) q[0];
sx q[0];
rz(-0.286239) q[0];
sx q[0];
rz(-2.2155227) q[0];
rz(-2.878506) q[2];
sx q[2];
rz(-1.7619942) q[2];
sx q[2];
rz(-2.4441602) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8728767) q[1];
sx q[1];
rz(-2.0375588) q[1];
sx q[1];
rz(-0.66653911) q[1];
rz(-pi) q[2];
rz(-0.29626366) q[3];
sx q[3];
rz(-0.54013541) q[3];
sx q[3];
rz(1.4512645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9033501) q[2];
sx q[2];
rz(-1.1652596) q[2];
sx q[2];
rz(1.2634574) q[2];
rz(-0.3271099) q[3];
sx q[3];
rz(-1.5771022) q[3];
sx q[3];
rz(1.9272778) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064421244) q[0];
sx q[0];
rz(-0.049296878) q[0];
sx q[0];
rz(-1.7984614) q[0];
rz(-2.893977) q[1];
sx q[1];
rz(-2.394948) q[1];
sx q[1];
rz(-2.6599191) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96268481) q[0];
sx q[0];
rz(-2.509153) q[0];
sx q[0];
rz(-2.2326367) q[0];
x q[1];
rz(0.33294296) q[2];
sx q[2];
rz(-2.1579086) q[2];
sx q[2];
rz(1.743403) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2783918) q[1];
sx q[1];
rz(-0.87676261) q[1];
sx q[1];
rz(-2.4064526) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2948202) q[3];
sx q[3];
rz(-2.0623042) q[3];
sx q[3];
rz(1.5566952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8032916) q[2];
sx q[2];
rz(-2.3241966) q[2];
sx q[2];
rz(-2.6417007) q[2];
rz(-2.5806184) q[3];
sx q[3];
rz(-1.2596954) q[3];
sx q[3];
rz(-1.6803754) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19701476) q[0];
sx q[0];
rz(-0.16600969) q[0];
sx q[0];
rz(-0.552185) q[0];
rz(1.5532956) q[1];
sx q[1];
rz(-0.89893666) q[1];
sx q[1];
rz(1.2447371) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9024076) q[0];
sx q[0];
rz(-0.33948487) q[0];
sx q[0];
rz(-0.0048588077) q[0];
x q[1];
rz(2.142698) q[2];
sx q[2];
rz(-2.9036387) q[2];
sx q[2];
rz(1.9902802) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1340449) q[1];
sx q[1];
rz(-1.3589077) q[1];
sx q[1];
rz(1.0210387) q[1];
rz(-2.7987715) q[3];
sx q[3];
rz(-2.5431513) q[3];
sx q[3];
rz(0.10859057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.84919471) q[2];
sx q[2];
rz(-1.8820102) q[2];
sx q[2];
rz(-1.1506895) q[2];
rz(-1.6644647) q[3];
sx q[3];
rz(-1.5090347) q[3];
sx q[3];
rz(0.48294827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078159049) q[0];
sx q[0];
rz(-0.76197356) q[0];
sx q[0];
rz(0.081469014) q[0];
rz(-0.062462656) q[1];
sx q[1];
rz(-2.000258) q[1];
sx q[1];
rz(1.6385471) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8836356) q[0];
sx q[0];
rz(-0.99146087) q[0];
sx q[0];
rz(2.9888319) q[0];
rz(0.083085255) q[2];
sx q[2];
rz(-1.462888) q[2];
sx q[2];
rz(-1.2667058) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.11322184) q[1];
sx q[1];
rz(-1.2591259) q[1];
sx q[1];
rz(0.10722864) q[1];
rz(-pi) q[2];
rz(-1.2522069) q[3];
sx q[3];
rz(-2.0584403) q[3];
sx q[3];
rz(0.82061758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.5082671) q[2];
sx q[2];
rz(-0.94255629) q[2];
sx q[2];
rz(1.903669) q[2];
rz(-1.1226908) q[3];
sx q[3];
rz(-0.676238) q[3];
sx q[3];
rz(-0.52156633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6102819) q[0];
sx q[0];
rz(-1.0045445) q[0];
sx q[0];
rz(-0.26671985) q[0];
rz(0.56089127) q[1];
sx q[1];
rz(-1.8436878) q[1];
sx q[1];
rz(0.7985324) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.610299) q[0];
sx q[0];
rz(-1.7113422) q[0];
sx q[0];
rz(-1.2929686) q[0];
x q[1];
rz(0.59136765) q[2];
sx q[2];
rz(-2.5632576) q[2];
sx q[2];
rz(2.8665198) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.27998566) q[1];
sx q[1];
rz(-2.1415188) q[1];
sx q[1];
rz(-0.062203783) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3409307) q[3];
sx q[3];
rz(-2.0293529) q[3];
sx q[3];
rz(0.82160219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16053998) q[2];
sx q[2];
rz(-1.8926228) q[2];
sx q[2];
rz(0.36671656) q[2];
rz(1.8803053) q[3];
sx q[3];
rz(-1.6882608) q[3];
sx q[3];
rz(3.0453851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7325608) q[0];
sx q[0];
rz(-0.92027396) q[0];
sx q[0];
rz(0.60638705) q[0];
rz(0.19730332) q[1];
sx q[1];
rz(-1.1261255) q[1];
sx q[1];
rz(-0.46404776) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2460829) q[0];
sx q[0];
rz(-1.9946788) q[0];
sx q[0];
rz(-0.8830107) q[0];
rz(-1.8115933) q[2];
sx q[2];
rz(-1.5903928) q[2];
sx q[2];
rz(-0.400825) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.11394994) q[1];
sx q[1];
rz(-1.0877891) q[1];
sx q[1];
rz(2.7999858) q[1];
rz(-0.32902284) q[3];
sx q[3];
rz(-2.892422) q[3];
sx q[3];
rz(2.3014625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.93418926) q[2];
sx q[2];
rz(-1.0031909) q[2];
sx q[2];
rz(2.8835473) q[2];
rz(1.1856273) q[3];
sx q[3];
rz(-1.5153171) q[3];
sx q[3];
rz(-0.0035088249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19514062) q[0];
sx q[0];
rz(-1.2807245) q[0];
sx q[0];
rz(2.7602957) q[0];
rz(3.0463468) q[1];
sx q[1];
rz(-0.97243273) q[1];
sx q[1];
rz(1.7000748) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9593175) q[0];
sx q[0];
rz(-1.5858985) q[0];
sx q[0];
rz(3.0781151) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0870073) q[2];
sx q[2];
rz(-0.61331257) q[2];
sx q[2];
rz(1.7577946) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.085233363) q[1];
sx q[1];
rz(-1.4409522) q[1];
sx q[1];
rz(-1.4443881) q[1];
x q[2];
rz(-2.8076595) q[3];
sx q[3];
rz(-1.6113558) q[3];
sx q[3];
rz(-3.0143152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9986481) q[2];
sx q[2];
rz(-0.412985) q[2];
sx q[2];
rz(-0.22658919) q[2];
rz(-0.44858027) q[3];
sx q[3];
rz(-1.6058763) q[3];
sx q[3];
rz(-0.8297689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7609693) q[0];
sx q[0];
rz(-2.3801104) q[0];
sx q[0];
rz(1.3990078) q[0];
rz(-2.8245068) q[1];
sx q[1];
rz(-1.4750907) q[1];
sx q[1];
rz(2.1549966) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6331659) q[0];
sx q[0];
rz(-1.0796483) q[0];
sx q[0];
rz(-1.7500061) q[0];
rz(-pi) q[1];
rz(-0.53194745) q[2];
sx q[2];
rz(-0.94594687) q[2];
sx q[2];
rz(-0.22235409) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.66227312) q[1];
sx q[1];
rz(-0.62367491) q[1];
sx q[1];
rz(-0.24346607) q[1];
rz(-1.9429728) q[3];
sx q[3];
rz(-2.6937006) q[3];
sx q[3];
rz(0.4263634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2150779) q[2];
sx q[2];
rz(-2.4145917) q[2];
sx q[2];
rz(-0.40965664) q[2];
rz(2.8783197) q[3];
sx q[3];
rz(-1.8299088) q[3];
sx q[3];
rz(2.6221361) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7423994) q[0];
sx q[0];
rz(-0.078646794) q[0];
sx q[0];
rz(1.7364527) q[0];
rz(0.82110226) q[1];
sx q[1];
rz(-2.2228873) q[1];
sx q[1];
rz(1.4155037) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1499407) q[0];
sx q[0];
rz(-1.2399925) q[0];
sx q[0];
rz(-1.3209692) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9344994) q[2];
sx q[2];
rz(-1.8038245) q[2];
sx q[2];
rz(1.9050913) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0714598) q[1];
sx q[1];
rz(-1.7839285) q[1];
sx q[1];
rz(0.15503426) q[1];
rz(-2.9726082) q[3];
sx q[3];
rz(-2.0770693) q[3];
sx q[3];
rz(-0.81912012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.71904174) q[2];
sx q[2];
rz(-2.8406403) q[2];
sx q[2];
rz(3.0174875) q[2];
rz(-2.1758046) q[3];
sx q[3];
rz(-1.6671168) q[3];
sx q[3];
rz(0.66108274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60349764) q[0];
sx q[0];
rz(-2.8932543) q[0];
sx q[0];
rz(2.2809991) q[0];
rz(-2.8339236) q[1];
sx q[1];
rz(-1.888702) q[1];
sx q[1];
rz(-1.9370334) q[1];
rz(-1.7514501) q[2];
sx q[2];
rz(-1.3144819) q[2];
sx q[2];
rz(-1.1983295) q[2];
rz(-0.23871213) q[3];
sx q[3];
rz(-0.56959) q[3];
sx q[3];
rz(0.068985229) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
