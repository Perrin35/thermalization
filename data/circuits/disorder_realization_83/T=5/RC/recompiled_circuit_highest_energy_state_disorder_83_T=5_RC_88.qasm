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
rz(1.243408) q[0];
sx q[0];
rz(4.7060407) q[0];
sx q[0];
rz(10.336182) q[0];
rz(-0.66943327) q[1];
sx q[1];
rz(-0.27639204) q[1];
sx q[1];
rz(0.82958329) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54691168) q[0];
sx q[0];
rz(-1.24238) q[0];
sx q[0];
rz(-1.1494779) q[0];
rz(-pi) q[1];
rz(0.24721036) q[2];
sx q[2];
rz(-1.6610378) q[2];
sx q[2];
rz(-1.6177141) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7698698) q[1];
sx q[1];
rz(-2.3300305) q[1];
sx q[1];
rz(0.4196336) q[1];
rz(-pi) q[2];
rz(-1.7237648) q[3];
sx q[3];
rz(-2.0024096) q[3];
sx q[3];
rz(-0.6237517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.36898819) q[2];
sx q[2];
rz(-2.172894) q[2];
sx q[2];
rz(1.1645092) q[2];
rz(-2.4146967) q[3];
sx q[3];
rz(-1.3317069) q[3];
sx q[3];
rz(0.070778457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1233391) q[0];
sx q[0];
rz(-0.66272074) q[0];
sx q[0];
rz(0.29139274) q[0];
rz(2.5413051) q[1];
sx q[1];
rz(-0.95708668) q[1];
sx q[1];
rz(0.16500638) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8383808) q[0];
sx q[0];
rz(-1.3086122) q[0];
sx q[0];
rz(-3.068046) q[0];
rz(3.0290452) q[2];
sx q[2];
rz(-1.1243382) q[2];
sx q[2];
rz(-3.1242862) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7126727) q[1];
sx q[1];
rz(-0.68376741) q[1];
sx q[1];
rz(-0.10989519) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8932472) q[3];
sx q[3];
rz(-1.5178866) q[3];
sx q[3];
rz(1.2179045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0642285) q[2];
sx q[2];
rz(-2.1891429) q[2];
sx q[2];
rz(-2.9020818) q[2];
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
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25642446) q[0];
sx q[0];
rz(-2.7056077) q[0];
sx q[0];
rz(-1.5727795) q[0];
rz(-2.7478711) q[1];
sx q[1];
rz(-2.5865793) q[1];
sx q[1];
rz(0.47600019) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18506348) q[0];
sx q[0];
rz(-2.1250884) q[0];
sx q[0];
rz(0.24538641) q[0];
rz(-1.3838686) q[2];
sx q[2];
rz(-1.5411959) q[2];
sx q[2];
rz(-1.6254978) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2326374) q[1];
sx q[1];
rz(-2.401315) q[1];
sx q[1];
rz(-2.2038682) q[1];
rz(-pi) q[2];
rz(-1.2331486) q[3];
sx q[3];
rz(-1.8051355) q[3];
sx q[3];
rz(-0.657814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8126882) q[2];
sx q[2];
rz(-2.4220971) q[2];
sx q[2];
rz(2.1997814) q[2];
rz(-1.719126) q[3];
sx q[3];
rz(-0.37426451) q[3];
sx q[3];
rz(-0.67371887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.068168966) q[0];
sx q[0];
rz(-2.9015559) q[0];
sx q[0];
rz(1.6085251) q[0];
rz(0.34670058) q[1];
sx q[1];
rz(-1.0418714) q[1];
sx q[1];
rz(0.18992058) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53341928) q[0];
sx q[0];
rz(-0.77866879) q[0];
sx q[0];
rz(0.37663583) q[0];
x q[1];
rz(1.0713351) q[2];
sx q[2];
rz(-2.609406) q[2];
sx q[2];
rz(2.4376873) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0611113) q[1];
sx q[1];
rz(-0.60855344) q[1];
sx q[1];
rz(1.670865) q[1];
x q[2];
rz(2.7330386) q[3];
sx q[3];
rz(-1.9667278) q[3];
sx q[3];
rz(2.8568639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7803663) q[2];
sx q[2];
rz(-2.2915338) q[2];
sx q[2];
rz(2.7596607) q[2];
rz(3.0319038) q[3];
sx q[3];
rz(-1.0332801) q[3];
sx q[3];
rz(-3.0221353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0383976) q[0];
sx q[0];
rz(-1.1829475) q[0];
sx q[0];
rz(-0.84684816) q[0];
rz(0.13941828) q[1];
sx q[1];
rz(-2.3250695) q[1];
sx q[1];
rz(2.3925508) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0160834) q[0];
sx q[0];
rz(-2.2575543) q[0];
sx q[0];
rz(1.3796877) q[0];
x q[1];
rz(1.7125569) q[2];
sx q[2];
rz(-2.0516011) q[2];
sx q[2];
rz(2.9713809) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3440123) q[1];
sx q[1];
rz(-1.3770119) q[1];
sx q[1];
rz(-1.4414201) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3588293) q[3];
sx q[3];
rz(-2.6304768) q[3];
sx q[3];
rz(-1.2455452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.9355115) q[2];
sx q[2];
rz(-1.2568018) q[2];
sx q[2];
rz(-1.067266) q[2];
rz(-2.4308128) q[3];
sx q[3];
rz(-1.4830517) q[3];
sx q[3];
rz(-1.425364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34273219) q[0];
sx q[0];
rz(-1.4787759) q[0];
sx q[0];
rz(0.96858281) q[0];
rz(2.4505278) q[1];
sx q[1];
rz(-1.7993118) q[1];
sx q[1];
rz(1.8341281) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7846061) q[0];
sx q[0];
rz(-2.1907867) q[0];
sx q[0];
rz(-1.2563353) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2910462) q[2];
sx q[2];
rz(-1.3851162) q[2];
sx q[2];
rz(3.1094375) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56388748) q[1];
sx q[1];
rz(-1.9576549) q[1];
sx q[1];
rz(3.0847286) q[1];
rz(-pi) q[2];
rz(1.8400165) q[3];
sx q[3];
rz(-1.3938187) q[3];
sx q[3];
rz(0.38401803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.4751733) q[2];
sx q[2];
rz(-2.9066777) q[2];
sx q[2];
rz(-1.2972181) q[2];
rz(-1.3951067) q[3];
sx q[3];
rz(-1.8078943) q[3];
sx q[3];
rz(1.5313799) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7995826) q[0];
sx q[0];
rz(-2.3260249) q[0];
sx q[0];
rz(0.94014257) q[0];
rz(-2.4932585) q[1];
sx q[1];
rz(-1.2038566) q[1];
sx q[1];
rz(0.26888332) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2139287) q[0];
sx q[0];
rz(-1.0335796) q[0];
sx q[0];
rz(1.8570333) q[0];
rz(-pi) q[1];
rz(-2.4388695) q[2];
sx q[2];
rz(-1.3303192) q[2];
sx q[2];
rz(0.81288494) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9536679) q[1];
sx q[1];
rz(-1.3410733) q[1];
sx q[1];
rz(-0.98120316) q[1];
rz(-2.712065) q[3];
sx q[3];
rz(-1.2614354) q[3];
sx q[3];
rz(1.7842899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1342643) q[2];
sx q[2];
rz(-1.8625883) q[2];
sx q[2];
rz(2.4556665) q[2];
rz(2.4053597) q[3];
sx q[3];
rz(-1.0400306) q[3];
sx q[3];
rz(-0.22549103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69970423) q[0];
sx q[0];
rz(-1.9206973) q[0];
sx q[0];
rz(-2.3928483) q[0];
rz(0.092983149) q[1];
sx q[1];
rz(-0.75983202) q[1];
sx q[1];
rz(-2.0704796) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94959983) q[0];
sx q[0];
rz(-1.9849536) q[0];
sx q[0];
rz(2.5903775) q[0];
rz(0.57141177) q[2];
sx q[2];
rz(-2.4021752) q[2];
sx q[2];
rz(-0.53958708) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0784356) q[1];
sx q[1];
rz(-0.74045762) q[1];
sx q[1];
rz(1.1073807) q[1];
rz(-pi) q[2];
x q[2];
rz(1.451501) q[3];
sx q[3];
rz(-0.41364851) q[3];
sx q[3];
rz(0.0095490015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6894655) q[2];
sx q[2];
rz(-2.5048544) q[2];
sx q[2];
rz(-2.9316736) q[2];
rz(0.12957761) q[3];
sx q[3];
rz(-0.94730535) q[3];
sx q[3];
rz(-1.5642222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.121948) q[0];
sx q[0];
rz(-2.8804998) q[0];
sx q[0];
rz(-0.013539465) q[0];
rz(2.6059222) q[1];
sx q[1];
rz(-0.56709138) q[1];
sx q[1];
rz(3.1279235) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5807728) q[0];
sx q[0];
rz(-0.81308621) q[0];
sx q[0];
rz(-1.8845425) q[0];
rz(-pi) q[1];
rz(-0.39054587) q[2];
sx q[2];
rz(-2.0261814) q[2];
sx q[2];
rz(-0.18965882) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7175258) q[1];
sx q[1];
rz(-2.2560511) q[1];
sx q[1];
rz(-1.1478893) q[1];
rz(-pi) q[2];
rz(-2.7444856) q[3];
sx q[3];
rz(-2.350961) q[3];
sx q[3];
rz(-1.9472029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.51836625) q[2];
sx q[2];
rz(-1.3806815) q[2];
sx q[2];
rz(-1.7784485) q[2];
rz(-2.3220883) q[3];
sx q[3];
rz(-2.6624811) q[3];
sx q[3];
rz(-2.4578186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.118947) q[0];
sx q[0];
rz(-0.85856694) q[0];
sx q[0];
rz(1.488142) q[0];
rz(-2.7986774) q[1];
sx q[1];
rz(-1.5838793) q[1];
sx q[1];
rz(0.75103474) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5540051) q[0];
sx q[0];
rz(-0.28592725) q[0];
sx q[0];
rz(-1.343973) q[0];
x q[1];
rz(-2.9411211) q[2];
sx q[2];
rz(-2.0256009) q[2];
sx q[2];
rz(-1.7998296) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.955525) q[1];
sx q[1];
rz(-0.97381089) q[1];
sx q[1];
rz(1.7138193) q[1];
rz(-pi) q[2];
rz(-1.6609825) q[3];
sx q[3];
rz(-1.6038365) q[3];
sx q[3];
rz(0.15109135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.47120961) q[2];
sx q[2];
rz(-0.94380108) q[2];
sx q[2];
rz(-2.75441) q[2];
rz(0.47532982) q[3];
sx q[3];
rz(-1.9742222) q[3];
sx q[3];
rz(0.79132426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8925856) q[0];
sx q[0];
rz(-2.1795166) q[0];
sx q[0];
rz(-2.4695061) q[0];
rz(0.36920209) q[1];
sx q[1];
rz(-0.75405706) q[1];
sx q[1];
rz(-1.3934607) q[1];
rz(1.7566219) q[2];
sx q[2];
rz(-0.86238774) q[2];
sx q[2];
rz(2.9288616) q[2];
rz(1.8440856) q[3];
sx q[3];
rz(-1.4617625) q[3];
sx q[3];
rz(1.1980496) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
