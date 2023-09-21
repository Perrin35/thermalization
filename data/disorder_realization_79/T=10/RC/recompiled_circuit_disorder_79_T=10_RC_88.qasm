OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(-2.2731279) q[0];
sx q[0];
rz(-2.948569) q[0];
rz(-1.9999737) q[1];
sx q[1];
rz(3.5715754) q[1];
sx q[1];
rz(6.9663098) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2642759) q[0];
sx q[0];
rz(-1.5125934) q[0];
sx q[0];
rz(-0.067406128) q[0];
x q[1];
rz(-0.97857742) q[2];
sx q[2];
rz(-0.41523146) q[2];
sx q[2];
rz(-0.65939553) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4343623) q[1];
sx q[1];
rz(-1.1303567) q[1];
sx q[1];
rz(-2.4643154) q[1];
rz(2.1300415) q[3];
sx q[3];
rz(-2.4992001) q[3];
sx q[3];
rz(-1.1005644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0156988) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(-2.0430298) q[2];
rz(2.0627608) q[3];
sx q[3];
rz(-2.1964985) q[3];
sx q[3];
rz(1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9803479) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(-0.20794491) q[0];
rz(-0.57693276) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(-1.4651441) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10509051) q[0];
sx q[0];
rz(-0.72421342) q[0];
sx q[0];
rz(3.1378531) q[0];
rz(-pi) q[1];
rz(0.25603489) q[2];
sx q[2];
rz(-1.0962152) q[2];
sx q[2];
rz(-1.7325967) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1190471) q[1];
sx q[1];
rz(-0.76821583) q[1];
sx q[1];
rz(2.4382298) q[1];
rz(-1.2259237) q[3];
sx q[3];
rz(-0.85540918) q[3];
sx q[3];
rz(-1.8638924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0186105) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(-0.1097651) q[2];
rz(-2.5189853) q[3];
sx q[3];
rz(-2.770335) q[3];
sx q[3];
rz(-1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1784172) q[0];
sx q[0];
rz(-0.85997471) q[0];
sx q[0];
rz(-2.6254568) q[0];
rz(-0.57488817) q[1];
sx q[1];
rz(-0.92620414) q[1];
sx q[1];
rz(0.80054545) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0101937) q[0];
sx q[0];
rz(-0.43361615) q[0];
sx q[0];
rz(-1.2244768) q[0];
rz(0.82528798) q[2];
sx q[2];
rz(-0.85916677) q[2];
sx q[2];
rz(-0.7691783) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.264818) q[1];
sx q[1];
rz(-0.41300981) q[1];
sx q[1];
rz(3.1289711) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9242026) q[3];
sx q[3];
rz(-2.4089775) q[3];
sx q[3];
rz(2.4917847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4642554) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(1.7819972) q[2];
rz(2.1740186) q[3];
sx q[3];
rz(-1.273497) q[3];
sx q[3];
rz(-1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99252218) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(0.46491369) q[0];
rz(-2.7930296) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(1.0850614) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1151428) q[0];
sx q[0];
rz(-1.0466252) q[0];
sx q[0];
rz(-1.8379704) q[0];
rz(-1.9374574) q[2];
sx q[2];
rz(-1.8030093) q[2];
sx q[2];
rz(-1.0166849) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9492053) q[1];
sx q[1];
rz(-2.9678223) q[1];
sx q[1];
rz(2.5237571) q[1];
rz(-2.5892341) q[3];
sx q[3];
rz(-2.4529152) q[3];
sx q[3];
rz(2.3104582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1365635) q[2];
sx q[2];
rz(-1.0483142) q[2];
sx q[2];
rz(-2.4604649) q[2];
rz(-2.629771) q[3];
sx q[3];
rz(-0.32326439) q[3];
sx q[3];
rz(-0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8624449) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(1.408668) q[0];
rz(2.7092343) q[1];
sx q[1];
rz(-0.84195781) q[1];
sx q[1];
rz(2.1599105) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22589332) q[0];
sx q[0];
rz(-0.75076538) q[0];
sx q[0];
rz(2.4322926) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1443411) q[2];
sx q[2];
rz(-2.855636) q[2];
sx q[2];
rz(0.85948247) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2770734) q[1];
sx q[1];
rz(-0.48950567) q[1];
sx q[1];
rz(2.8237052) q[1];
x q[2];
rz(1.6793628) q[3];
sx q[3];
rz(-2.541399) q[3];
sx q[3];
rz(1.997228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0354054) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(2.6521818) q[2];
rz(-1.0148467) q[3];
sx q[3];
rz(-2.615052) q[3];
sx q[3];
rz(-1.813252) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4846102) q[0];
sx q[0];
rz(-2.9841612) q[0];
sx q[0];
rz(-0.37242517) q[0];
rz(1.3308446) q[1];
sx q[1];
rz(-1.0354038) q[1];
sx q[1];
rz(-0.16528027) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80711354) q[0];
sx q[0];
rz(-1.4757336) q[0];
sx q[0];
rz(-1.5414184) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68768244) q[2];
sx q[2];
rz(-1.624794) q[2];
sx q[2];
rz(0.24535594) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2051516) q[1];
sx q[1];
rz(-1.9736104) q[1];
sx q[1];
rz(-2.5564297) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9061756) q[3];
sx q[3];
rz(-1.132292) q[3];
sx q[3];
rz(0.32270839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2281987) q[2];
sx q[2];
rz(-1.0281576) q[2];
sx q[2];
rz(1.4432663) q[2];
rz(-0.36744395) q[3];
sx q[3];
rz(-1.2747217) q[3];
sx q[3];
rz(-2.929556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7845602) q[0];
sx q[0];
rz(-1.0914047) q[0];
sx q[0];
rz(2.7923287) q[0];
rz(2.3941984) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(-2.4051037) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7613477) q[0];
sx q[0];
rz(-1.5537098) q[0];
sx q[0];
rz(-1.6197617) q[0];
x q[1];
rz(-0.24918208) q[2];
sx q[2];
rz(-1.2417925) q[2];
sx q[2];
rz(-2.5788384) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6474825) q[1];
sx q[1];
rz(-0.78299114) q[1];
sx q[1];
rz(-2.8857735) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1398846) q[3];
sx q[3];
rz(-1.1093372) q[3];
sx q[3];
rz(1.0678837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.13654576) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(2.6100256) q[2];
rz(2.452204) q[3];
sx q[3];
rz(-1.6770984) q[3];
sx q[3];
rz(1.2954856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-3.0625967) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(1.2980365) q[0];
rz(-2.334306) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(0.92179006) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73496504) q[0];
sx q[0];
rz(-1.7056744) q[0];
sx q[0];
rz(-0.69358967) q[0];
rz(-pi) q[1];
rz(2.3085262) q[2];
sx q[2];
rz(-1.3563915) q[2];
sx q[2];
rz(0.094878541) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31729749) q[1];
sx q[1];
rz(-1.3166787) q[1];
sx q[1];
rz(1.3586678) q[1];
x q[2];
rz(-0.86254085) q[3];
sx q[3];
rz(-1.4338014) q[3];
sx q[3];
rz(2.9555637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2395997) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(1.1716589) q[2];
rz(-1.7840067) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(1.6931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25306025) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(-1.6660447) q[0];
rz(-1.6015923) q[1];
sx q[1];
rz(-1.6801291) q[1];
sx q[1];
rz(-0.67970651) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1595759) q[0];
sx q[0];
rz(-1.7734818) q[0];
sx q[0];
rz(0.52635877) q[0];
rz(-2.26608) q[2];
sx q[2];
rz(-1.4317703) q[2];
sx q[2];
rz(-1.306844) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9437127) q[1];
sx q[1];
rz(-2.3730952) q[1];
sx q[1];
rz(2.6430623) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4623975) q[3];
sx q[3];
rz(-2.0742356) q[3];
sx q[3];
rz(-0.60590832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0150962) q[2];
sx q[2];
rz(-1.482778) q[2];
sx q[2];
rz(0.91040197) q[2];
rz(-0.67534584) q[3];
sx q[3];
rz(-2.2048435) q[3];
sx q[3];
rz(2.2909686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0891721) q[0];
sx q[0];
rz(-1.2344673) q[0];
sx q[0];
rz(-2.4269379) q[0];
rz(0.71406281) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(-1.9649327) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.048621) q[0];
sx q[0];
rz(-1.7341359) q[0];
sx q[0];
rz(3.0924348) q[0];
rz(-pi) q[1];
rz(-0.87877019) q[2];
sx q[2];
rz(-2.3792017) q[2];
sx q[2];
rz(0.15904418) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9709819) q[1];
sx q[1];
rz(-2.5787528) q[1];
sx q[1];
rz(-1.4303722) q[1];
x q[2];
rz(-0.51104607) q[3];
sx q[3];
rz(-1.647445) q[3];
sx q[3];
rz(2.6443036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.24361336) q[2];
sx q[2];
rz(-1.672013) q[2];
sx q[2];
rz(-1.127355) q[2];
rz(2.7838498) q[3];
sx q[3];
rz(-2.2704411) q[3];
sx q[3];
rz(2.0991142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42416278) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(-2.7453616) q[1];
sx q[1];
rz(-0.025066499) q[1];
sx q[1];
rz(0.33096663) q[1];
rz(2.836543) q[2];
sx q[2];
rz(-0.76022824) q[2];
sx q[2];
rz(-0.40679731) q[2];
rz(1.5789923) q[3];
sx q[3];
rz(-1.9215487) q[3];
sx q[3];
rz(2.4196845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
