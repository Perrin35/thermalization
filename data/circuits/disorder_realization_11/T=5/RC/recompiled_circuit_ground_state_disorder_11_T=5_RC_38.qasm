OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.017555822) q[0];
sx q[0];
rz(-2.9889844) q[0];
sx q[0];
rz(0.18327644) q[0];
rz(0.89219379) q[1];
sx q[1];
rz(5.1434864) q[1];
sx q[1];
rz(8.8465717) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59435049) q[0];
sx q[0];
rz(-3.0495803) q[0];
sx q[0];
rz(3.0334183) q[0];
rz(3.1084539) q[2];
sx q[2];
rz(-1.6553179) q[2];
sx q[2];
rz(-0.13730857) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3751341) q[1];
sx q[1];
rz(-2.441808) q[1];
sx q[1];
rz(2.3699939) q[1];
rz(-1.6434604) q[3];
sx q[3];
rz(-2.6748219) q[3];
sx q[3];
rz(2.6456986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9649428) q[2];
sx q[2];
rz(-0.86805934) q[2];
sx q[2];
rz(3.0435666) q[2];
rz(-0.81884223) q[3];
sx q[3];
rz(-3.0374073) q[3];
sx q[3];
rz(-0.63481832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.1408511) q[0];
sx q[0];
rz(-0.31743693) q[0];
sx q[0];
rz(0.66824085) q[0];
rz(-2.5375598) q[1];
sx q[1];
rz(-3.1112473) q[1];
sx q[1];
rz(-0.86413962) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5174422) q[0];
sx q[0];
rz(-1.6675341) q[0];
sx q[0];
rz(-2.5521432) q[0];
rz(-pi) q[1];
rz(-0.89189826) q[2];
sx q[2];
rz(-1.9111655) q[2];
sx q[2];
rz(-0.75770411) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.171794) q[1];
sx q[1];
rz(-1.329833) q[1];
sx q[1];
rz(1.765224) q[1];
rz(2.941468) q[3];
sx q[3];
rz(-1.8299377) q[3];
sx q[3];
rz(-2.5166496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.981367) q[2];
sx q[2];
rz(-1.9619433) q[2];
sx q[2];
rz(-2.3571864) q[2];
rz(1.9085599) q[3];
sx q[3];
rz(-2.8630856) q[3];
sx q[3];
rz(1.2109141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3311555) q[0];
sx q[0];
rz(-1.5657319) q[0];
sx q[0];
rz(-2.1203777) q[0];
rz(-1.5039697) q[1];
sx q[1];
rz(-0.31871381) q[1];
sx q[1];
rz(2.5655897) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5344328) q[0];
sx q[0];
rz(-2.7078621) q[0];
sx q[0];
rz(2.8612616) q[0];
rz(-pi) q[1];
rz(1.0025315) q[2];
sx q[2];
rz(-2.6883548) q[2];
sx q[2];
rz(-1.9788048) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0124023) q[1];
sx q[1];
rz(-0.92817143) q[1];
sx q[1];
rz(-1.9985649) q[1];
rz(1.7873476) q[3];
sx q[3];
rz(-1.8365905) q[3];
sx q[3];
rz(-2.1590421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7716498) q[2];
sx q[2];
rz(-0.45054951) q[2];
sx q[2];
rz(3.0453299) q[2];
rz(2.7104968) q[3];
sx q[3];
rz(-2.1754706) q[3];
sx q[3];
rz(-0.14804429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3189321) q[0];
sx q[0];
rz(-3.0156101) q[0];
sx q[0];
rz(-0.2952964) q[0];
rz(0.88116208) q[1];
sx q[1];
rz(-2.9144139) q[1];
sx q[1];
rz(-2.6491902) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2797564) q[0];
sx q[0];
rz(-1.4040213) q[0];
sx q[0];
rz(1.357973) q[0];
x q[1];
rz(-1.0855851) q[2];
sx q[2];
rz(-3.036068) q[2];
sx q[2];
rz(-0.72805041) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7192235) q[1];
sx q[1];
rz(-1.2417267) q[1];
sx q[1];
rz(-0.13767875) q[1];
rz(-2.4851068) q[3];
sx q[3];
rz(-0.53561775) q[3];
sx q[3];
rz(-1.519062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.015335036) q[2];
sx q[2];
rz(-2.8110562) q[2];
sx q[2];
rz(2.7988561) q[2];
rz(2.1526509) q[3];
sx q[3];
rz(-1.9804695) q[3];
sx q[3];
rz(-0.69800085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2547176) q[0];
sx q[0];
rz(-2.8509792) q[0];
sx q[0];
rz(1.2683723) q[0];
rz(0.36240029) q[1];
sx q[1];
rz(-2.2762894) q[1];
sx q[1];
rz(1.5579582) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11626205) q[0];
sx q[0];
rz(-1.1209348) q[0];
sx q[0];
rz(2.5309358) q[0];
x q[1];
rz(-0.98619104) q[2];
sx q[2];
rz(-2.8132943) q[2];
sx q[2];
rz(2.3165348) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0456699) q[1];
sx q[1];
rz(-1.5858923) q[1];
sx q[1];
rz(-2.5445055) q[1];
rz(-0.5823808) q[3];
sx q[3];
rz(-0.60027585) q[3];
sx q[3];
rz(1.511661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.6543988) q[2];
sx q[2];
rz(-1.0728269) q[2];
sx q[2];
rz(2.2844592) q[2];
rz(2.1060139) q[3];
sx q[3];
rz(-2.3643957) q[3];
sx q[3];
rz(-2.5100759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71984464) q[0];
sx q[0];
rz(-0.48374614) q[0];
sx q[0];
rz(-2.8068722) q[0];
rz(2.1597629) q[1];
sx q[1];
rz(-1.5900759) q[1];
sx q[1];
rz(-0.88012153) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80106884) q[0];
sx q[0];
rz(-0.33310091) q[0];
sx q[0];
rz(1.0794425) q[0];
x q[1];
rz(2.1507929) q[2];
sx q[2];
rz(-2.8282165) q[2];
sx q[2];
rz(2.6277781) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6374911) q[1];
sx q[1];
rz(-1.7107297) q[1];
sx q[1];
rz(-1.0969694) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7699748) q[3];
sx q[3];
rz(-1.1303969) q[3];
sx q[3];
rz(-0.49345106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6191787) q[2];
sx q[2];
rz(-2.652707) q[2];
sx q[2];
rz(1.6467113) q[2];
rz(1.304168) q[3];
sx q[3];
rz(-2.2924278) q[3];
sx q[3];
rz(-2.9554534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78231597) q[0];
sx q[0];
rz(-2.9501811) q[0];
sx q[0];
rz(0.070847832) q[0];
rz(-2.518636) q[1];
sx q[1];
rz(-0.40095913) q[1];
sx q[1];
rz(-0.43359044) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98862749) q[0];
sx q[0];
rz(-2.835243) q[0];
sx q[0];
rz(3.0393837) q[0];
rz(-1.5135399) q[2];
sx q[2];
rz(-1.861915) q[2];
sx q[2];
rz(1.2886485) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5035491) q[1];
sx q[1];
rz(-1.8942327) q[1];
sx q[1];
rz(1.4694609) q[1];
rz(0.23728063) q[3];
sx q[3];
rz(-10/(7*pi)) q[3];
sx q[3];
rz(0.90516289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.6314466) q[2];
sx q[2];
rz(-0.53312174) q[2];
sx q[2];
rz(-0.637429) q[2];
rz(-2.0006477) q[3];
sx q[3];
rz(-0.44064042) q[3];
sx q[3];
rz(0.91658896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3797853) q[0];
sx q[0];
rz(-2.8768235) q[0];
sx q[0];
rz(-2.8763212) q[0];
rz(0.028566407) q[1];
sx q[1];
rz(-2.9539689) q[1];
sx q[1];
rz(-2.2144337) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2664911) q[0];
sx q[0];
rz(-1.3463839) q[0];
sx q[0];
rz(1.5427179) q[0];
rz(2.7716088) q[2];
sx q[2];
rz(-1.4346037) q[2];
sx q[2];
rz(-0.16384928) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4614951) q[1];
sx q[1];
rz(-2.3465183) q[1];
sx q[1];
rz(0.14393385) q[1];
rz(-1.4577436) q[3];
sx q[3];
rz(-1.9917484) q[3];
sx q[3];
rz(2.497626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.11904968) q[2];
sx q[2];
rz(-2.0930585) q[2];
sx q[2];
rz(-1.0724462) q[2];
rz(-0.33506814) q[3];
sx q[3];
rz(-0.11490331) q[3];
sx q[3];
rz(2.9160299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2304614) q[0];
sx q[0];
rz(-1.1536396) q[0];
sx q[0];
rz(-2.352584) q[0];
rz(-0.59793961) q[1];
sx q[1];
rz(-1.1006678) q[1];
sx q[1];
rz(0.71802968) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1624915) q[0];
sx q[0];
rz(-1.4153061) q[0];
sx q[0];
rz(3.0999712) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4903551) q[2];
sx q[2];
rz(-0.1165963) q[2];
sx q[2];
rz(-1.9060022) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2788417) q[1];
sx q[1];
rz(-0.92581823) q[1];
sx q[1];
rz(0.74826333) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17810693) q[3];
sx q[3];
rz(-1.4796039) q[3];
sx q[3];
rz(-0.9739463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.614552) q[2];
sx q[2];
rz(-2.8895832) q[2];
sx q[2];
rz(1.4177812) q[2];
rz(1.1082015) q[3];
sx q[3];
rz(-2.7751444) q[3];
sx q[3];
rz(2.7537856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80477667) q[0];
sx q[0];
rz(-0.62938654) q[0];
sx q[0];
rz(-0.25548536) q[0];
rz(-2.3707223) q[1];
sx q[1];
rz(-1.5522542) q[1];
sx q[1];
rz(1.593387) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1436018) q[0];
sx q[0];
rz(-1.6842457) q[0];
sx q[0];
rz(-1.1313637) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.911627) q[2];
sx q[2];
rz(-0.88092025) q[2];
sx q[2];
rz(-1.752319) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.76308668) q[1];
sx q[1];
rz(-0.47927472) q[1];
sx q[1];
rz(2.5162426) q[1];
rz(0.58695768) q[3];
sx q[3];
rz(-1.2511593) q[3];
sx q[3];
rz(2.2037821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8678681) q[2];
sx q[2];
rz(-0.74213433) q[2];
sx q[2];
rz(1.8871657) q[2];
rz(-2.6015688) q[3];
sx q[3];
rz(-0.050914474) q[3];
sx q[3];
rz(0.77903265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.222432) q[0];
sx q[0];
rz(-1.5615015) q[0];
sx q[0];
rz(-1.1175565) q[0];
rz(-2.9149105) q[1];
sx q[1];
rz(-3.0381028) q[1];
sx q[1];
rz(1.4729952) q[1];
rz(2.1098577) q[2];
sx q[2];
rz(-2.4895) q[2];
sx q[2];
rz(-1.1498462) q[2];
rz(0.72262121) q[3];
sx q[3];
rz(-1.327805) q[3];
sx q[3];
rz(-1.1205825) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
