OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.1240368) q[0];
sx q[0];
rz(-0.15260829) q[0];
sx q[0];
rz(-0.18327644) q[0];
rz(0.89219379) q[1];
sx q[1];
rz(-1.1396989) q[1];
sx q[1];
rz(-0.5782063) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86872549) q[0];
sx q[0];
rz(-1.5807165) q[0];
sx q[0];
rz(-3.0501151) q[0];
rz(-pi) q[1];
x q[1];
rz(0.033138795) q[2];
sx q[2];
rz(-1.4862747) q[2];
sx q[2];
rz(-0.13730857) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.44408724) q[1];
sx q[1];
rz(-2.0365391) q[1];
sx q[1];
rz(0.5429786) q[1];
rz(2.036506) q[3];
sx q[3];
rz(-1.5381201) q[3];
sx q[3];
rz(1.0099883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1766498) q[2];
sx q[2];
rz(-0.86805934) q[2];
sx q[2];
rz(3.0435666) q[2];
rz(-0.81884223) q[3];
sx q[3];
rz(-3.0374073) q[3];
sx q[3];
rz(2.5067743) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0007415) q[0];
sx q[0];
rz(-2.8241557) q[0];
sx q[0];
rz(0.66824085) q[0];
rz(-0.60403281) q[1];
sx q[1];
rz(-3.1112473) q[1];
sx q[1];
rz(-2.277453) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62415045) q[0];
sx q[0];
rz(-1.4740586) q[0];
sx q[0];
rz(0.58944943) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42704849) q[2];
sx q[2];
rz(-2.2042254) q[2];
sx q[2];
rz(1.0761997) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64795463) q[1];
sx q[1];
rz(-1.3820547) q[1];
sx q[1];
rz(-0.24540875) q[1];
rz(-pi) q[2];
rz(2.941468) q[3];
sx q[3];
rz(-1.311655) q[3];
sx q[3];
rz(-0.6249431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.981367) q[2];
sx q[2];
rz(-1.1796494) q[2];
sx q[2];
rz(-2.3571864) q[2];
rz(-1.9085599) q[3];
sx q[3];
rz(-2.8630856) q[3];
sx q[3];
rz(1.9306785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8104372) q[0];
sx q[0];
rz(-1.5657319) q[0];
sx q[0];
rz(-1.021215) q[0];
rz(-1.637623) q[1];
sx q[1];
rz(-0.31871381) q[1];
sx q[1];
rz(-2.5655897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2191781) q[0];
sx q[0];
rz(-1.4542581) q[0];
sx q[0];
rz(-2.7228505) q[0];
rz(0.25635135) q[2];
sx q[2];
rz(-1.9487938) q[2];
sx q[2];
rz(-2.5963654) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4788379) q[1];
sx q[1];
rz(-2.38677) q[1];
sx q[1];
rz(0.50601064) q[1];
rz(2.473439) q[3];
sx q[3];
rz(-0.34121554) q[3];
sx q[3];
rz(1.4618272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7716498) q[2];
sx q[2];
rz(-0.45054951) q[2];
sx q[2];
rz(-3.0453299) q[2];
rz(-2.7104968) q[3];
sx q[3];
rz(-0.96612203) q[3];
sx q[3];
rz(2.9935484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8226606) q[0];
sx q[0];
rz(-3.0156101) q[0];
sx q[0];
rz(0.2952964) q[0];
rz(2.2604306) q[1];
sx q[1];
rz(-0.22717871) q[1];
sx q[1];
rz(0.49240246) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8146956) q[0];
sx q[0];
rz(-1.7806223) q[0];
sx q[0];
rz(-2.9710415) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4773764) q[2];
sx q[2];
rz(-1.619941) q[2];
sx q[2];
rz(-1.8159332) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.82741) q[1];
sx q[1];
rz(-0.35574177) q[1];
sx q[1];
rz(1.9529424) q[1];
rz(-pi) q[2];
rz(-1.9183228) q[3];
sx q[3];
rz(-1.9869958) q[3];
sx q[3];
rz(-0.89194854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1262576) q[2];
sx q[2];
rz(-2.8110562) q[2];
sx q[2];
rz(-2.7988561) q[2];
rz(-2.1526509) q[3];
sx q[3];
rz(-1.1611232) q[3];
sx q[3];
rz(-0.69800085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.2547176) q[0];
sx q[0];
rz(-2.8509792) q[0];
sx q[0];
rz(-1.8732204) q[0];
rz(-2.7791924) q[1];
sx q[1];
rz(-2.2762894) q[1];
sx q[1];
rz(-1.5836345) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1311296) q[0];
sx q[0];
rz(-0.74105019) q[0];
sx q[0];
rz(-0.69990943) q[0];
x q[1];
rz(-2.1554016) q[2];
sx q[2];
rz(-0.3282983) q[2];
sx q[2];
rz(-0.82505783) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4526738) q[1];
sx q[1];
rz(-2.5443379) q[1];
sx q[1];
rz(-3.1147472) q[1];
rz(-pi) q[2];
rz(1.2107047) q[3];
sx q[3];
rz(-1.0795169) q[3];
sx q[3];
rz(-2.1852428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4871939) q[2];
sx q[2];
rz(-2.0687658) q[2];
sx q[2];
rz(2.2844592) q[2];
rz(-2.1060139) q[3];
sx q[3];
rz(-2.3643957) q[3];
sx q[3];
rz(-0.63151675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71984464) q[0];
sx q[0];
rz(-0.48374614) q[0];
sx q[0];
rz(0.33472043) q[0];
rz(-0.98182976) q[1];
sx q[1];
rz(-1.5515168) q[1];
sx q[1];
rz(0.88012153) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9036569) q[0];
sx q[0];
rz(-1.7256883) q[0];
sx q[0];
rz(1.2747033) q[0];
x q[1];
rz(-1.8354957) q[2];
sx q[2];
rz(-1.4010426) q[2];
sx q[2];
rz(1.614326) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6374911) q[1];
sx q[1];
rz(-1.7107297) q[1];
sx q[1];
rz(1.0969694) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7699748) q[3];
sx q[3];
rz(-2.0111957) q[3];
sx q[3];
rz(0.49345106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.52241391) q[2];
sx q[2];
rz(-0.48888561) q[2];
sx q[2];
rz(1.6467113) q[2];
rz(1.8374247) q[3];
sx q[3];
rz(-0.84916484) q[3];
sx q[3];
rz(-2.9554534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.78231597) q[0];
sx q[0];
rz(-2.9501811) q[0];
sx q[0];
rz(-0.070847832) q[0];
rz(0.62295667) q[1];
sx q[1];
rz(-0.40095913) q[1];
sx q[1];
rz(-0.43359044) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6569048) q[0];
sx q[0];
rz(-1.6015717) q[0];
sx q[0];
rz(-2.8367443) q[0];
rz(-2.8500227) q[2];
sx q[2];
rz(-1.6256412) q[2];
sx q[2];
rz(0.26569732) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.63804352) q[1];
sx q[1];
rz(-1.2473599) q[1];
sx q[1];
rz(1.4694609) q[1];
rz(-pi) q[2];
rz(-1.456377) q[3];
sx q[3];
rz(-1.1297207) q[3];
sx q[3];
rz(2.499388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.6314466) q[2];
sx q[2];
rz(-0.53312174) q[2];
sx q[2];
rz(2.5041637) q[2];
rz(-1.140945) q[3];
sx q[3];
rz(-2.7009522) q[3];
sx q[3];
rz(-2.2250037) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7618074) q[0];
sx q[0];
rz(-2.8768235) q[0];
sx q[0];
rz(-0.26527143) q[0];
rz(-0.028566407) q[1];
sx q[1];
rz(-0.18762372) q[1];
sx q[1];
rz(0.92715895) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4521479) q[0];
sx q[0];
rz(-1.5434221) q[0];
sx q[0];
rz(2.9170947) q[0];
x q[1];
rz(-2.7716088) q[2];
sx q[2];
rz(-1.4346037) q[2];
sx q[2];
rz(-2.9777434) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3520364) q[1];
sx q[1];
rz(-1.4682143) q[1];
sx q[1];
rz(2.351715) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42334564) q[3];
sx q[3];
rz(-1.6739427) q[3];
sx q[3];
rz(-2.1684017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.022543) q[2];
sx q[2];
rz(-2.0930585) q[2];
sx q[2];
rz(-1.0724462) q[2];
rz(2.8065245) q[3];
sx q[3];
rz(-0.11490331) q[3];
sx q[3];
rz(2.9160299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2304614) q[0];
sx q[0];
rz(-1.1536396) q[0];
sx q[0];
rz(-0.78900868) q[0];
rz(-0.59793961) q[1];
sx q[1];
rz(-2.0409248) q[1];
sx q[1];
rz(2.423563) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58524581) q[0];
sx q[0];
rz(-1.5296773) q[0];
sx q[0];
rz(1.7264191) q[0];
x q[1];
rz(-1.687019) q[2];
sx q[2];
rz(-1.5801443) q[2];
sx q[2];
rz(2.7264894) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.78290547) q[1];
sx q[1];
rz(-0.99596874) q[1];
sx q[1];
rz(2.3692822) q[1];
x q[2];
rz(-1.6634462) q[3];
sx q[3];
rz(-1.7481553) q[3];
sx q[3];
rz(2.5611344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52704063) q[2];
sx q[2];
rz(-2.8895832) q[2];
sx q[2];
rz(-1.4177812) q[2];
rz(2.0333911) q[3];
sx q[3];
rz(-2.7751444) q[3];
sx q[3];
rz(0.38780701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.80477667) q[0];
sx q[0];
rz(-2.5122061) q[0];
sx q[0];
rz(0.25548536) q[0];
rz(2.3707223) q[1];
sx q[1];
rz(-1.5522542) q[1];
sx q[1];
rz(1.5482056) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3324677) q[0];
sx q[0];
rz(-2.6886779) q[0];
sx q[0];
rz(1.8324773) q[0];
rz(-1.8403271) q[2];
sx q[2];
rz(-0.72119207) q[2];
sx q[2];
rz(1.0367971) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0615427) q[1];
sx q[1];
rz(-1.953974) q[1];
sx q[1];
rz(1.2754759) q[1];
rz(-0.53867619) q[3];
sx q[3];
rz(-0.65924257) q[3];
sx q[3];
rz(-1.074312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2737246) q[2];
sx q[2];
rz(-2.3994583) q[2];
sx q[2];
rz(1.8871657) q[2];
rz(-0.54002386) q[3];
sx q[3];
rz(-3.0906782) q[3];
sx q[3];
rz(-2.36256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9191606) q[0];
sx q[0];
rz(-1.5615015) q[0];
sx q[0];
rz(-1.1175565) q[0];
rz(-0.22668214) q[1];
sx q[1];
rz(-0.10348987) q[1];
sx q[1];
rz(-1.6685974) q[1];
rz(-0.37353362) q[2];
sx q[2];
rz(-2.1185771) q[2];
sx q[2];
rz(2.6368903) q[2];
rz(0.35861438) q[3];
sx q[3];
rz(-0.75533854) q[3];
sx q[3];
rz(0.71661267) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
