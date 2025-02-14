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
rz(2.9583162) q[0];
rz(0.89219379) q[1];
sx q[1];
rz(5.1434864) q[1];
sx q[1];
rz(8.8465717) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86872549) q[0];
sx q[0];
rz(-1.5608762) q[0];
sx q[0];
rz(3.0501151) q[0];
rz(-pi) q[1];
rz(0.033138795) q[2];
sx q[2];
rz(-1.4862747) q[2];
sx q[2];
rz(3.0042841) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7664585) q[1];
sx q[1];
rz(-0.69978461) q[1];
sx q[1];
rz(-2.3699939) q[1];
x q[2];
rz(0.036567612) q[3];
sx q[3];
rz(-1.1053549) q[3];
sx q[3];
rz(0.57722604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1766498) q[2];
sx q[2];
rz(-2.2735333) q[2];
sx q[2];
rz(-3.0435666) q[2];
rz(-2.3227504) q[3];
sx q[3];
rz(-3.0374073) q[3];
sx q[3];
rz(0.63481832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1408511) q[0];
sx q[0];
rz(-0.31743693) q[0];
sx q[0];
rz(-0.66824085) q[0];
rz(-0.60403281) q[1];
sx q[1];
rz(-3.1112473) q[1];
sx q[1];
rz(0.86413962) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0111496) q[0];
sx q[0];
rz(-2.1571223) q[0];
sx q[0];
rz(-1.6870113) q[0];
rz(-pi) q[1];
rz(2.2496944) q[2];
sx q[2];
rz(-1.2304272) q[2];
sx q[2];
rz(0.75770411) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.171794) q[1];
sx q[1];
rz(-1.329833) q[1];
sx q[1];
rz(-1.3763687) q[1];
x q[2];
rz(-1.3066176) q[3];
sx q[3];
rz(-1.3774398) q[3];
sx q[3];
rz(0.99778432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.981367) q[2];
sx q[2];
rz(-1.1796494) q[2];
sx q[2];
rz(-2.3571864) q[2];
rz(1.2330327) q[3];
sx q[3];
rz(-0.27850702) q[3];
sx q[3];
rz(-1.9306785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8104372) q[0];
sx q[0];
rz(-1.5657319) q[0];
sx q[0];
rz(-2.1203777) q[0];
rz(-1.637623) q[1];
sx q[1];
rz(-2.8228788) q[1];
sx q[1];
rz(2.5655897) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9224145) q[0];
sx q[0];
rz(-1.6873345) q[0];
sx q[0];
rz(2.7228505) q[0];
rz(-pi) q[1];
rz(-1.9603278) q[2];
sx q[2];
rz(-1.8086872) q[2];
sx q[2];
rz(2.212461) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6627548) q[1];
sx q[1];
rz(-0.75482268) q[1];
sx q[1];
rz(0.50601064) q[1];
x q[2];
rz(1.7873476) q[3];
sx q[3];
rz(-1.8365905) q[3];
sx q[3];
rz(0.98255052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7716498) q[2];
sx q[2];
rz(-0.45054951) q[2];
sx q[2];
rz(0.096262781) q[2];
rz(2.7104968) q[3];
sx q[3];
rz(-2.1754706) q[3];
sx q[3];
rz(-0.14804429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8226606) q[0];
sx q[0];
rz(-0.12598251) q[0];
sx q[0];
rz(0.2952964) q[0];
rz(-0.88116208) q[1];
sx q[1];
rz(-0.22717871) q[1];
sx q[1];
rz(-2.6491902) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7775901) q[0];
sx q[0];
rz(-2.8719898) q[0];
sx q[0];
rz(-2.2436938) q[0];
rz(-pi) q[1];
rz(-2.0560076) q[2];
sx q[2];
rz(-0.10552465) q[2];
sx q[2];
rz(2.4135422) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.82741) q[1];
sx q[1];
rz(-2.7858509) q[1];
sx q[1];
rz(1.9529424) q[1];
rz(-pi) q[2];
rz(-0.6564859) q[3];
sx q[3];
rz(-0.53561775) q[3];
sx q[3];
rz(-1.6225306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.015335036) q[2];
sx q[2];
rz(-2.8110562) q[2];
sx q[2];
rz(0.34273657) q[2];
rz(0.98894173) q[3];
sx q[3];
rz(-1.1611232) q[3];
sx q[3];
rz(2.4435918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88687503) q[0];
sx q[0];
rz(-0.29061341) q[0];
sx q[0];
rz(1.2683723) q[0];
rz(-0.36240029) q[1];
sx q[1];
rz(-2.2762894) q[1];
sx q[1];
rz(-1.5579582) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11626205) q[0];
sx q[0];
rz(-1.1209348) q[0];
sx q[0];
rz(0.6106569) q[0];
rz(1.8475632) q[2];
sx q[2];
rz(-1.74969) q[2];
sx q[2];
rz(-1.8362294) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.09592274) q[1];
sx q[1];
rz(-1.5557003) q[1];
sx q[1];
rz(-0.59708718) q[1];
x q[2];
rz(-2.5592119) q[3];
sx q[3];
rz(-2.5413168) q[3];
sx q[3];
rz(1.511661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.6543988) q[2];
sx q[2];
rz(-2.0687658) q[2];
sx q[2];
rz(-0.85713345) q[2];
rz(2.1060139) q[3];
sx q[3];
rz(-0.77719694) q[3];
sx q[3];
rz(2.5100759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71984464) q[0];
sx q[0];
rz(-2.6578465) q[0];
sx q[0];
rz(-0.33472043) q[0];
rz(2.1597629) q[1];
sx q[1];
rz(-1.5515168) q[1];
sx q[1];
rz(0.88012153) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80106884) q[0];
sx q[0];
rz(-2.8084917) q[0];
sx q[0];
rz(2.0621501) q[0];
x q[1];
rz(2.9658365) q[2];
sx q[2];
rz(-1.309991) q[2];
sx q[2];
rz(-0.089290451) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5041015) q[1];
sx q[1];
rz(-1.7107297) q[1];
sx q[1];
rz(-1.0969694) q[1];
rz(-0.37161785) q[3];
sx q[3];
rz(-2.0111957) q[3];
sx q[3];
rz(-0.49345106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.52241391) q[2];
sx q[2];
rz(-0.48888561) q[2];
sx q[2];
rz(1.4948814) q[2];
rz(1.304168) q[3];
sx q[3];
rz(-0.84916484) q[3];
sx q[3];
rz(-0.18613923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78231597) q[0];
sx q[0];
rz(-0.19141153) q[0];
sx q[0];
rz(-3.0707448) q[0];
rz(2.518636) q[1];
sx q[1];
rz(-2.7406335) q[1];
sx q[1];
rz(2.7080022) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0957903) q[0];
sx q[0];
rz(-1.266097) q[0];
sx q[0];
rz(1.5385344) q[0];
x q[1];
rz(-1.5135399) q[2];
sx q[2];
rz(-1.2796776) q[2];
sx q[2];
rz(-1.2886485) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2411464) q[1];
sx q[1];
rz(-1.4747319) q[1];
sx q[1];
rz(-0.32498951) q[1];
x q[2];
rz(1.456377) q[3];
sx q[3];
rz(-2.0118719) q[3];
sx q[3];
rz(2.499388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.6314466) q[2];
sx q[2];
rz(-2.6084709) q[2];
sx q[2];
rz(2.5041637) q[2];
rz(2.0006477) q[3];
sx q[3];
rz(-2.7009522) q[3];
sx q[3];
rz(-2.2250037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7618074) q[0];
sx q[0];
rz(-2.8768235) q[0];
sx q[0];
rz(2.8763212) q[0];
rz(0.028566407) q[1];
sx q[1];
rz(-0.18762372) q[1];
sx q[1];
rz(2.2144337) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4521479) q[0];
sx q[0];
rz(-1.5434221) q[0];
sx q[0];
rz(-2.9170947) q[0];
rz(-2.7793332) q[2];
sx q[2];
rz(-2.7484244) q[2];
sx q[2];
rz(-2.071381) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2573942) q[1];
sx q[1];
rz(-0.78620877) q[1];
sx q[1];
rz(1.4255852) q[1];
rz(0.24683471) q[3];
sx q[3];
rz(-2.7066019) q[3];
sx q[3];
rz(0.37295476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.022543) q[2];
sx q[2];
rz(-2.0930585) q[2];
sx q[2];
rz(1.0724462) q[2];
rz(0.33506814) q[3];
sx q[3];
rz(-0.11490331) q[3];
sx q[3];
rz(-2.9160299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91113126) q[0];
sx q[0];
rz(-1.9879531) q[0];
sx q[0];
rz(-0.78900868) q[0];
rz(2.543653) q[1];
sx q[1];
rz(-1.1006678) q[1];
sx q[1];
rz(-2.423563) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58524581) q[0];
sx q[0];
rz(-1.5296773) q[0];
sx q[0];
rz(-1.7264191) q[0];
rz(0.0094114509) q[2];
sx q[2];
rz(-1.6870139) q[2];
sx q[2];
rz(1.1546017) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3586872) q[1];
sx q[1];
rz(-0.99596874) q[1];
sx q[1];
rz(2.3692822) q[1];
rz(-pi) q[2];
rz(-2.6651023) q[3];
sx q[3];
rz(-2.9417178) q[3];
sx q[3];
rz(2.0764004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.614552) q[2];
sx q[2];
rz(-2.8895832) q[2];
sx q[2];
rz(1.4177812) q[2];
rz(-1.1082015) q[3];
sx q[3];
rz(-2.7751444) q[3];
sx q[3];
rz(0.38780701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80477667) q[0];
sx q[0];
rz(-2.5122061) q[0];
sx q[0];
rz(-2.8861073) q[0];
rz(2.3707223) q[1];
sx q[1];
rz(-1.5522542) q[1];
sx q[1];
rz(1.5482056) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99799081) q[0];
sx q[0];
rz(-1.6842457) q[0];
sx q[0];
rz(2.0102289) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8403271) q[2];
sx q[2];
rz(-0.72119207) q[2];
sx q[2];
rz(2.1047956) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7640969) q[1];
sx q[1];
rz(-1.2974655) q[1];
sx q[1];
rz(-0.39876858) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9491626) q[3];
sx q[3];
rz(-1.017166) q[3];
sx q[3];
rz(2.7146482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8678681) q[2];
sx q[2];
rz(-2.3994583) q[2];
sx q[2];
rz(-1.254427) q[2];
rz(-0.54002386) q[3];
sx q[3];
rz(-0.050914474) q[3];
sx q[3];
rz(2.36256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-1.9191606) q[0];
sx q[0];
rz(-1.5800911) q[0];
sx q[0];
rz(2.0240361) q[0];
rz(0.22668214) q[1];
sx q[1];
rz(-3.0381028) q[1];
sx q[1];
rz(1.4729952) q[1];
rz(2.1098577) q[2];
sx q[2];
rz(-2.4895) q[2];
sx q[2];
rz(-1.1498462) q[2];
rz(-1.2516102) q[3];
sx q[3];
rz(-2.2678015) q[3];
sx q[3];
rz(0.24116596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
