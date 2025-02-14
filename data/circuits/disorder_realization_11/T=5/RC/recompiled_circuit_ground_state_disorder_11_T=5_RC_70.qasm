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
rz(3.2942009) q[0];
sx q[0];
rz(9.6080544) q[0];
rz(0.89219379) q[1];
sx q[1];
rz(-1.1396989) q[1];
sx q[1];
rz(-0.5782063) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70298083) q[0];
sx q[0];
rz(-1.6622694) q[0];
sx q[0];
rz(1.5807581) q[0];
rz(-pi) q[1];
rz(1.9435801) q[2];
sx q[2];
rz(-0.090771505) q[2];
sx q[2];
rz(-2.6300989) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2795434) q[1];
sx q[1];
rz(-1.0909362) q[1];
sx q[1];
rz(-1.0399489) q[1];
rz(1.6434604) q[3];
sx q[3];
rz(-0.46677073) q[3];
sx q[3];
rz(-0.49589402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9649428) q[2];
sx q[2];
rz(-0.86805934) q[2];
sx q[2];
rz(-3.0435666) q[2];
rz(-0.81884223) q[3];
sx q[3];
rz(-0.10418532) q[3];
sx q[3];
rz(-2.5067743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0007415) q[0];
sx q[0];
rz(-2.8241557) q[0];
sx q[0];
rz(-2.4733518) q[0];
rz(-2.5375598) q[1];
sx q[1];
rz(-0.030345358) q[1];
sx q[1];
rz(0.86413962) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5174422) q[0];
sx q[0];
rz(-1.4740586) q[0];
sx q[0];
rz(0.58944943) q[0];
rz(-pi) q[1];
rz(-0.89189826) q[2];
sx q[2];
rz(-1.9111655) q[2];
sx q[2];
rz(-0.75770411) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.64795463) q[1];
sx q[1];
rz(-1.3820547) q[1];
sx q[1];
rz(-2.8961839) q[1];
rz(-pi) q[2];
rz(2.941468) q[3];
sx q[3];
rz(-1.311655) q[3];
sx q[3];
rz(-0.6249431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.16022564) q[2];
sx q[2];
rz(-1.1796494) q[2];
sx q[2];
rz(-2.3571864) q[2];
rz(-1.9085599) q[3];
sx q[3];
rz(-0.27850702) q[3];
sx q[3];
rz(-1.9306785) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3311555) q[0];
sx q[0];
rz(-1.5657319) q[0];
sx q[0];
rz(1.021215) q[0];
rz(-1.5039697) q[1];
sx q[1];
rz(-0.31871381) q[1];
sx q[1];
rz(2.5655897) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60715985) q[0];
sx q[0];
rz(-2.7078621) q[0];
sx q[0];
rz(2.8612616) q[0];
x q[1];
rz(-1.1812649) q[2];
sx q[2];
rz(-1.8086872) q[2];
sx q[2];
rz(-2.212461) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6627548) q[1];
sx q[1];
rz(-2.38677) q[1];
sx q[1];
rz(-0.50601064) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8697474) q[3];
sx q[3];
rz(-1.3619641) q[3];
sx q[3];
rz(2.6110716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7716498) q[2];
sx q[2];
rz(-0.45054951) q[2];
sx q[2];
rz(3.0453299) q[2];
rz(2.7104968) q[3];
sx q[3];
rz(-2.1754706) q[3];
sx q[3];
rz(2.9935484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3189321) q[0];
sx q[0];
rz(-3.0156101) q[0];
sx q[0];
rz(2.8462963) q[0];
rz(2.2604306) q[1];
sx q[1];
rz(-2.9144139) q[1];
sx q[1];
rz(2.6491902) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7775901) q[0];
sx q[0];
rz(-0.26960281) q[0];
sx q[0];
rz(-0.89789884) q[0];
rz(-pi) q[1];
rz(3.0922331) q[2];
sx q[2];
rz(-1.6641031) q[2];
sx q[2];
rz(2.9010584) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7192235) q[1];
sx q[1];
rz(-1.8998659) q[1];
sx q[1];
rz(0.13767875) q[1];
rz(-pi) q[2];
rz(-0.6564859) q[3];
sx q[3];
rz(-0.53561775) q[3];
sx q[3];
rz(1.519062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.015335036) q[2];
sx q[2];
rz(-0.33053645) q[2];
sx q[2];
rz(2.7988561) q[2];
rz(0.98894173) q[3];
sx q[3];
rz(-1.1611232) q[3];
sx q[3];
rz(2.4435918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-0.86530322) q[1];
sx q[1];
rz(1.5836345) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9824955) q[0];
sx q[0];
rz(-1.0282059) q[0];
sx q[0];
rz(1.0382023) q[0];
rz(-pi) q[1];
rz(-2.1554016) q[2];
sx q[2];
rz(-2.8132943) q[2];
sx q[2];
rz(0.82505783) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.09592274) q[1];
sx q[1];
rz(-1.5557003) q[1];
sx q[1];
rz(-2.5445055) q[1];
rz(0.5823808) q[3];
sx q[3];
rz(-0.60027585) q[3];
sx q[3];
rz(1.6299316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4871939) q[2];
sx q[2];
rz(-1.0728269) q[2];
sx q[2];
rz(-0.85713345) q[2];
rz(-2.1060139) q[3];
sx q[3];
rz(-2.3643957) q[3];
sx q[3];
rz(2.5100759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.421748) q[0];
sx q[0];
rz(-2.6578465) q[0];
sx q[0];
rz(-0.33472043) q[0];
rz(0.98182976) q[1];
sx q[1];
rz(-1.5515168) q[1];
sx q[1];
rz(2.2614711) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80106884) q[0];
sx q[0];
rz(-0.33310091) q[0];
sx q[0];
rz(-1.0794425) q[0];
x q[1];
rz(2.1507929) q[2];
sx q[2];
rz(-2.8282165) q[2];
sx q[2];
rz(-0.51381451) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.33226686) q[1];
sx q[1];
rz(-2.6490445) q[1];
sx q[1];
rz(-1.2713856) q[1];
rz(-pi) q[2];
rz(-0.91428925) q[3];
sx q[3];
rz(-0.56826545) q[3];
sx q[3];
rz(-1.9075346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.52241391) q[2];
sx q[2];
rz(-0.48888561) q[2];
sx q[2];
rz(1.4948814) q[2];
rz(1.8374247) q[3];
sx q[3];
rz(-0.84916484) q[3];
sx q[3];
rz(0.18613923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78231597) q[0];
sx q[0];
rz(-2.9501811) q[0];
sx q[0];
rz(-0.070847832) q[0];
rz(-2.518636) q[1];
sx q[1];
rz(-2.7406335) q[1];
sx q[1];
rz(-2.7080022) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0458023) q[0];
sx q[0];
rz(-1.266097) q[0];
sx q[0];
rz(-1.6030583) q[0];
rz(-pi) q[1];
rz(-1.6280528) q[2];
sx q[2];
rz(-1.861915) q[2];
sx q[2];
rz(1.8529441) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9476871) q[1];
sx q[1];
rz(-0.33840678) q[1];
sx q[1];
rz(2.848495) q[1];
rz(-pi) q[2];
x q[2];
rz(1.456377) q[3];
sx q[3];
rz(-1.1297207) q[3];
sx q[3];
rz(0.64220465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.6314466) q[2];
sx q[2];
rz(-0.53312174) q[2];
sx q[2];
rz(2.5041637) q[2];
rz(-2.0006477) q[3];
sx q[3];
rz(-2.7009522) q[3];
sx q[3];
rz(2.2250037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7618074) q[0];
sx q[0];
rz(-0.26476911) q[0];
sx q[0];
rz(-2.8763212) q[0];
rz(3.1130262) q[1];
sx q[1];
rz(-0.18762372) q[1];
sx q[1];
rz(0.92715895) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68944478) q[0];
sx q[0];
rz(-1.5981705) q[0];
sx q[0];
rz(2.9170947) q[0];
rz(-pi) q[1];
rz(0.36998387) q[2];
sx q[2];
rz(-1.4346037) q[2];
sx q[2];
rz(-2.9777434) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4614951) q[1];
sx q[1];
rz(-0.79507438) q[1];
sx q[1];
rz(0.14393385) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6838491) q[3];
sx q[3];
rz(-1.9917484) q[3];
sx q[3];
rz(-0.64396665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.022543) q[2];
sx q[2];
rz(-2.0930585) q[2];
sx q[2];
rz(-1.0724462) q[2];
rz(-0.33506814) q[3];
sx q[3];
rz(-0.11490331) q[3];
sx q[3];
rz(-0.2255628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2304614) q[0];
sx q[0];
rz(-1.1536396) q[0];
sx q[0];
rz(0.78900868) q[0];
rz(-2.543653) q[1];
sx q[1];
rz(-2.0409248) q[1];
sx q[1];
rz(-2.423563) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1624915) q[0];
sx q[0];
rz(-1.7262865) q[0];
sx q[0];
rz(-3.0999712) q[0];
rz(-pi) q[1];
rz(-1.6512376) q[2];
sx q[2];
rz(-3.0249964) q[2];
sx q[2];
rz(1.2355905) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2749697) q[1];
sx q[1];
rz(-2.1964745) q[1];
sx q[1];
rz(-2.3060206) q[1];
rz(2.6651023) q[3];
sx q[3];
rz(-0.19987488) q[3];
sx q[3];
rz(-1.0651922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52704063) q[2];
sx q[2];
rz(-0.25200945) q[2];
sx q[2];
rz(1.7238114) q[2];
rz(2.0333911) q[3];
sx q[3];
rz(-0.36644822) q[3];
sx q[3];
rz(-0.38780701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80477667) q[0];
sx q[0];
rz(-2.5122061) q[0];
sx q[0];
rz(-2.8861073) q[0];
rz(0.77087036) q[1];
sx q[1];
rz(-1.5893385) q[1];
sx q[1];
rz(1.5482056) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6219538) q[0];
sx q[0];
rz(-2.0072091) q[0];
sx q[0];
rz(-0.12524097) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86779197) q[2];
sx q[2];
rz(-1.3940587) q[2];
sx q[2];
rz(2.8121626) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.378506) q[1];
sx q[1];
rz(-0.47927472) q[1];
sx q[1];
rz(-2.5162426) q[1];
rz(-2.554635) q[3];
sx q[3];
rz(-1.8904333) q[3];
sx q[3];
rz(0.93781059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8678681) q[2];
sx q[2];
rz(-2.3994583) q[2];
sx q[2];
rz(1.8871657) q[2];
rz(2.6015688) q[3];
sx q[3];
rz(-0.050914474) q[3];
sx q[3];
rz(2.36256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.222432) q[0];
sx q[0];
rz(-1.5800911) q[0];
sx q[0];
rz(2.0240361) q[0];
rz(2.9149105) q[1];
sx q[1];
rz(-0.10348987) q[1];
sx q[1];
rz(-1.6685974) q[1];
rz(1.031735) q[2];
sx q[2];
rz(-0.65209263) q[2];
sx q[2];
rz(1.9917464) q[2];
rz(2.7829783) q[3];
sx q[3];
rz(-2.3862541) q[3];
sx q[3];
rz(-2.42498) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
