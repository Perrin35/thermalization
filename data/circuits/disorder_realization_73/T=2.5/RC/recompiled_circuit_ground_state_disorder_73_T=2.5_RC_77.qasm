OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.29326987) q[0];
sx q[0];
rz(-2.8889416) q[0];
sx q[0];
rz(1.2055612) q[0];
rz(2.0134917) q[1];
sx q[1];
rz(-1.3265346) q[1];
sx q[1];
rz(-2.3958652) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099961258) q[0];
sx q[0];
rz(-1.3486762) q[0];
sx q[0];
rz(3.0801386) q[0];
rz(-3.1118561) q[2];
sx q[2];
rz(-0.22448891) q[2];
sx q[2];
rz(0.19408524) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8516525) q[1];
sx q[1];
rz(-1.9059685) q[1];
sx q[1];
rz(-0.22214684) q[1];
rz(-pi) q[2];
rz(-1.9609591) q[3];
sx q[3];
rz(-0.50658617) q[3];
sx q[3];
rz(-0.49662874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2241609) q[2];
sx q[2];
rz(-1.5459583) q[2];
sx q[2];
rz(1.6708299) q[2];
rz(-1.6338232) q[3];
sx q[3];
rz(-0.016851146) q[3];
sx q[3];
rz(-0.9233709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725752) q[0];
sx q[0];
rz(-1.9439531) q[0];
sx q[0];
rz(1.5731328) q[0];
rz(-0.16954999) q[1];
sx q[1];
rz(-3.0270271) q[1];
sx q[1];
rz(0.13470185) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41551521) q[0];
sx q[0];
rz(-1.9163791) q[0];
sx q[0];
rz(-0.4536566) q[0];
x q[1];
rz(1.3943761) q[2];
sx q[2];
rz(-3.0718832) q[2];
sx q[2];
rz(-1.4331872) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1253648) q[1];
sx q[1];
rz(-1.1332336) q[1];
sx q[1];
rz(-1.9047649) q[1];
x q[2];
rz(-3.0613741) q[3];
sx q[3];
rz(-1.748198) q[3];
sx q[3];
rz(-2.3631848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1255101) q[2];
sx q[2];
rz(-1.6566015) q[2];
sx q[2];
rz(-2.9888195) q[2];
rz(-1.7759391) q[3];
sx q[3];
rz(-3.1047265) q[3];
sx q[3];
rz(-0.16807817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.633054) q[0];
sx q[0];
rz(-2.3336053) q[0];
sx q[0];
rz(-2.659944) q[0];
rz(-2.956849) q[1];
sx q[1];
rz(-1.7748723) q[1];
sx q[1];
rz(2.1462323) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99670519) q[0];
sx q[0];
rz(-1.2249682) q[0];
sx q[0];
rz(-2.8976827) q[0];
x q[1];
rz(-1.6312509) q[2];
sx q[2];
rz(-1.5222856) q[2];
sx q[2];
rz(-1.8497576) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0479483) q[1];
sx q[1];
rz(-1.8694832) q[1];
sx q[1];
rz(0.68409749) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97484346) q[3];
sx q[3];
rz(-1.5454834) q[3];
sx q[3];
rz(1.5060177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.92622009) q[2];
sx q[2];
rz(-3.0871349) q[2];
sx q[2];
rz(-0.064662956) q[2];
rz(1.0489382) q[3];
sx q[3];
rz(-3.1147396) q[3];
sx q[3];
rz(1.2803199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8365086) q[0];
sx q[0];
rz(-3.0236112) q[0];
sx q[0];
rz(-2.3205561) q[0];
rz(-3.0631284) q[1];
sx q[1];
rz(-1.4469701) q[1];
sx q[1];
rz(-0.99748126) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067046384) q[0];
sx q[0];
rz(-2.1830478) q[0];
sx q[0];
rz(-1.9590098) q[0];
x q[1];
rz(0.036254739) q[2];
sx q[2];
rz(-1.6299575) q[2];
sx q[2];
rz(1.7922557) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9033244) q[1];
sx q[1];
rz(-2.0966171) q[1];
sx q[1];
rz(0.30491288) q[1];
x q[2];
rz(-0.11028744) q[3];
sx q[3];
rz(-0.48088851) q[3];
sx q[3];
rz(0.65910027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0923826) q[2];
sx q[2];
rz(-0.02928484) q[2];
sx q[2];
rz(1.1537665) q[2];
rz(-2.9554101) q[3];
sx q[3];
rz(-0.081705339) q[3];
sx q[3];
rz(2.8103099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1368197) q[0];
sx q[0];
rz(-0.81757075) q[0];
sx q[0];
rz(1.9983043) q[0];
rz(1.1413057) q[1];
sx q[1];
rz(-0.80991304) q[1];
sx q[1];
rz(0.51796651) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9310936) q[0];
sx q[0];
rz(-2.2552975) q[0];
sx q[0];
rz(0.71523358) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1406691) q[2];
sx q[2];
rz(-1.5858272) q[2];
sx q[2];
rz(1.817944) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8570366) q[1];
sx q[1];
rz(-1.9003881) q[1];
sx q[1];
rz(2.8524141) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6916396) q[3];
sx q[3];
rz(-0.3287238) q[3];
sx q[3];
rz(0.69531073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3190069) q[2];
sx q[2];
rz(-2.1876882) q[2];
sx q[2];
rz(0.32773584) q[2];
rz(-1.94708) q[3];
sx q[3];
rz(-3.0142398) q[3];
sx q[3];
rz(-2.2849042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1389393) q[0];
sx q[0];
rz(-2.8728573) q[0];
sx q[0];
rz(-0.56185454) q[0];
rz(1.6506763) q[1];
sx q[1];
rz(-1.6249388) q[1];
sx q[1];
rz(0.098310016) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5033103) q[0];
sx q[0];
rz(-1.1489963) q[0];
sx q[0];
rz(3.1027334) q[0];
rz(-pi) q[1];
rz(0.00064050015) q[2];
sx q[2];
rz(-1.5709236) q[2];
sx q[2];
rz(1.887111) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1001818) q[1];
sx q[1];
rz(-1.5242002) q[1];
sx q[1];
rz(-0.46943922) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0344072) q[3];
sx q[3];
rz(-2.0258198) q[3];
sx q[3];
rz(0.95038271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.085122846) q[2];
sx q[2];
rz(-0.19530185) q[2];
sx q[2];
rz(1.0657715) q[2];
rz(0.39984518) q[3];
sx q[3];
rz(-0.53374922) q[3];
sx q[3];
rz(1.8736418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9999076) q[0];
sx q[0];
rz(-3.0596924) q[0];
sx q[0];
rz(1.6837233) q[0];
rz(-1.1204002) q[1];
sx q[1];
rz(-0.1381865) q[1];
sx q[1];
rz(-2.8021326) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65050478) q[0];
sx q[0];
rz(-1.5575214) q[0];
sx q[0];
rz(-0.076939452) q[0];
rz(-pi) q[1];
rz(3.1192084) q[2];
sx q[2];
rz(-2.5867992) q[2];
sx q[2];
rz(0.024352976) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.124954) q[1];
sx q[1];
rz(-1.5691688) q[1];
sx q[1];
rz(-1.485199) q[1];
rz(1.6704329) q[3];
sx q[3];
rz(-1.2335586) q[3];
sx q[3];
rz(-1.8189614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3642984) q[2];
sx q[2];
rz(-3.1396781) q[2];
sx q[2];
rz(2.778229) q[2];
rz(-2.0511138) q[3];
sx q[3];
rz(-0.57991475) q[3];
sx q[3];
rz(2.0232078) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61984396) q[0];
sx q[0];
rz(-0.32829568) q[0];
sx q[0];
rz(-2.2186665) q[0];
rz(1.482831) q[1];
sx q[1];
rz(-0.62983477) q[1];
sx q[1];
rz(-0.0042075687) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9516966) q[0];
sx q[0];
rz(-1.0605264) q[0];
sx q[0];
rz(2.4416591) q[0];
rz(-pi) q[1];
rz(-1.5571874) q[2];
sx q[2];
rz(-0.41092349) q[2];
sx q[2];
rz(0.027337242) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.15531047) q[1];
sx q[1];
rz(-2.661507) q[1];
sx q[1];
rz(1.7106777) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61242731) q[3];
sx q[3];
rz(-1.6086626) q[3];
sx q[3];
rz(-1.7671536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5795472) q[2];
sx q[2];
rz(-1.6041218) q[2];
sx q[2];
rz(-1.9407678) q[2];
rz(-1.3935401) q[3];
sx q[3];
rz(-3.1378919) q[3];
sx q[3];
rz(-0.70011955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8414108) q[0];
sx q[0];
rz(-0.44898471) q[0];
sx q[0];
rz(1.9218943) q[0];
rz(1.8118743) q[1];
sx q[1];
rz(-1.1390353) q[1];
sx q[1];
rz(-2.9776998) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4224098) q[0];
sx q[0];
rz(-1.9253732) q[0];
sx q[0];
rz(0.18300458) q[0];
rz(-1.8437949) q[2];
sx q[2];
rz(-2.1381209) q[2];
sx q[2];
rz(1.8253872) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1023952) q[1];
sx q[1];
rz(-1.6886687) q[1];
sx q[1];
rz(-3.0815691) q[1];
x q[2];
rz(0.74045351) q[3];
sx q[3];
rz(-1.0426842) q[3];
sx q[3];
rz(-1.5453366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71358744) q[2];
sx q[2];
rz(-3.0516477) q[2];
sx q[2];
rz(-0.62291992) q[2];
rz(-3.0981787) q[3];
sx q[3];
rz(-2.2446938) q[3];
sx q[3];
rz(2.3626732) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49122214) q[0];
sx q[0];
rz(-0.0064042052) q[0];
sx q[0];
rz(-0.48625913) q[0];
rz(0.69475118) q[1];
sx q[1];
rz(-0.34725747) q[1];
sx q[1];
rz(-0.4756701) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18494517) q[0];
sx q[0];
rz(-1.534919) q[0];
sx q[0];
rz(0.033680276) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6108405) q[2];
sx q[2];
rz(-0.70654987) q[2];
sx q[2];
rz(1.5020811) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.116058) q[1];
sx q[1];
rz(-2.5339911) q[1];
sx q[1];
rz(3.1302384) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3856412) q[3];
sx q[3];
rz(-1.2655711) q[3];
sx q[3];
rz(-0.023630649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4645369) q[2];
sx q[2];
rz(-0.060364351) q[2];
sx q[2];
rz(-1.2976868) q[2];
rz(-1.248598) q[3];
sx q[3];
rz(-2.5864351) q[3];
sx q[3];
rz(0.36778522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0155335) q[0];
sx q[0];
rz(-1.5801237) q[0];
sx q[0];
rz(1.7042241) q[0];
rz(0.88232782) q[1];
sx q[1];
rz(-3.075141) q[1];
sx q[1];
rz(-2.3022423) q[1];
rz(2.7693314) q[2];
sx q[2];
rz(-0.098975565) q[2];
sx q[2];
rz(3.0437058) q[2];
rz(0.23918693) q[3];
sx q[3];
rz(-1.5755972) q[3];
sx q[3];
rz(1.5920873) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
