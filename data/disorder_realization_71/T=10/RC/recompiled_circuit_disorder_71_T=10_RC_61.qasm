OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.527737) q[0];
sx q[0];
rz(-1.4976488) q[0];
sx q[0];
rz(-2.3117476) q[0];
rz(0.78015503) q[1];
sx q[1];
rz(-2.0766139) q[1];
sx q[1];
rz(2.2652664) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5397545) q[0];
sx q[0];
rz(-1.9617426) q[0];
sx q[0];
rz(1.3215617) q[0];
x q[1];
rz(2.1562188) q[2];
sx q[2];
rz(-1.4101763) q[2];
sx q[2];
rz(-1.0759682) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7180011) q[1];
sx q[1];
rz(-0.33034409) q[1];
sx q[1];
rz(1.6172536) q[1];
x q[2];
rz(1.417744) q[3];
sx q[3];
rz(-1.9733841) q[3];
sx q[3];
rz(3.0229085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.17065389) q[2];
sx q[2];
rz(-1.2761513) q[2];
sx q[2];
rz(2.412964) q[2];
rz(-0.5209926) q[3];
sx q[3];
rz(-0.96121585) q[3];
sx q[3];
rz(-2.9339824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.306863) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(1.1215425) q[0];
rz(-2.8858378) q[1];
sx q[1];
rz(-1.6633727) q[1];
sx q[1];
rz(-2.2671525) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6073608) q[0];
sx q[0];
rz(-0.53593862) q[0];
sx q[0];
rz(0.94632728) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15408709) q[2];
sx q[2];
rz(-1.0006957) q[2];
sx q[2];
rz(-0.64683435) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0622105) q[1];
sx q[1];
rz(-2.1746238) q[1];
sx q[1];
rz(1.0152597) q[1];
rz(-1.3014684) q[3];
sx q[3];
rz(-1.6225617) q[3];
sx q[3];
rz(2.2477637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.740739) q[2];
sx q[2];
rz(-0.48745552) q[2];
sx q[2];
rz(0.43593105) q[2];
rz(0.68108264) q[3];
sx q[3];
rz(-0.77107945) q[3];
sx q[3];
rz(-0.40288231) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9044559) q[0];
sx q[0];
rz(-0.28770068) q[0];
sx q[0];
rz(1.0748192) q[0];
rz(2.3020321) q[1];
sx q[1];
rz(-2.3222175) q[1];
sx q[1];
rz(0.39594617) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2303282) q[0];
sx q[0];
rz(-2.0718144) q[0];
sx q[0];
rz(-2.688174) q[0];
rz(-pi) q[1];
rz(0.01404889) q[2];
sx q[2];
rz(-1.8648246) q[2];
sx q[2];
rz(-0.77169466) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.46519854) q[1];
sx q[1];
rz(-2.4583543) q[1];
sx q[1];
rz(-0.015124358) q[1];
x q[2];
rz(1.5941761) q[3];
sx q[3];
rz(-1.3153207) q[3];
sx q[3];
rz(-1.3083096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5376771) q[2];
sx q[2];
rz(-1.2664412) q[2];
sx q[2];
rz(-2.9023857) q[2];
rz(0.075332969) q[3];
sx q[3];
rz(-1.1139261) q[3];
sx q[3];
rz(-1.3031561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72162119) q[0];
sx q[0];
rz(-2.2890838) q[0];
sx q[0];
rz(0.81992942) q[0];
rz(0.48768249) q[1];
sx q[1];
rz(-0.90355021) q[1];
sx q[1];
rz(-0.23342361) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.061302) q[0];
sx q[0];
rz(-1.6675321) q[0];
sx q[0];
rz(0.060782766) q[0];
rz(-1.8065622) q[2];
sx q[2];
rz(-2.1996017) q[2];
sx q[2];
rz(0.51970383) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.213821) q[1];
sx q[1];
rz(-1.2328316) q[1];
sx q[1];
rz(-0.73422276) q[1];
rz(2.7531284) q[3];
sx q[3];
rz(-0.82287517) q[3];
sx q[3];
rz(0.32271656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0507811) q[2];
sx q[2];
rz(-1.2819141) q[2];
sx q[2];
rz(0.74742571) q[2];
rz(-2.9181972) q[3];
sx q[3];
rz(-0.59745336) q[3];
sx q[3];
rz(2.8994765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6915879) q[0];
sx q[0];
rz(-1.7084028) q[0];
sx q[0];
rz(3.0773556) q[0];
rz(-0.94379395) q[1];
sx q[1];
rz(-2.4143024) q[1];
sx q[1];
rz(-0.76104004) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1431664) q[0];
sx q[0];
rz(-2.879062) q[0];
sx q[0];
rz(-1.6869998) q[0];
rz(-2.9391187) q[2];
sx q[2];
rz(-1.379181) q[2];
sx q[2];
rz(-1.3445878) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3984822) q[1];
sx q[1];
rz(-2.8201436) q[1];
sx q[1];
rz(0.076996315) q[1];
x q[2];
rz(1.8499591) q[3];
sx q[3];
rz(-1.2984635) q[3];
sx q[3];
rz(-2.2945822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3349907) q[2];
sx q[2];
rz(-1.0706173) q[2];
sx q[2];
rz(-0.09207329) q[2];
rz(2.4798685) q[3];
sx q[3];
rz(-2.3501553) q[3];
sx q[3];
rz(-1.7592643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6732366) q[0];
sx q[0];
rz(-1.3156923) q[0];
sx q[0];
rz(-3.1150505) q[0];
rz(-2.2684855) q[1];
sx q[1];
rz(-1.1353506) q[1];
sx q[1];
rz(0.32593265) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3098785) q[0];
sx q[0];
rz(-1.1967812) q[0];
sx q[0];
rz(-2.5783587) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6961156) q[2];
sx q[2];
rz(-1.2947086) q[2];
sx q[2];
rz(2.7518227) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50800486) q[1];
sx q[1];
rz(-1.8691917) q[1];
sx q[1];
rz(-1.7119346) q[1];
rz(1.553922) q[3];
sx q[3];
rz(-1.4697187) q[3];
sx q[3];
rz(0.24731393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59763336) q[2];
sx q[2];
rz(-1.8240857) q[2];
sx q[2];
rz(2.4874172) q[2];
rz(1.7116961) q[3];
sx q[3];
rz(-1.9734029) q[3];
sx q[3];
rz(-0.54106075) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27286801) q[0];
sx q[0];
rz(-1.6690212) q[0];
sx q[0];
rz(0.72189271) q[0];
rz(1.7294653) q[1];
sx q[1];
rz(-2.3528603) q[1];
sx q[1];
rz(0.11925764) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0605474) q[0];
sx q[0];
rz(-1.8348872) q[0];
sx q[0];
rz(-2.413785) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1642045) q[2];
sx q[2];
rz(-0.19492976) q[2];
sx q[2];
rz(0.56602636) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6357928) q[1];
sx q[1];
rz(-1.074082) q[1];
sx q[1];
rz(-2.990681) q[1];
rz(-pi) q[2];
rz(1.9414385) q[3];
sx q[3];
rz(-0.89768411) q[3];
sx q[3];
rz(0.17236575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6922336) q[2];
sx q[2];
rz(-1.8211775) q[2];
sx q[2];
rz(-1.7822441) q[2];
rz(-2.3826777) q[3];
sx q[3];
rz(-0.24154285) q[3];
sx q[3];
rz(-0.53708491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83157241) q[0];
sx q[0];
rz(-2.4825403) q[0];
sx q[0];
rz(1.0634364) q[0];
rz(0.27451441) q[1];
sx q[1];
rz(-1.9332705) q[1];
sx q[1];
rz(-0.88561052) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0622334) q[0];
sx q[0];
rz(-2.057312) q[0];
sx q[0];
rz(1.3079206) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1259414) q[2];
sx q[2];
rz(-0.99132292) q[2];
sx q[2];
rz(1.5660812) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.209621) q[1];
sx q[1];
rz(-1.069427) q[1];
sx q[1];
rz(0.45100905) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86730154) q[3];
sx q[3];
rz(-0.9356381) q[3];
sx q[3];
rz(-0.22658539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.72835913) q[2];
sx q[2];
rz(-0.76247549) q[2];
sx q[2];
rz(-1.1317066) q[2];
rz(2.0570095) q[3];
sx q[3];
rz(-1.0794493) q[3];
sx q[3];
rz(-1.926698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8383012) q[0];
sx q[0];
rz(-1.6475995) q[0];
sx q[0];
rz(-1.0820748) q[0];
rz(-1.2754296) q[1];
sx q[1];
rz(-2.137303) q[1];
sx q[1];
rz(-1.1358322) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4275883) q[0];
sx q[0];
rz(-1.4052608) q[0];
sx q[0];
rz(-1.0323314) q[0];
rz(3.0938221) q[2];
sx q[2];
rz(-0.43521817) q[2];
sx q[2];
rz(1.6795295) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.218704) q[1];
sx q[1];
rz(-1.6195546) q[1];
sx q[1];
rz(2.7152275) q[1];
rz(2.7157852) q[3];
sx q[3];
rz(-2.1771181) q[3];
sx q[3];
rz(3.1283034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11848005) q[2];
sx q[2];
rz(-1.8818405) q[2];
sx q[2];
rz(-2.2686968) q[2];
rz(-2.2980799) q[3];
sx q[3];
rz(-0.25965634) q[3];
sx q[3];
rz(0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89235598) q[0];
sx q[0];
rz(-1.7974239) q[0];
sx q[0];
rz(-1.8632442) q[0];
rz(-1.0247914) q[1];
sx q[1];
rz(-2.0147851) q[1];
sx q[1];
rz(-1.9445673) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68162936) q[0];
sx q[0];
rz(-0.75461331) q[0];
sx q[0];
rz(2.2370536) q[0];
rz(-1.808666) q[2];
sx q[2];
rz(-1.4537721) q[2];
sx q[2];
rz(-0.34840096) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.69233209) q[1];
sx q[1];
rz(-2.3093866) q[1];
sx q[1];
rz(2.2663121) q[1];
x q[2];
rz(0.50935575) q[3];
sx q[3];
rz(-1.9037814) q[3];
sx q[3];
rz(1.2437362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8355576) q[2];
sx q[2];
rz(-0.71283895) q[2];
sx q[2];
rz(2.3416134) q[2];
rz(1.1768613) q[3];
sx q[3];
rz(-1.7088944) q[3];
sx q[3];
rz(-0.95361382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4326614) q[0];
sx q[0];
rz(-0.14871696) q[0];
sx q[0];
rz(0.8014252) q[0];
rz(2.6196383) q[1];
sx q[1];
rz(-2.3028761) q[1];
sx q[1];
rz(0.16470673) q[1];
rz(-0.37408806) q[2];
sx q[2];
rz(-1.2384407) q[2];
sx q[2];
rz(0.25638914) q[2];
rz(2.5064777) q[3];
sx q[3];
rz(-1.0233581) q[3];
sx q[3];
rz(-0.81627853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
