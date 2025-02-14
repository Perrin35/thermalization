OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.26389709) q[0];
sx q[0];
rz(2.1030302) q[0];
sx q[0];
rz(9.8197688) q[0];
rz(0.78139961) q[1];
sx q[1];
rz(7.6548792) q[1];
sx q[1];
rz(10.21488) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.794644) q[0];
sx q[0];
rz(-0.99604368) q[0];
sx q[0];
rz(-1.5331506) q[0];
rz(2.6955881) q[2];
sx q[2];
rz(-0.82740192) q[2];
sx q[2];
rz(-2.3349886) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.61617) q[1];
sx q[1];
rz(-1.425263) q[1];
sx q[1];
rz(2.1110181) q[1];
rz(-pi) q[2];
rz(-1.9856307) q[3];
sx q[3];
rz(-1.1006163) q[3];
sx q[3];
rz(-3.1347203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.52749741) q[2];
sx q[2];
rz(-3.0588394) q[2];
sx q[2];
rz(-2.5772074) q[2];
rz(0.91974059) q[3];
sx q[3];
rz(-2.4904833) q[3];
sx q[3];
rz(1.9981599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98519242) q[0];
sx q[0];
rz(-0.99034482) q[0];
sx q[0];
rz(1.9957969) q[0];
rz(2.1339553) q[1];
sx q[1];
rz(-0.8877019) q[1];
sx q[1];
rz(0.064858286) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6332209) q[0];
sx q[0];
rz(-2.457058) q[0];
sx q[0];
rz(0.9100911) q[0];
x q[1];
rz(-2.1880661) q[2];
sx q[2];
rz(-2.2489886) q[2];
sx q[2];
rz(1.8770521) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.25003281) q[1];
sx q[1];
rz(-1.7549101) q[1];
sx q[1];
rz(1.5450983) q[1];
rz(-2.6363144) q[3];
sx q[3];
rz(-2.0483096) q[3];
sx q[3];
rz(1.8846412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8811983) q[2];
sx q[2];
rz(-0.10513267) q[2];
sx q[2];
rz(-2.6283188) q[2];
rz(1.5848292) q[3];
sx q[3];
rz(-1.1791752) q[3];
sx q[3];
rz(1.5796327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8455115) q[0];
sx q[0];
rz(-0.054231461) q[0];
sx q[0];
rz(-0.84501141) q[0];
rz(-0.71146479) q[1];
sx q[1];
rz(-0.80159801) q[1];
sx q[1];
rz(0.62666384) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.959528) q[0];
sx q[0];
rz(-1.7668889) q[0];
sx q[0];
rz(-0.0053868731) q[0];
x q[1];
rz(0.22876744) q[2];
sx q[2];
rz(-1.1950777) q[2];
sx q[2];
rz(1.2496333) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9756445) q[1];
sx q[1];
rz(-2.2774146) q[1];
sx q[1];
rz(-2.9099275) q[1];
x q[2];
rz(1.5845039) q[3];
sx q[3];
rz(-2.6506069) q[3];
sx q[3];
rz(1.0233228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6381548) q[2];
sx q[2];
rz(-1.902709) q[2];
sx q[2];
rz(0.16866355) q[2];
rz(-0.1564129) q[3];
sx q[3];
rz(-0.63786879) q[3];
sx q[3];
rz(-2.8580581) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4446568) q[0];
sx q[0];
rz(-1.1035656) q[0];
sx q[0];
rz(-0.37100938) q[0];
rz(-1.3320529) q[1];
sx q[1];
rz(-1.5468372) q[1];
sx q[1];
rz(0.36088774) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9661117) q[0];
sx q[0];
rz(-1.4752441) q[0];
sx q[0];
rz(0.086009228) q[0];
x q[1];
rz(-0.12973327) q[2];
sx q[2];
rz(-3.0184944) q[2];
sx q[2];
rz(2.77746) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2672878) q[1];
sx q[1];
rz(-2.1308769) q[1];
sx q[1];
rz(-0.31927424) q[1];
x q[2];
rz(0.64334433) q[3];
sx q[3];
rz(-0.78668252) q[3];
sx q[3];
rz(-2.1718882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0767625) q[2];
sx q[2];
rz(-0.37134376) q[2];
sx q[2];
rz(-1.4999464) q[2];
rz(-2.0867945) q[3];
sx q[3];
rz(-2.0513556) q[3];
sx q[3];
rz(1.9487901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37831369) q[0];
sx q[0];
rz(-1.9739456) q[0];
sx q[0];
rz(1.8462697) q[0];
rz(0.12652215) q[1];
sx q[1];
rz(-0.69252068) q[1];
sx q[1];
rz(-1.2476547) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037654951) q[0];
sx q[0];
rz(-1.7154365) q[0];
sx q[0];
rz(-2.8580998) q[0];
x q[1];
rz(-1.1791897) q[2];
sx q[2];
rz(-0.84184781) q[2];
sx q[2];
rz(0.1192418) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3646255) q[1];
sx q[1];
rz(-2.7086621) q[1];
sx q[1];
rz(2.6294699) q[1];
x q[2];
rz(1.2887824) q[3];
sx q[3];
rz(-1.3913034) q[3];
sx q[3];
rz(0.13403758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.36914545) q[2];
sx q[2];
rz(-0.36448604) q[2];
sx q[2];
rz(0.44949731) q[2];
rz(-0.43826023) q[3];
sx q[3];
rz(-2.0354383) q[3];
sx q[3];
rz(1.0615758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19675572) q[0];
sx q[0];
rz(-1.0642835) q[0];
sx q[0];
rz(0.34360316) q[0];
rz(-0.65394872) q[1];
sx q[1];
rz(-2.3908354) q[1];
sx q[1];
rz(-0.59069815) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78498392) q[0];
sx q[0];
rz(-1.3425047) q[0];
sx q[0];
rz(-2.0697631) q[0];
rz(-pi) q[1];
rz(-2.3542064) q[2];
sx q[2];
rz(-0.91286589) q[2];
sx q[2];
rz(-1.3546561) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0456344) q[1];
sx q[1];
rz(-1.5812629) q[1];
sx q[1];
rz(-1.4399308) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77883522) q[3];
sx q[3];
rz(-0.86495728) q[3];
sx q[3];
rz(-0.51222982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.58064738) q[2];
sx q[2];
rz(-2.2717768) q[2];
sx q[2];
rz(-2.297164) q[2];
rz(-1.0405509) q[3];
sx q[3];
rz(-1.8620164) q[3];
sx q[3];
rz(-1.1738663) q[3];
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
rz(pi/2) q[3];
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
rz(1.3867253) q[0];
sx q[0];
rz(-2.6193021) q[0];
sx q[0];
rz(2.0544384) q[0];
rz(-0.48671752) q[1];
sx q[1];
rz(-1.646128) q[1];
sx q[1];
rz(1.4901935) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8069902) q[0];
sx q[0];
rz(-2.2740472) q[0];
sx q[0];
rz(-3.0038805) q[0];
rz(-pi) q[1];
rz(-1.0045687) q[2];
sx q[2];
rz(-0.90551584) q[2];
sx q[2];
rz(0.69349223) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9864038) q[1];
sx q[1];
rz(-1.1935368) q[1];
sx q[1];
rz(-0.44694408) q[1];
x q[2];
rz(2.9341615) q[3];
sx q[3];
rz(-1.279502) q[3];
sx q[3];
rz(2.8136611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3460441) q[2];
sx q[2];
rz(-2.2717924) q[2];
sx q[2];
rz(-2.8022433) q[2];
rz(-0.36335534) q[3];
sx q[3];
rz(-0.27000579) q[3];
sx q[3];
rz(1.5615162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.617008) q[0];
sx q[0];
rz(-1.120485) q[0];
sx q[0];
rz(-2.8560915) q[0];
rz(0.78954804) q[1];
sx q[1];
rz(-1.6098846) q[1];
sx q[1];
rz(1.550536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68398482) q[0];
sx q[0];
rz(-1.8681861) q[0];
sx q[0];
rz(-1.7575197) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3197081) q[2];
sx q[2];
rz(-0.54080117) q[2];
sx q[2];
rz(3.1201517) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62666368) q[1];
sx q[1];
rz(-2.0928934) q[1];
sx q[1];
rz(1.0842647) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66724369) q[3];
sx q[3];
rz(-1.6341471) q[3];
sx q[3];
rz(-0.9982619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8153814) q[2];
sx q[2];
rz(-2.9823494) q[2];
sx q[2];
rz(-0.86686575) q[2];
rz(2.4032812) q[3];
sx q[3];
rz(-1.5574995) q[3];
sx q[3];
rz(-2.2332634) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0758122) q[0];
sx q[0];
rz(-0.74956885) q[0];
sx q[0];
rz(2.442389) q[0];
rz(2.5993644) q[1];
sx q[1];
rz(-1.7991853) q[1];
sx q[1];
rz(0.29621616) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2767773) q[0];
sx q[0];
rz(-2.3010198) q[0];
sx q[0];
rz(-0.093412799) q[0];
rz(-pi) q[1];
rz(2.2454295) q[2];
sx q[2];
rz(-0.8783717) q[2];
sx q[2];
rz(-1.2330173) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51402873) q[1];
sx q[1];
rz(-0.96935287) q[1];
sx q[1];
rz(2.264643) q[1];
rz(1.7858539) q[3];
sx q[3];
rz(-0.38588215) q[3];
sx q[3];
rz(2.1467211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2518623) q[2];
sx q[2];
rz(-2.6984062) q[2];
sx q[2];
rz(2.6207793) q[2];
rz(-3.1098747) q[3];
sx q[3];
rz(-2.7364276) q[3];
sx q[3];
rz(-0.051137663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.48145914) q[0];
sx q[0];
rz(-1.0986468) q[0];
sx q[0];
rz(0.81137401) q[0];
rz(2.7157281) q[1];
sx q[1];
rz(-0.90161294) q[1];
sx q[1];
rz(-1.3070377) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9872914) q[0];
sx q[0];
rz(-1.5736921) q[0];
sx q[0];
rz(-1.5737246) q[0];
rz(-pi) q[1];
rz(-1.6128236) q[2];
sx q[2];
rz(-2.6767455) q[2];
sx q[2];
rz(0.56715172) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.20504925) q[1];
sx q[1];
rz(-1.8920415) q[1];
sx q[1];
rz(0.46830362) q[1];
x q[2];
rz(-1.3974992) q[3];
sx q[3];
rz(-1.2504645) q[3];
sx q[3];
rz(-0.082435247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9797152) q[2];
sx q[2];
rz(-0.76354176) q[2];
sx q[2];
rz(2.7757576) q[2];
rz(-0.93519917) q[3];
sx q[3];
rz(-1.5949944) q[3];
sx q[3];
rz(2.759815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.81015051) q[0];
sx q[0];
rz(-1.3642385) q[0];
sx q[0];
rz(1.4862899) q[0];
rz(1.4797795) q[1];
sx q[1];
rz(-1.231709) q[1];
sx q[1];
rz(-0.92373893) q[1];
rz(0.059546373) q[2];
sx q[2];
rz(-2.1179143) q[2];
sx q[2];
rz(2.2512022) q[2];
rz(-1.4167841) q[3];
sx q[3];
rz(-1.3091675) q[3];
sx q[3];
rz(-3.0021734) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
