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
rz(0.44952196) q[0];
sx q[0];
rz(4.8604453) q[0];
sx q[0];
rz(7.3331375) q[0];
rz(2.6226251) q[1];
sx q[1];
rz(-1.608404) q[1];
sx q[1];
rz(3.0906711) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.19456) q[0];
sx q[0];
rz(-1.2480535) q[0];
sx q[0];
rz(2.6850924) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21587431) q[2];
sx q[2];
rz(-1.8980489) q[2];
sx q[2];
rz(0.33282166) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5320961) q[1];
sx q[1];
rz(-1.8724964) q[1];
sx q[1];
rz(-2.8068551) q[1];
rz(-pi) q[2];
rz(-1.2818358) q[3];
sx q[3];
rz(-0.70631344) q[3];
sx q[3];
rz(0.044580288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8672987) q[2];
sx q[2];
rz(-0.87624246) q[2];
sx q[2];
rz(1.743861) q[2];
rz(-2.6869669) q[3];
sx q[3];
rz(-2.2635098) q[3];
sx q[3];
rz(1.4286058) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4947263) q[0];
sx q[0];
rz(-0.88748256) q[0];
sx q[0];
rz(-2.5383762) q[0];
rz(-1.2546722) q[1];
sx q[1];
rz(-0.61379495) q[1];
sx q[1];
rz(-2.5616554) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3361084) q[0];
sx q[0];
rz(-1.0382129) q[0];
sx q[0];
rz(-1.3773772) q[0];
rz(1.1179148) q[2];
sx q[2];
rz(-2.09399) q[2];
sx q[2];
rz(0.69566899) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7407826) q[1];
sx q[1];
rz(-2.1745958) q[1];
sx q[1];
rz(-3.0933988) q[1];
rz(-pi) q[2];
rz(-0.55508695) q[3];
sx q[3];
rz(-2.4984341) q[3];
sx q[3];
rz(-0.023338524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4846399) q[2];
sx q[2];
rz(-2.3886949) q[2];
sx q[2];
rz(2.7638655) q[2];
rz(1.7083302) q[3];
sx q[3];
rz(-1.6323615) q[3];
sx q[3];
rz(-1.1908971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.3483873) q[0];
sx q[0];
rz(-0.054940104) q[0];
sx q[0];
rz(-1.4980263) q[0];
rz(1.161423) q[1];
sx q[1];
rz(-1.8893416) q[1];
sx q[1];
rz(0.84322554) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0204165) q[0];
sx q[0];
rz(-1.5112229) q[0];
sx q[0];
rz(-0.94734933) q[0];
rz(-pi) q[1];
rz(1.9887504) q[2];
sx q[2];
rz(-1.770263) q[2];
sx q[2];
rz(-2.0355088) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9084042) q[1];
sx q[1];
rz(-0.68906765) q[1];
sx q[1];
rz(1.3702964) q[1];
rz(-pi) q[2];
rz(-0.70103443) q[3];
sx q[3];
rz(-1.8314519) q[3];
sx q[3];
rz(3.1150027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8680385) q[2];
sx q[2];
rz(-0.66323438) q[2];
sx q[2];
rz(-0.794945) q[2];
rz(-2.3718209) q[3];
sx q[3];
rz(-0.92307463) q[3];
sx q[3];
rz(2.9845089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5539826) q[0];
sx q[0];
rz(-1.8164604) q[0];
sx q[0];
rz(-0.28833589) q[0];
rz(0.68710697) q[1];
sx q[1];
rz(-1.4930875) q[1];
sx q[1];
rz(1.4289325) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16609678) q[0];
sx q[0];
rz(-1.5642484) q[0];
sx q[0];
rz(3.1295144) q[0];
x q[1];
rz(1.7111139) q[2];
sx q[2];
rz(-1.9694917) q[2];
sx q[2];
rz(1.2600419) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7894788) q[1];
sx q[1];
rz(-1.6360456) q[1];
sx q[1];
rz(-0.25635527) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5023584) q[3];
sx q[3];
rz(-1.0022638) q[3];
sx q[3];
rz(2.6871339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5700506) q[2];
sx q[2];
rz(-2.1752581) q[2];
sx q[2];
rz(-0.16560444) q[2];
rz(-0.064519493) q[3];
sx q[3];
rz(-0.12972984) q[3];
sx q[3];
rz(1.8701514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22873838) q[0];
sx q[0];
rz(-1.9485291) q[0];
sx q[0];
rz(-2.2304529) q[0];
rz(0.97081026) q[1];
sx q[1];
rz(-1.3263005) q[1];
sx q[1];
rz(0.86311805) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1939094) q[0];
sx q[0];
rz(-1.860431) q[0];
sx q[0];
rz(-1.0632443) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8392856) q[2];
sx q[2];
rz(-1.5483861) q[2];
sx q[2];
rz(0.42699285) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6533901) q[1];
sx q[1];
rz(-1.0428671) q[1];
sx q[1];
rz(-1.3009682) q[1];
x q[2];
rz(1.8831598) q[3];
sx q[3];
rz(-1.527597) q[3];
sx q[3];
rz(0.13052065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0866278) q[2];
sx q[2];
rz(-1.8885771) q[2];
sx q[2];
rz(-2.6524554) q[2];
rz(-0.9984115) q[3];
sx q[3];
rz(-0.17397927) q[3];
sx q[3];
rz(-0.099099549) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2973706) q[0];
sx q[0];
rz(-0.64592823) q[0];
sx q[0];
rz(-0.55322629) q[0];
rz(0.79752254) q[1];
sx q[1];
rz(-0.99634606) q[1];
sx q[1];
rz(1.1513938) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11565514) q[0];
sx q[0];
rz(-2.135072) q[0];
sx q[0];
rz(-1.7274117) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55091605) q[2];
sx q[2];
rz(-1.7168432) q[2];
sx q[2];
rz(-2.3178326) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7094671) q[1];
sx q[1];
rz(-1.8399939) q[1];
sx q[1];
rz(1.7749191) q[1];
rz(0.90535935) q[3];
sx q[3];
rz(-1.8784154) q[3];
sx q[3];
rz(0.24505982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3449751) q[2];
sx q[2];
rz(-1.4130219) q[2];
sx q[2];
rz(2.3100992) q[2];
rz(-1.8031395) q[3];
sx q[3];
rz(-2.0667388) q[3];
sx q[3];
rz(-0.89053806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3727342) q[0];
sx q[0];
rz(-1.9356198) q[0];
sx q[0];
rz(0.15705577) q[0];
rz(-1.318469) q[1];
sx q[1];
rz(-2.1993957) q[1];
sx q[1];
rz(2.7838321) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46808711) q[0];
sx q[0];
rz(-2.0131677) q[0];
sx q[0];
rz(2.6297188) q[0];
rz(-2.116647) q[2];
sx q[2];
rz(-2.3966925) q[2];
sx q[2];
rz(-2.2226983) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.50510943) q[1];
sx q[1];
rz(-0.92527522) q[1];
sx q[1];
rz(-2.4277975) q[1];
rz(-0.24202204) q[3];
sx q[3];
rz(-0.58804711) q[3];
sx q[3];
rz(-0.96554276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1498108) q[2];
sx q[2];
rz(-1.8842183) q[2];
sx q[2];
rz(-0.020616654) q[2];
rz(-1.2425544) q[3];
sx q[3];
rz(-0.81762448) q[3];
sx q[3];
rz(-2.8500565) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5307411) q[0];
sx q[0];
rz(-1.1854956) q[0];
sx q[0];
rz(0.32671842) q[0];
rz(0.84699455) q[1];
sx q[1];
rz(-1.0495443) q[1];
sx q[1];
rz(-2.0832031) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0857226) q[0];
sx q[0];
rz(-0.51858178) q[0];
sx q[0];
rz(-1.4297755) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63916529) q[2];
sx q[2];
rz(-1.4286388) q[2];
sx q[2];
rz(-0.41691142) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10524532) q[1];
sx q[1];
rz(-0.39064841) q[1];
sx q[1];
rz(1.7772783) q[1];
rz(-pi) q[2];
x q[2];
rz(2.136738) q[3];
sx q[3];
rz(-1.3019239) q[3];
sx q[3];
rz(0.0026800935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1191001) q[2];
sx q[2];
rz(-0.34467003) q[2];
sx q[2];
rz(1.1723088) q[2];
rz(-3.1398224) q[3];
sx q[3];
rz(-0.5032731) q[3];
sx q[3];
rz(-2.1625904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86972791) q[0];
sx q[0];
rz(-2.3379022) q[0];
sx q[0];
rz(-1.3386238) q[0];
rz(1.9610693) q[1];
sx q[1];
rz(-1.0931284) q[1];
sx q[1];
rz(-2.6611633) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.521489) q[0];
sx q[0];
rz(-2.0248981) q[0];
sx q[0];
rz(0.2257077) q[0];
x q[1];
rz(0.24238853) q[2];
sx q[2];
rz(-0.17504642) q[2];
sx q[2];
rz(2.4180061) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.980396) q[1];
sx q[1];
rz(-0.56884407) q[1];
sx q[1];
rz(-1.3854762) q[1];
rz(-pi) q[2];
rz(-0.35828405) q[3];
sx q[3];
rz(-1.7833774) q[3];
sx q[3];
rz(2.6019118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7052445) q[2];
sx q[2];
rz(-0.99506012) q[2];
sx q[2];
rz(-2.056541) q[2];
rz(-1.1121701) q[3];
sx q[3];
rz(-0.89434904) q[3];
sx q[3];
rz(-0.81356847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28806624) q[0];
sx q[0];
rz(-2.1948094) q[0];
sx q[0];
rz(2.1296401) q[0];
rz(2.2400253) q[1];
sx q[1];
rz(-0.26600599) q[1];
sx q[1];
rz(0.56545767) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31663592) q[0];
sx q[0];
rz(-0.68616435) q[0];
sx q[0];
rz(1.9030722) q[0];
rz(-pi) q[1];
rz(2.1109963) q[2];
sx q[2];
rz(-2.6517617) q[2];
sx q[2];
rz(-0.67661197) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8248837) q[1];
sx q[1];
rz(-1.663534) q[1];
sx q[1];
rz(0.85272654) q[1];
rz(-0.10528414) q[3];
sx q[3];
rz(-1.8859047) q[3];
sx q[3];
rz(1.9721102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3507583) q[2];
sx q[2];
rz(-1.5459205) q[2];
sx q[2];
rz(1.2606384) q[2];
rz(0.44220051) q[3];
sx q[3];
rz(-1.0298157) q[3];
sx q[3];
rz(-1.376576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.5762536) q[0];
sx q[0];
rz(-0.88247846) q[0];
sx q[0];
rz(-0.89378617) q[0];
rz(-1.5743938) q[1];
sx q[1];
rz(-1.6849453) q[1];
sx q[1];
rz(-1.4243855) q[1];
rz(0.5657351) q[2];
sx q[2];
rz(-2.6489352) q[2];
sx q[2];
rz(0.87633662) q[2];
rz(-0.16408739) q[3];
sx q[3];
rz(-1.224095) q[3];
sx q[3];
rz(-1.7188354) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
