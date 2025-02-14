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
rz(0.32938862) q[0];
sx q[0];
rz(3.7731054) q[0];
sx q[0];
rz(9.4256529) q[0];
rz(2.5007091) q[1];
sx q[1];
rz(-2.1407318) q[1];
sx q[1];
rz(-0.34520087) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76287718) q[0];
sx q[0];
rz(-1.595531) q[0];
sx q[0];
rz(1.4799165) q[0];
rz(-2.7617747) q[2];
sx q[2];
rz(-1.4336515) q[2];
sx q[2];
rz(-1.6452546) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7789014) q[1];
sx q[1];
rz(-0.76299113) q[1];
sx q[1];
rz(-2.1982026) q[1];
rz(-1.4963989) q[3];
sx q[3];
rz(-1.5545903) q[3];
sx q[3];
rz(0.61283535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9125646) q[2];
sx q[2];
rz(-1.3338858) q[2];
sx q[2];
rz(2.933617) q[2];
rz(0.29911706) q[3];
sx q[3];
rz(-2.5695473) q[3];
sx q[3];
rz(-2.0007029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3308554) q[0];
sx q[0];
rz(-0.24092291) q[0];
sx q[0];
rz(-2.148707) q[0];
rz(-1.317124) q[1];
sx q[1];
rz(-0.31610745) q[1];
sx q[1];
rz(-2.3670926) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68607722) q[0];
sx q[0];
rz(-0.17855274) q[0];
sx q[0];
rz(1.6327842) q[0];
rz(-2.108091) q[2];
sx q[2];
rz(-2.3024493) q[2];
sx q[2];
rz(3.0947859) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.61770536) q[1];
sx q[1];
rz(-1.7951709) q[1];
sx q[1];
rz(2.7588506) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22529545) q[3];
sx q[3];
rz(-2.0301798) q[3];
sx q[3];
rz(-2.6090906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2449067) q[2];
sx q[2];
rz(-1.2865571) q[2];
sx q[2];
rz(-1.0150821) q[2];
rz(2.8924938) q[3];
sx q[3];
rz(-2.2738012) q[3];
sx q[3];
rz(0.36756137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-2.2772813) q[0];
sx q[0];
rz(-0.69121498) q[0];
sx q[0];
rz(-0.15750289) q[0];
rz(-1.0109673) q[1];
sx q[1];
rz(-0.33282655) q[1];
sx q[1];
rz(0.13883042) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.39108) q[0];
sx q[0];
rz(-1.8939928) q[0];
sx q[0];
rz(1.9531519) q[0];
rz(-1.1558258) q[2];
sx q[2];
rz(-0.89902821) q[2];
sx q[2];
rz(2.6443554) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.638354) q[1];
sx q[1];
rz(-1.279083) q[1];
sx q[1];
rz(0.20511638) q[1];
rz(-pi) q[2];
rz(2.5957727) q[3];
sx q[3];
rz(-0.32836093) q[3];
sx q[3];
rz(-2.2203746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6802754) q[2];
sx q[2];
rz(-2.2291144) q[2];
sx q[2];
rz(-0.3581363) q[2];
rz(2.4628468) q[3];
sx q[3];
rz(-0.94921422) q[3];
sx q[3];
rz(2.1024735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045227483) q[0];
sx q[0];
rz(-0.97310936) q[0];
sx q[0];
rz(2.8367693) q[0];
rz(-2.7335956) q[1];
sx q[1];
rz(-1.7034737) q[1];
sx q[1];
rz(-0.92794424) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5281677) q[0];
sx q[0];
rz(-1.6990464) q[0];
sx q[0];
rz(-1.7665461) q[0];
rz(-pi) q[1];
rz(-2.5408542) q[2];
sx q[2];
rz(-1.033412) q[2];
sx q[2];
rz(1.635765) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.4907704) q[1];
sx q[1];
rz(-2.4374552) q[1];
sx q[1];
rz(2.9068374) q[1];
rz(0.66189142) q[3];
sx q[3];
rz(-0.74569476) q[3];
sx q[3];
rz(2.6357366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1912332) q[2];
sx q[2];
rz(-2.5706036) q[2];
sx q[2];
rz(2.9534269) q[2];
rz(1.865271) q[3];
sx q[3];
rz(-1.3151582) q[3];
sx q[3];
rz(1.7108542) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6998049) q[0];
sx q[0];
rz(-2.7975174) q[0];
sx q[0];
rz(-3.1222043) q[0];
rz(-2.0806606) q[1];
sx q[1];
rz(-1.414199) q[1];
sx q[1];
rz(1.2394989) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083154924) q[0];
sx q[0];
rz(-2.4920336) q[0];
sx q[0];
rz(2.044528) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73915787) q[2];
sx q[2];
rz(-2.0817167) q[2];
sx q[2];
rz(2.7834653) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2693138) q[1];
sx q[1];
rz(-1.2904823) q[1];
sx q[1];
rz(-2.6107236) q[1];
x q[2];
rz(-1.0723713) q[3];
sx q[3];
rz(-1.2066505) q[3];
sx q[3];
rz(2.2143603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1218607) q[2];
sx q[2];
rz(-1.3132361) q[2];
sx q[2];
rz(1.3266374) q[2];
rz(2.895368) q[3];
sx q[3];
rz(-1.1028057) q[3];
sx q[3];
rz(-0.8518014) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5354079) q[0];
sx q[0];
rz(-2.1562205) q[0];
sx q[0];
rz(2.4024409) q[0];
rz(2.7958561) q[1];
sx q[1];
rz(-1.812499) q[1];
sx q[1];
rz(-0.87127042) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3987253) q[0];
sx q[0];
rz(-2.350432) q[0];
sx q[0];
rz(-1.4732811) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61001444) q[2];
sx q[2];
rz(-0.84976174) q[2];
sx q[2];
rz(-0.80243669) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4193486) q[1];
sx q[1];
rz(-2.7461395) q[1];
sx q[1];
rz(1.761318) q[1];
rz(-pi) q[2];
rz(0.61136758) q[3];
sx q[3];
rz(-1.1725559) q[3];
sx q[3];
rz(-2.7643725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8295916) q[2];
sx q[2];
rz(-2.5025554) q[2];
sx q[2];
rz(-2.269022) q[2];
rz(-1.5729337) q[3];
sx q[3];
rz(-0.60397732) q[3];
sx q[3];
rz(-0.93305552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(2.0912112) q[0];
sx q[0];
rz(-0.25930431) q[0];
sx q[0];
rz(0.76164371) q[0];
rz(-0.42539445) q[1];
sx q[1];
rz(-2.0014706) q[1];
sx q[1];
rz(0.92686191) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0480014) q[0];
sx q[0];
rz(-2.8107852) q[0];
sx q[0];
rz(0.76865102) q[0];
rz(-0.71304597) q[2];
sx q[2];
rz(-2.2879061) q[2];
sx q[2];
rz(-0.46592679) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.35173479) q[1];
sx q[1];
rz(-2.7119293) q[1];
sx q[1];
rz(-1.0142782) q[1];
x q[2];
rz(0.048153444) q[3];
sx q[3];
rz(-2.1934436) q[3];
sx q[3];
rz(0.66648713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0285792) q[2];
sx q[2];
rz(-2.4455652) q[2];
sx q[2];
rz(2.2704303) q[2];
rz(2.8431559) q[3];
sx q[3];
rz(-1.8961743) q[3];
sx q[3];
rz(-1.8106073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7262064) q[0];
sx q[0];
rz(-2.6415249) q[0];
sx q[0];
rz(-0.4183847) q[0];
rz(-2.8437974) q[1];
sx q[1];
rz(-1.9458385) q[1];
sx q[1];
rz(-0.13430886) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.031536438) q[0];
sx q[0];
rz(-0.93288619) q[0];
sx q[0];
rz(-0.67968926) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9842582) q[2];
sx q[2];
rz(-0.43817876) q[2];
sx q[2];
rz(-2.7587121) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.93523504) q[1];
sx q[1];
rz(-0.22769732) q[1];
sx q[1];
rz(0.63657324) q[1];
rz(-1.1014492) q[3];
sx q[3];
rz(-1.3162287) q[3];
sx q[3];
rz(1.001844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7316651) q[2];
sx q[2];
rz(-1.8525367) q[2];
sx q[2];
rz(0.64201391) q[2];
rz(-1.0632473) q[3];
sx q[3];
rz(-0.65170538) q[3];
sx q[3];
rz(2.1003907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6006271) q[0];
sx q[0];
rz(-3.0525115) q[0];
sx q[0];
rz(-3.0294321) q[0];
rz(1.3345831) q[1];
sx q[1];
rz(-2.5490675) q[1];
sx q[1];
rz(1.0293915) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4698668) q[0];
sx q[0];
rz(-0.70505667) q[0];
sx q[0];
rz(1.5211578) q[0];
rz(-pi) q[1];
rz(3.0156187) q[2];
sx q[2];
rz(-2.8354037) q[2];
sx q[2];
rz(0.15737113) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2160714) q[1];
sx q[1];
rz(-0.23434429) q[1];
sx q[1];
rz(1.7611124) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46799) q[3];
sx q[3];
rz(-1.6525998) q[3];
sx q[3];
rz(1.3036659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7182497) q[2];
sx q[2];
rz(-3.0666879) q[2];
sx q[2];
rz(-3.1035799) q[2];
rz(-0.612261) q[3];
sx q[3];
rz(-0.93658787) q[3];
sx q[3];
rz(3.0003701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8925979) q[0];
sx q[0];
rz(-1.4709512) q[0];
sx q[0];
rz(-0.41879642) q[0];
rz(1.8839802) q[1];
sx q[1];
rz(-1.6323171) q[1];
sx q[1];
rz(0.14258252) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0973546) q[0];
sx q[0];
rz(-1.729106) q[0];
sx q[0];
rz(1.5990785) q[0];
rz(0.37303961) q[2];
sx q[2];
rz(-2.4786886) q[2];
sx q[2];
rz(-2.474624) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.31227641) q[1];
sx q[1];
rz(-1.8354776) q[1];
sx q[1];
rz(-0.97977248) q[1];
rz(-pi) q[2];
rz(1.110465) q[3];
sx q[3];
rz(-1.8813731) q[3];
sx q[3];
rz(-0.47720695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3451781) q[2];
sx q[2];
rz(-3.0609481) q[2];
sx q[2];
rz(-2.8734015) q[2];
rz(-1.7372519) q[3];
sx q[3];
rz(-1.0920478) q[3];
sx q[3];
rz(1.0442737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0161229) q[0];
sx q[0];
rz(-2.7653427) q[0];
sx q[0];
rz(0.34633037) q[0];
rz(-0.12915962) q[1];
sx q[1];
rz(-1.8864514) q[1];
sx q[1];
rz(1.4056978) q[1];
rz(-0.47507349) q[2];
sx q[2];
rz(-1.8865841) q[2];
sx q[2];
rz(0.24812698) q[2];
rz(1.3553452) q[3];
sx q[3];
rz(-0.30134311) q[3];
sx q[3];
rz(-2.338196) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
